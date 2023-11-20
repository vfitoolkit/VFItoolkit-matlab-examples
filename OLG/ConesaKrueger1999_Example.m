% Example based on Conesa & Krueger (1999) - Social Security Reform with Heterogeneous Agents

% Solves an OLG model for a general equilibrium transition path

transpathoptions.fastOLG=0;
% fastOLG=1 means that the code parallelizes over the finite-horizon (J)
% while computing the transition path. It is faster but requires more memory.
% fastOLG=0 will be slower but use less memory

% You can choose from three options for the exogenous shock process
IndivProdShock=1; % 0 for 'no heterogeneity', 1 for Diaz-Gimenez et al (1997) the 'asymmetric case', 2 for Storesletten et al (1998) the 'symmetric case'.
% You can choose from three options for what reform to consider
PolicyReform=1; % 1=Policy Reform A,  2=Policy Reform B,  3=Policy Reform C

% Note: CK1999 model can actually be solved without decision variable. (You can calculate analytic value for labor supply based on a and aprime, see
% CK1999 paper. I model it explicitly so that code is easier to modify for other uses.)
n_d=51; % fraction of time worked
n_a=301; % assets
% n_z % labour efficiency units, depends on 'IndivProdShock', so set below.
N_j=66; % age (number of periods, to be precise this is age-19)

% A few lines I needed for running on the Server
addpath(genpath('./MatlabToolkits/'))

% Set simoptions.ncores
simoptions.ncores=feature('numcores'); % Number of CPU cores

vfoptions.policy_forceintegertype=2; % Policy was not being treated as integers (one of the elements was 10^(-15) different from an integer)
heteroagentoptions.verbose=1;

%% Set model parameters

%Preference parameters
Params.sigma=2;
Params.beta=0.97;
Params.gamma=0.42;

%Technology parameters
Params.alpha=0.36;
Params.delta=0.06;
Params.theta=1;

%People become economically active (j=1) at age 20, retire at 65, and max age is 85
Params.J=N_j; %=85-19
Params.Jr=46; % =65-19

%Population growth of 1.1%
Params.n=0.011;
%Social Security replacement rate
Params.b_initial=0.5; Params.b=Params.b_initial;

% This commented out part is how I has originally did the survival
% probabilities. I later received the codes of Conesa & Krueger by email from them and these are copy-pasted below.
% Probability of surviving, conditional on being alive (Conesa & Krueger 1999 just say that they get them from Faber, 1982)
% The commented out lines below contain my original version using more recent numbers from 'same' source.
% %As this document is not available online I instead take the values from the updated version of the document, namely from 
% %F. C. Bell and M. L. Miller (2005), Life Tables for the United States Social Security Area 1900-2100, Actuarial Study No. 120, Office of the Chief Actuary
% %http://www.socialsecurity.gov/oact/NOTES/s2000s.html
% %Table 8 — Period Probabilities of Death Within One Year (qx) at Selected Exact Ages, by Sex and Calendar Year (Cont.)
% %The raw data from there is
% %          Sex and Exact Age
% %    |  Male                                                Female
% %Year| [0 30 60 65 70 100]                                  [0 30 60 65 70 100]
% %2010| [0.00587,0.00116,0.01086,0.01753,0.02785,0.39134]    [0.00495,0.00060,0.00734,0.01201,0.01912,0.34031]
% %I just take the numbers for Males, and then set my actual values based on a linear interpolation of the data.
% dj_temp=interp1([0,30,60,65,70,100],[0.00587,0.00116,0.01086,0.01753,0.02785,0.39134],0:1:100,'linear');
% Params.sj=1-dj_temp(20:85);
% Params.sj(1)=1;
% Params.sj(end)=0;
% % I use linear interpolation to fill in the numbers inbetween those reported by Bell & Miller (2005).
% % I have aditionally imposed that the prob of death at age 20 be zero and that prob of death at age 85 is one.

% Following is taken from original codes that were sent to me by email by Conesa & Krueger (demographics.f90):
% Population Numbers from Faber (1982) (this line was their comment)
PopulationData=[196729, 196546, 196357, 196161, 195958, 195747, 195530, 195306, 195082, 194864, 194655, 194457, 194266, 194080, 193893, 193699, 193497, 193283, 193054, 192806,...
    192533, 192232, 191897, 191527, 191120, 190671, 190177, 189633, 189031, 188358, 187610, 186774, 185846, 184836, 183752, 182599, 181373, 180058, 178634, 177077,...
    175367, 173494, 171452, 169223, 166789, 164138, 161257, 158158, 154847, 151369, 147740, 143958, 140011, 135897, 131617, 127168, 122543, 117741, 112772, 107648,...
	102379, 96966, 91412, 85731, 79939, 74066, 68144, 62209, 56302, 50467, 44755, 39220, 33920, 28915, 24259, 20005, 16209, 12908, 10108, 7792, 5922];
% Survival probabilities: phi(i)=prob(alive in i+1|alive in i) (their comment; phi(i) is what I here call sj)
Params.sj=PopulationData(2:end)./PopulationData(1:end-1);
% Note, this is more ages than we need/want. So just grab those which are relevant
Params.sj=Params.sj(1:Params.J);
Params.sj(end)=0; % Fill in a zero for final period (note that the value of zero here is actually irrelevant to codes, I do it purely as CK1999 do in their codes: demographics.f90)

% Age profile of productivity (based on lifetime profile of earnings, Hansen 1993), Epsilon_j
% I have to interpolate this in some manner to get the values for each age. I am not sure how CK1999 did the interpolation as they never mention it.
% Again this is same source as Imrohoroglu, Imrohoroglu & Joines (1995) and
% while replicating that paper they sent me a copy of the their interpolated numbers which I will assume are same as CK1999.
%First the raw data of Hansen 1993:
%Table II. Weights assigned to age-sex
%         Males       Females
%Age      Weight      Weight
%16-19    0.56        0.52
%20-24    0.78        0.69
%25-34    1.14        0.89
%35-44    1.37        0.90
%45-54    1.39        0.87
%55-64    1.33        0.84
%65 +     0.89        0.66
%14-17    0.56        0.52
%14-19    0.56        0.52
%18-24    0.78        0.69
%25-44    1.24        0.89
%45-64    1.37        0.86
%epsilon_j is set based on the first entries of the data for males
% The following commented out part is my original attempt
% epsilon_j=zeros(Params.J,1); 
% epsilon_j(1:5)=0.78; epsilon_j(6:15)=1.14; epsilon_j(16:25)=1.37; 
% epsilon_j(26:35)=1.39; epsilon_j(36:45)=1.33; epsilon_j(46:66)=0.89;
% Params.epsilon_j=epsilon_j/epsilon_j(1); %Conesa & Krueger refer to labour efficiency as 1 at age 20, so presumably they make some renormalization like this
% Following lines are copy-pasted from the codes of Conesa & Krueger (which they sent me by email; demographics.f90)
Params.epsilon_j=zeros(Params.J,1); % CK1999 code fills in everything from period 46 on with zeros
Params.epsilon_j(1:45)=[1.0000,1.0719,1.1438, 1.2158, 1.2842, 1.3527, 1.4212, 1.4897, 1.5582, 1.6267, 1.6952, 1.7217, 1.7438, 1.7748, 1.8014, 1.8279, 1.8545, 1.8810, 1.9075, 1.9341,...
    1.9606, 1.9623, 1.9640, 1.9658, 1.9675, 1.9692, 1.9709, 1.9726, 1.9743, 1.9760, 1.9777, 1.9700, 1.9623, 1.9546, 1.9469, 1.9392, 1.9315, 1.9238, 1.9161, 1.9084, 1.9007, 1.8354, 1.7701, 1.7048, 1.6396]';
% Note that this already includes the normalization of labor efficiency to 1 at age 20
% Conesa & Krueger ban retirees from working. One way to do this is just to
% set epsilon_j=0 for retirees. Rather than do this via epsilon_j it is
% done by I_j which in practice is multiplied by epsilon_j.
Params.I_j=[ones(Params.Jr-1,1);zeros(Params.J-Params.Jr+1,1)];
% Note that with my original version I_j was needed, but is redundant using the original epsilon_j numbers from CK1999


vfoptions.verbose=0;
vfoptions.policy_forceintegertype=2; % Policy was not being treated as integers (one of the elements was 10^(-15) different from an integer)
vfoptions.lowmemory=0; % This is changed to other values for some of the model variants.
% simoptions.nsims=4*10^5;
% simoptions.iterate=1;
heteroagentoptions.verbose=1;


%% Set exogneous shock process based on IndivProdShock
% Parameters in Table IV: 1 for Storesletten et al (1998) the 'symmetric case', 2 for Diaz-Gimenez et al (1997) the 'asymmetric case'.
if IndivProdShock==0 % 0 for the 'no heterogeneity'
    % CK1999 never explicitly describe what 'no hetergeneity' means. But
    % given that idiosyncratic productivity shocks eta are the source of
    % (ex-post) heterogeneity presumably it is simply a matter of setting them to their mean value of 1.
    eta_grid=1;
    pi_eta=1;
elseif IndivProdShock==1 % 1 for Storesletten et al (1998) the 'symmetric case'.
    eta_grid=[0.73; 1.27];
    pi_eta=[0.82, 0.18; 0.18, 0.82];
    % Following lines are needed just to create Table 4
    Params.eta1=eta_grid(1); Params.eta2=eta_grid(2);
    Params.pi1=pi_eta(1,1); Params.pi2=pi_eta(2,2);
elseif IndivProdShock==2 % 2 for Diaz-Gimenez et al (1997) the 'asymmetric case'
    eta_grid=[0.5;3];
    pi_eta=[0.9811, 1-0.9811; 1-0.9261, 0.9261];
    % Following lines are needed just to create Table 4
    Params.eta1=eta_grid(1); Params.eta2=eta_grid(2);
    Params.pi1=pi_eta(1,1); Params.pi2=pi_eta(2,2);
end
n_z=length(eta_grid);
z_grid=eta_grid;
pi_z=pi_eta;
N_z=prod(n_z);

%% Grid for assets
k_grid=55*(linspace(0,1,n_a).^3)';

%% Get problem into format for using the toolkit
d_grid=linspace(0,1,n_d)'; %fraction of time worked
a_grid=k_grid; %assets

%% Calculate the population distribution across the ages (population at time
% t is made stationary by dividing it by \Pi_{i=1}^{t} (1+n_{i}) (product of all
% growth rates up till current period))
% CK2009 do not appear to give it an name (they instead say in
% Computational Appendix that they simply work with the model after making
% all the adjustments for population growth). The VFI Toolkit needs to give
% it a name so that it can be automatically used when calculating model
% outputs.
Params.mewj=ones(Params.J,1);
for jj=2:Params.J
    Params.mewj(jj)=Params.mewj(jj-1)*(1/(1+Params.n))*Params.sj(jj-1);
end
Params.mewj=Params.mewj/sum(Params.mewj); % normalize to measure one
Params.mewj=Params.mewj'; % Age weights must be a row vector.

AgeWeightsParamNames={'mewj'}; % Many finite horizon models apply different weights to different 'ages'; eg., due to survival rates or population growth rates.

% This population distribution across the ages is also used to back out initial guess for the calibrated value for tau.
Params.tau_initial=Params.b*sum(Params.mewj(:,Params.Jr:Params.J))/sum(Params.mewj(:,1:Params.Jr-1)); % From eqn 2.13 and 2.12 we get: tau=b*(frac of population retired)/(fraction of population of working age)
% Note that rather than use this directly as the value for tau I instead
% compute it as a general equilbrium condition. This is done just as a double check on rounding errors.

% Stationary distribution of eta (the exogenous state z, the efficiency labour process)
statdist_z=((ones(1,N_z)/N_z)*(pi_z^(10^6)))'; % Not really needed for anything. Just calculate it out of interest.

%% Now, create the return function
DiscountFactorParamNames={'beta','sj'};

ReturnFn=@(l,aprime,a,eta,r,tau,epsilon_j,I_j,Tr,SS,theta,alpha,delta,gamma,sigma) ConesaKrueger1999_ReturnFn(l,aprime,a,eta,r,tau,epsilon_j,I_j,Tr,SS,theta,alpha,delta,gamma,sigma)

%% Initial distribution of agents at birth (j=1)
jequaloneDist=zeros([n_a,n_z]);
jequaloneDist(1,:)=statdist_z; % All agents born with zero assets and with based on stationary distribution of the exogenous process on labour productivity units (comments elsewhere in CK1999 suggest that this is what they did, is not explicit in the paper)

%% Some things for final general equilibrium
% All three policy reforms have the same end point, namely completely terminate the social security systen:
Params.b_final=0;
% We know what the General Eqm values for tau and SS will be, so may as well set this as our initial guess.
Params.tau_final=0;
Params.SS_final=0;

%% FnsToEvaluate
FnsToEvaluate.K = @(d,aprime,a,z) a; % Aggregate assets K
FnsToEvaluate.N = @(d,aprime,a,z,I_j,epsilon_j) I_j*epsilon_j*z*d; % Aggregate labour supply (in efficiency units), CK1999 call this N
FnsToEvaluate.Tr_left = @(d,aprime,a,z,sj,n) (1-sj)*aprime; % Tr, accidental bequest transfers % This will later get divided by (1+n) due to population growth
FnsToEvaluate.FractionWorkingAge = @(d,aprime,a,z,I_j) I_j; % Fraction of population of working age (note that in principle this can calculated directly, but am computing it here as a way to double-check)
FnsToEvaluate.L = @(d,aprime,a,z) d; % Hours worked L (this is not needed to find the general eqm, but we will want it later)

%% General eqm variables: give some initial values
GEPriceParamNames={'r','tau','SS','Tr'};
Params.r=0.06; % interest rate on assets
Params.tau=Params.tau_initial; % Payroll tax rate.
Params.SS=0.4; % Benefits level for retirees
Params.Tr=0.3; % lumpsum transfers made out of the accidental bequests
% These were originally r=0.06, SS=1.2, Tr=1. Changed to present values as
% these are much closer to the solution and so substantially reduce the run
% time for finding General Eqm.

GeneralEqmEqns.CapitalMarket = @(r,K,N,theta,alpha,delta) r-(theta*(alpha)*(K^(alpha-1))*(N^(1-alpha))-delta); % Rate of return on assets is related to Marginal Product of Capital
GeneralEqmEqns.taxbalance = @(tau,FractionWorkingAge,b) tau-b*(1-FractionWorkingAge)/FractionWorkingAge; % From eqn 2.13 and 2.12 we get: tau=b*(frac of population retired)/(fraction of population of working age)
GeneralEqmEqns.SSbudgetbalance = @(SS,K,N,FractionWorkingAge,b,theta,alpha) SS-b*(theta*(1-alpha)*((K/N)^alpha))*N/FractionWorkingAge; % Social Security adds up, based on eqn (2.12) b*w*N/(fraction working age) [Note that because of how tau is calibrated (from b) this also means eqn (2.13) will hold.]
GeneralEqmEqns.AccBequests = @(Tr_left,Tr,n) Tr-Tr_left/(1+n); % Accidental bequests (adjusted for population growth) are equal to transfers received

%% The policy reforms for which we want to find the transition paths (note they are also needed before for finding the final general equilibrium)

% Number of time periods to allow for the transition (if you set T too low
% it will cause problems, too high just means run-time will be longer).
% (Problems in the sense that solution would be 'incorrect', it would likely still solve)
T=150; % CK1999 use T=150

if PolicyReform==1 % 1=Policy Reform A
    % "Beginning with period 2 the replacement rate is set equal to 0 and
    % stays there forever, i.e., b_1=0.5, b_t=0 for all t>1. This reform
    % terminates the social security system immediately and does not honor
    % entitlements to social security payments."
    ParamPath.b=zeros(T,1); % ParamPath is matrix of size T-by-'number of parameters that change over path'
elseif PolicyReform==2 % 2=Policy Reform B
    % "The replacement rate, 50%, is linearly reduced by one percentage
    % point over 50 periods, i.e., bt=0.5-0.01(t-1), t=1,2,...,50, bt=0 for
    % all t>50. This reform terminates the social security system gradually
    % so that entitlements are partially honored and payroll taxes are
    % accordingly reduced to finance progressively smaller benefits."
    ParamPath.b=[linspace(0.49,0,50)'; zeros(T-50,1)]; % ParamPath is matrix of size T-by-'number of parameters that change over path'
    % Note: Start at 0.49 as the "50 periods" of CK1999 includes the initial steady state prior to announcement.
elseif PolicyReform==3 %  3=Policy Reform C
    % "The replacement rate is fixed for 20 years at 50% and set to 0
    % thereafter, i.e., bt=0.5, t=1,2,...,20, bt=0 for all t>20. Therefore,
    % all individuals retired or about to retire keep their social security
    % benefits, but future retirees anticipate that they will receive only
    % part or no social security benefits. This reform allows agents to
    % readjust their plans for the anticipated reform in 20 years."
    ParamPath.b=[0.5*ones(19,1); zeros(T-19,1)]; % ParamPath is matrix of size T-by-'number of parameters that change over path'
    % Note: 19 as the "20 periods" of CK1999 includes the initial steady state prior to announcement.
end

% Because we are going to calculate the general equilibrium transition
% paths lets set it so that the general eqm (initial and final) are
% calculated to very high (likely excessive) levels of accuracy.
heteroagentoptions.toleranceGEprices=10^(-5); % Final eqm prices need to be highly accurate when using transition paths
heteroagentoptions.toleranceGEcondns=10^(-5); % Final eqm condns need to be highly accurate when using transition paths


%% Compute the initial general equilibrium and output based on it
[p_eqm_init,~, GeneralEqmEqnsValues_init]=HeteroAgentStationaryEqm_Case1_FHorz(jequaloneDist,AgeWeightsParamNames,n_d, n_a, n_z, N_j, 0, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, [],[],[], GEPriceParamNames, heteroagentoptions, simoptions, vfoptions);

Params.r=p_eqm_init.r;
Params.tau=p_eqm_init.tau;
Params.SS=p_eqm_init.SS;
Params.Tr=p_eqm_init.Tr;

% Some things about the initial general equilibrium  we will need for the transition path.
[V_init, Policy_init]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j, d_grid, a_grid, z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, [],vfoptions);
StationaryDist_init=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightsParamNames,Policy_init,n_d,n_a,n_z,N_j,pi_z,Params,simoptions);
AggVars_init=EvalFnOnAgentDist_AggVars_FHorz_Case1(StationaryDist_init, Policy_init, FnsToEvaluate, Params, [], n_d, n_a, n_z,N_j, d_grid, a_grid, z_grid); % The 2 is for Parallel (use GPU)

% Calculate the relevant entries for Table 5.
AggVars=EvalFnOnAgentDist_AggVars_FHorz_Case1(StationaryDist_init, Policy_init, FnsToEvaluate, Params, [], n_d, n_a, n_z,N_j, d_grid, a_grid, z_grid); 
AllStats=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDist_init, Policy_init, FnsToEvaluate, Params, [], n_d, n_a, n_z,N_j, d_grid, a_grid, z_grid); 

StationaryEqmStats(1).b=Params.b;
StationaryEqmStats(1).r=Params.r;
StationaryEqmStats(1).K=AggVars.K.Mean;
StationaryEqmStats(1).N1=AggVars.N.Mean;
KdivN1=AggVars.K.Mean/AggVars.N.Mean;
KdivN2=((Params.r+Params.delta)/(Params.theta*Params.alpha))^(1/(Params.alpha-1)); % This is just to double check
StationaryEqmStats(1).w=Params.theta*(1-Params.alpha)*(KdivN1^Params.alpha);
StationaryEqmStats(1).h=100*AggVars.L.Mean/AggVars.FractionWorkingAge.Mean; % total hours by working age divided by fraction working age. (multiply by 100 to express as percent rather than as fraction)
y=Params.theta*(AggVars.K.Mean^Params.alpha)*(AggVars.N.Mean^(1-Params.alpha));
StationaryEqmStats(1).y=y;
StationaryEqmStats(1).Kdivy=AggVars.K.Mean/y;
StationaryEqmStats(1).SSdivy=Params.SS/y;
StationaryEqmStats(1).HoursWorked=AggVars.L.Mean;
StationaryEqmStats(1).CoeffOfVariance_a=AllStats.K.StdDev/AllStats.K.Mean; % Standard deviation/Mean

% Statistics just for working age
simoptions.agegroupings=[1,Params.Jr]; % Working age, Retired (not actually interested in the numbers for retired)
AllEmployedStats=LifeCycleProfiles_FHorz_Case1(StationaryDist_init,Policy_init,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions);

% Not sure from reading CK1999, but seems to be labour hours?
StationaryEqmStats(1).CoeffOfVariance_l=sqrt(AllEmployedStats.L.Variance(1))/AllEmployedStats.L.Mean(1); % the 1 refers to the first 'agegroupings' which is working age

% For Figure 1 & 2, need life cycle profile of average hours worked and assets.
% [Actual CK1999 figure looks like it may in fact use 5 year age bins, but I do individual ages here.]
simoptions.agegroupings=1:1:Params.J; % This is actually the default, but need to set it as previously changed it to working age
LifeCycleProfiles_init=LifeCycleProfiles_FHorz_Case1(StationaryDist_init,Policy_init,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions);

%% Compute the final General Equilibrium
% All three policy reforms have the same end point, namely completely terminate the social security systen:
Params.b=Params.b_final;
% We know what the General Eqm value for SS will be, so may as well set this as our initial guess.
Params.tau=Params.tau_final;
Params.SS=Params.SS_final;

[p_eqm_final,~, GeneralEqmEqnsValues_final]=HeteroAgentStationaryEqm_Case1_FHorz(jequaloneDist,AgeWeightsParamNames,n_d, n_a, n_z, N_j, 0, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, [], [], [], GEPriceParamNames, heteroagentoptions, simoptions, vfoptions);

Params.r=p_eqm_final.r;
Params.tau=p_eqm_final.tau;
Params.SS=p_eqm_final.SS;
Params.Tr=p_eqm_final.Tr;

% tau and ss should be zero, so if they are close to it then just set them
% to exactly zero (I only do this if close so that otherwise it is clear there is a problem)
if abs(Params.tau)<10^(-3)
    Params.tau=0;
end
if abs(Params.SS)<10^(-3)
    Params.SS=0;
end

% Some things about the final general equilibrium  we will need for the transition path.
[V_final, Policy_final]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j, d_grid, a_grid, z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, [],vfoptions);
StationaryDist_final=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightsParamNames,Policy_final,n_d,n_a,n_z,N_j,pi_z,Params,simoptions);
AggVars_final=EvalFnOnAgentDist_AggVars_FHorz_Case1(StationaryDist_final, Policy_final, FnsToEvaluate, Params, [], n_d, n_a, n_z,N_j, d_grid, a_grid, z_grid); 

% Calculate the relevant entries for Table 5.
AggVars=EvalFnOnAgentDist_AggVars_FHorz_Case1(StationaryDist_final, Policy_final, FnsToEvaluate, Params, [], n_d, n_a, n_z,N_j, d_grid, a_grid, z_grid);
AllStats=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDist_final, Policy_final, FnsToEvaluate, Params, [], n_d, n_a, n_z,N_j, d_grid, a_grid, z_grid);

StationaryEqmStats(2).b=Params.b;
StationaryEqmStats(2).r=Params.r;
StationaryEqmStats(2).K=AggVars.K.Mean;
KdivN1=AggVars.K.Mean/AggVars.N.Mean;
KdivN2=((Params.r+Params.delta)/(Params.theta*Params.alpha))^(1/(Params.alpha-1)); % This is just to double check
StationaryEqmStats(2).w=Params.theta*(1-Params.alpha)*(KdivN1^Params.alpha);
StationaryEqmStats(2).h=100*AggVars.L.Mean/AggVars.FractionWorkingAge.Mean; % total hours by working age divided by fraction working age. (multiply by 100 to express as percent rather than as fraction)
y=Params.theta*(AggVars.K.Mean^Params.alpha)*(AggVars.N.Mean^(1-Params.alpha));
StationaryEqmStats(2).y=y;
StationaryEqmStats(2).Kdivy=AggVars.K.Mean/y;
StationaryEqmStats(2).SSdivy=Params.SS/y;
StationaryEqmStats(2).HoursWorked=AggVars.L.Mean;

StationaryEqmStats(2).CoeffOfVariance_a=AllStats.K.StdDev/AllStats.K.Mean; % Standard deviation/Mean

% Statistics just for working age
simoptions.agegroupings=[1,Params.Jr]; % Working age, Retired (not actually interested in the numbers for retired)
AllEmployedStats=LifeCycleProfiles_FHorz_Case1(StationaryDist_final,Policy_final,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions);

% Not sure from reading CK1999, but seems to be labour hours?
StationaryEqmStats(2).CoeffOfVariance_l=sqrt(AllEmployedStats.L.Variance(1))/AllEmployedStats.L.Mean(1); % One refers to the first 'agegroupings' which is working age

EV_SS_numerator=sum(V_final(1,:,1).*StationaryDist_final(1,:,1)); % 0 assets and age 1
EV_SS_denominator=sum(V_init(1,:,1).*StationaryDist_init(1,:,1)); % 0 assets and age 1

StationaryEqmStats(2).EV_SS=(EV_SS_numerator/EV_SS_denominator)^(1/(Params.gamma*(1-Params.sigma)))-1;

% For Figure 1 & 2, need life cycle profile of average hours worked and assets.
% [Actual CK1999 figure looks like it may in fact use 5 year age bins, but I do individual ages here.]
simoptions.agegroupings=1:1:Params.J; % This is actually the default, but am setting it here explicitly anyway.
LifeCycleProfiles_final=LifeCycleProfiles_FHorz_Case1(StationaryDist_final,Policy_final,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions);


% % %% Figure of Life-Cycle Profiles
% % % Relate to Figures 1 and 2 of the paper
% % figure(1)
% % age_axis_Jr=19+(1:1:Params.Jr); % Just working age
% % subplot(2,1,1); plot(age_axis_Jr, LifeCycleProfiles_init.L.Mean(1:Params.Jr)) 
% % hold on
% % plot(age_axis_Jr, LifeCycleProfiles_final.L.Mean(1:Params.Jr),'--')
% % hold off
% % title('Life-Cycle Profile of hours worked')
% % ylim([0,0.5])
% % xlabel('Age')
% % legend('initial','final')
% % subplot(2,1,2); plot(age_axis_J, LifeCycleProfiles_init.K.Mean) 
% % hold on
% % plot(age_axis_J, LifeCycleProfiles_final.K.Mean,'--')
% % hold off
% % title('Life-Cycle Profile of assets')
% % ylim([0,15])
% % xlabel('Age')


%% Now, the transition path itself
% For this we need the following extra objects: PricePathOld, PriceParamNames, ParamPath, ParamPathNames, T, V_final, StationaryDist_init
% (already calculated V_final & StationaryDist_init above)
Params.r=p_eqm_init.r;
Params.tau=p_eqm_init.tau;
Params.SS=p_eqm_init.SS;
Params.Tr=p_eqm_init.Tr;

p_init=p_eqm_init; %0.75*p_eqm_init(2)]; % Ignoring the behavioural responses we would expect revenues (and thus spending) to fall to 0.75 of what they were, so lets use this as a 'starting guess'
p_final=p_eqm_final;

% We need to give an initial guess for the price path on interest rates
PricePath0.r=[linspace(p_init.r, p_final.r, floor(T/2))'; p_final.r*ones(T-floor(T/2),1)];
PricePath0.tau=p_init.tau*(2*ParamPath.b); % Because tau is so closely related to b this seems a reasonable guess (the 2* is so the init value of b of 0.5 is instead 1)
PricePath0.SS=p_init.SS*(2*ParamPath.b); % Because SS is so closely related to b this seems a reasonable guess (the 2* is so the init value of b of 0.5 is instead 1)PricePath0.SS=[linspace(p_init.SS, p_final.SS, floor(T/2))'; p_final.SS*ones(T-floor(T/2),1)];
PricePath0.Tr=[linspace(p_init.Tr, p_final.Tr, floor(T/2))'; p_final.Tr*ones(T-floor(T/2),1)];

% The following transition path general equilibrium conditions for are just
% the same as those the stationary model, but in some modeles they would differ
GeneralEqmEqns_Transition.CapitalMarket = @(r,K,N,theta,alpha,delta) r-(theta*(alpha)*(K^(alpha-1))*(N^(1-alpha))-delta); % Rate of return on assets is related to Marginal Product of Capital
GeneralEqmEqns_Transition.taxbalance = @(tau,FractionWorkingAge,b) tau-b*(1-FractionWorkingAge)/FractionWorkingAge; % From eqn 2.13 and 2.12 we get: tau=b*(frac of population retired)/(fraction of population of working age)
GeneralEqmEqns_Transition.SSbudgetbalance = @(SS,K,N,FractionWorkingAge,b,theta,alpha) SS-b*(theta*(1-alpha)*((K/N)^alpha))*N/FractionWorkingAge; % Social Security adds up, based on eqn (2.12) b*w*N/(fraction working age) [Note that because of how tau is calibrated (from b) this also means eqn (2.13) will hold.]
GeneralEqmEqns_Transition.AccBequests = @(Tr_left_tminus1,Tr,n) Tr-Tr_left_tminus1/(1+n); % Accidental bequests (from last period, and adjusted for population growth) are equal to transfers received
% VFI Toolkit understands '_tminus1' to mean the value from the previous period.

% To use a '_tminus1' variable we must include its initial value
% To be able to use a '_tminus1' as input, we have to provide the initial value
transpathoptions.initialvalues.Tr_left=p_eqm_init.Tr*(1+Params.n); % Tr is received, so what was left last period is (1+n) times this.


% Setup the options relating to the transition path
transpathoptions.verbose=1;
transpathoptions.maxiterations=200; % default is ???
vfoptions.policy_forceintegertype=2;

transpathoptions.GEnewprice=3;
% Need to explain to transpathoptions how to use the GeneralEqmEqns to
% update the general eqm transition prices (in PricePath).
transpathoptions.GEnewprice3.howtoupdate=... % a row is: GEcondn, price, add, factor
    {'CapitalMarket','r',0,0.1;...  % CaptialMarket GE condition will be positive if r is too big, so subtract
    'taxbalance','tau',0,0.1;... % taxbalance GE condition will be positive if tau is too big, so subtract
    'SSbudgetbalance','SS',0,0.1;... % taxbalance GE condition will be positive if tau is too big, so subtract
    'AccBequests','Tr',0,0.1;... % taxbalance GE condition will be positive if tau is too big, so subtract
    };
% Note: the update is essentially new_price=price+factor*add*GEcondn_value-factor*(1-add)*GEcondn_value
% Notice that this adds factor*GEcondn_value when add=1 and subtracts it what add=0
% A small 'factor' will make the convergence to solution take longer, but too large a value will make it 
% unstable (fail to converge). Technically this is the damping factor in a shooting algorithm.


fprintf('Now computing the Transition path itself \n')
PricePath=TransitionPath_Case1_FHorz(PricePath0, ParamPath, T, V_final, StationaryDist_init, n_d, n_a, n_z, N_j, pi_z, d_grid,a_grid,z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns_Transition, Params, DiscountFactorParamNames, AgeWeightsParamNames, transpathoptions, simoptions, vfoptions);

%% Now calculate some things about the transition path (The path for Value fn, Policy fn, and Agent Distribution)

% You can calculate the value and policy functions for the transition path
[VPath,PolicyPath]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, V_final, Policy_final, Params, n_d, n_a, n_z, N_j, pi_z, d_grid, a_grid,z_grid, DiscountFactorParamNames, ReturnFn, transpathoptions);

% You can then use these to calculate the agent distribution for the transition path
AgentDistPath=AgentDistOnTransPath_Case1_FHorz(StationaryDist_init,PricePath, ParamPath, PolicyPath, AgeWeightsParamNames,n_d,n_a,n_z,N_j,pi_z,T, Params);


%% For Figures 3, 9, and 12
% Create graphs of output per capita, interest rates, capital-output ratio,
% and labor supply over the transition path.

Figures_TransitionAggregates.interestrate=[p_eqm_init.r; PricePath.r];

% FOLLOWING TWO LINES WERE ABOUT DEBUGGING
fprintf('Check if this is zero %8.8f \n',max(max(max(max(max(abs(Policy_final-round(Policy_final))))))))
fprintf('Check if this is nonzero %8.8f \n',min(min(min(min(Policy_final)))))

AggVarsPath=EvalFnOnTransPath_AggVars_Case1_FHorz(FnsToEvaluate, AgentDistPath,PolicyPath, PricePath, ParamPath, Params, T, n_d, n_a, n_z, N_j, pi_z, d_grid, a_grid,z_grid, DiscountFactorParamNames, transpathoptions);

K_Path=[StationaryEqmStats(1).K; AggVarsPath.K.Mean];
N1_Path=[StationaryEqmStats(1).N1; AggVarsPath.N.Mean];
L_Path=[StationaryEqmStats(1).HoursWorked; AggVarsPath.L.Mean];
y_Path=Params.theta*(K_Path.^Params.alpha).*(N1_Path.^(1-Params.alpha));

Figures_TransitionAggregates.outputpercapita=y_Path; % Output per capita is just equal to output as mass of agents is normalized to one.
Figures_TransitionAggregates.capitaloutputratio=K_Path./y_Path;
Figures_TransitionAggregates.laborsupply=N1_Path;
Figures_TransitionAggregates.hoursworked=L_Path;

%% Create graphs of output per capita, interest rates, capital-output ratio, and labor supply over the transition path.
% Relates to figures 3, 9, 12 of the paper

% % figure(2)
% % subplot(2,2,1)
% % plot(0:1:T, Figures_TransitionAggregates.outputpercapita)
% % title('Evolution of Output per capita')
% % ylim([1,2])
% % xlabel('Time')
% % subplot(2,2,2)
% % plot(0:1:T, Figures_TransitionAggregates.interestrate)
% % title('Evolution of Interest Rate')
% % ylim([0,0.08])
% % xlabel('Time')
% % subplot(2,2,3)
% % plot(0:1:T, TransitionAggregates.capitaloutputratio)
% % title('Evolution of Capital-Output Ratio')
% % ylim([2.5,5])
% % xlabel('Time')
% % subplot(2,2,4)
% % plot(0:1:T, 100*Figures_TransitionAggregates.hoursworked)
% % title('Evolution of Labor Supply (hours worked)')
% % % ylim([22,32])
% % xlabel('Time')


%% Create graph of the share of the population that 'votes' in favour (that is, the share for whom Equivalent Variation of the reform is positive)
% Relates to figures 4, 10, 13 of the paper
% Note that figures 5, 6, 7 of the paper would also just graph the EquivVariation calculated here but in different ways

% (consumption equivalent variation) for each point on the grid between the initial period and the first period of the transition (CK1999, bottom pg
% 764). This is eqn (3.1) of CK1999. We already have V_init, which in the notation of eqn (3.1) is v_1
% So we just need to compute v_2, the value fn in the first period of the transition.

% Note: we actually already calculated these, but I will do it again anyway
[VPath,PolicyPath]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, V_final, Policy_final, Params, n_d, n_a, n_z, N_j, pi_z, d_grid, a_grid,z_grid, DiscountFactorParamNames, ReturnFn, transpathoptions);
EquivVariation=(VPath(:,:,:,1)./V_init).^(1/(Params.gamma*(1-Params.sigma)))-1;

% % figure(3)
% % % Start with the 'population share' line.
% % populationshares=shiftdim(sum(sum(StationaryDist_init,1),2),2); % Note, in principle could just use mewj, but lets do this in following manner
% % plot(age_axis_J, populationshares,':');
% % % double-check: should equal Params.(AgeWeightsParamNames{1})
% % VoteInFavor_nohet=(EquivVariation>=0); % agents who are indifferent are assumed to vote in favour
% % VoteInFavor_symhet=(EquivVariation>=0); % agents who are indifferent are assumed to vote in favour
% % VoteInFavor_asymhet=(EquivVariation>=0); % agents who are indifferent are assumed to vote in favour
% % ShareOfVotesInFavor_nohet=shiftdim(sum(sum(VoteInFavor_nohet.*StationaryDist_init,1),2),2);
% % ShareOfVotesInFavor_symhet=shiftdim(sum(sum(VoteInFavor_symhet.*StationaryDist_init,1),2),2);
% % ShareOfVotesInFavor_asymhet=shiftdim(sum(sum(VoteInFavor_asymhet.*StationaryDist_init,1),2),2);
% % hold on
% % plot(age_axis_J, ShareOfVotesInFavor_nohet, age_axis_J, ShareOfVotesInFavor_symhet,'--',age_axis_J, ShareOfVotesInFavor_asymhet,'-.')
% % hold off
% % legend('Pop. Share','No Het','Sym Het','Asym Het')
% % xlabel('age')
% % ylabel('Percentage of total votes')

%% Create graph of the share of the population that 'votes' in favour in partial equilbirium (if we do not allow for the prices to change as they do in the general equilibrium analysis)
% Partial equilibrium transition path, with prices fixed at their inital values.
PricePath_Fixed.r=p_init.r*ones(T,1);
PricePath_Fixed.tau=p_final.tau*ones(T,1); % Note that this is will be zero
PricePath_Fixed.SS=p_final.SS*ones(T,1); % Note that this is will be zero
PricePath_Fixed.Tr=p_init.Tr*ones(T,1);
% Note that will uses the ParamPath the reflects the policy changes.

[VPath,PolicyPath]=ValueFnOnTransPath_Case1_FHorz(PricePath_Fixed, ParamPath, T, V_final, Policy_final, Params, n_d, n_a, n_z, N_j, pi_z, d_grid, a_grid,z_grid, DiscountFactorParamNames, ReturnFn, transpathoptions);
EquivVariation_FixedPrices=(VPath(:,:,:,1)./V_init).^(1/(Params.gamma*(1-Params.sigma)))-1;

% % figure(4)
% % % The paper plots this in Figure 8
% % hold on
% % Y1=plot(a_grid, EquivVariation_FixedPrices(:,1,20-19),'r');
% % Y2=plot(a_grid, EquivVariation_FixedPrices(:,1,30-19),'b');
% % Y3=plot(a_grid, EquivVariation_FixedPrices(:,1,60-19),'g');
% % plot(a_grid, zeros(size(a_grid)),'k')
% % plot(a_grid, EquivVariation_FixedPrices(:,2,20-19),'--r', a_grid, EquivVariation_FixedPrices(:,2,30-19),'--b', a_grid, EquivVariation_FixedPrices(:,2,60-19),'--g')
% % % Requires setting up hidden lines with the appropriate style
% % H1 = plot(a_grid,EquivVariation_FixedPrices(:,1,20-19), '-', 'LineWidth', 2, 'Color', 'k', 'Visible', 'off'); % solid line (invisible and black)
% % H2 = plot(a_grid,EquivVariation_FixedPrices(:,1,20-19), '--', 'LineWidth', 2, 'Color', 'k', 'Visible', 'off'); % dashed line (invisible and black)
% % hold off
% % legend([Y1,Y2,Y3,H1,H2],'Age 20', 'Age 30', 'Age 60','Bad Shock','Good Shock');
% % title('Welfare effects of Reform A: Asymm. heterogeneity')
% % ylim([-0.8,0.2])
% % xlabel('Asset Position')
% % ylabel('Cons. Equiv. Var.')


