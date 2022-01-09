% Example based on Imrohoroglu, Imrohoroglu & Joines (1995) - A life cycle analysis of social security
%
% Solves the baseline model and extensions of it to include 'medical shocks' and 'deterministic growth'.
% See 'Example' for a slightly simpler code that just solves the baseline model.
%
% No decision variable: IIJ1995 has inelastic labor supply
% One endogenous state variable: assets (k)
% One stochastic exogenous state variable (z, which is either 'employed' or 'unemployed')
% Age

Params.J=65; % Number of period in life-cycle (represents ages 21 to 85)

% Grid sizes to use
% % labour supply is exogenous (IIJ1995 refers to this as 'perfectly inelastic labor supply')
n_k=1251; % Assets
n_z=2; % employed or unemployed (IIJ1995 refers to this shock as Markow, but in calibration actually imposes that it is iid).
N_j=Params.J; % Number of periods in finite horizon

%% Set model parameters

%Preference parameters
Params.gamma=2; % IIJ1995, pg 91: "1/gamma=0.5"
Params.beta=1.011; % Note that once this is combined with the conditional probability of survival it gives numbers that are typically less than 1.
Params.h=0.45; % IIJ1995, pg 90

%Technology parameters
Params.A=1.3193; % What IIJ1995 call 'B'
Params.alpha=0.36; % Note that I use this as the parameter on capital in the Cobb-Douglas production fn while IIJ1995 used 1-alpha (they use alpha as the parameter on labor input)
Params.delta=0.08;

% People become economically active (j=1) at age 21, retire at 65, and max age is 85
Params.J=N_j; %=85-21
Params.Jr=45; % =65-21 % I1998 refers to this as retirement age j*
Params.I_j=[1:1:Params.Jr,zeros(1,Params.J-Params.Jr)]; % Indicator for retirement

% % Population growth of 1.2% IIJ1995, pg 90
Params.n=0.012;
% % Social Security replacement rate (which IIJ1995 refer to as theta)
Params.b=0.3; % IIJ1995, pg 92: "We search over values of [b between 0 and 1]", I set baseline as 0.3 which is they value they find to be optimal. IIJ1995 calls this theta. (what IIJ1995 calls b, is here called 'SS', the actual social security pension benefits received)
% % Note: social security is based on "average lifetime employed wage", not
% % the actual average lifetime earnings which would account for whether or
% % not you were employed (and would need to be kept track of as an extra endogenous state).
% Unemployment replacement rate
Params.zeta=0.4; % IIJ1995, pg 90

% tau_s and tau_u are set below as part of general eqm

%Probability of surviving, conditional on being alive (Imrohoroglu, Imrohoroglu & Joines 1995 mentions getting them from Faber, 1982)
% The commented out lines below contain my original version using more recent numbers from 'same' source.
% I was sent the original survival probabilities by email by Selahattin
% Imrohoroglu and the contents of that 'psurv.dat' file are used below.
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
dj_temp=interp1([0,30,60,65,70,100],[0.00587,0.00116,0.01086,0.01753,0.02785,0.39134],0:1:100,'linear');
Params.sj=1-dj_temp(21:85);
Params.sj(1)=1;
Params.sj(end)=0;
% % I use linear interpolation to fill in the numbers inbetween those reported by Bell & Miller (2005).
% % I have aditionally imposed that the prob of death at age 20 be zero and that prob of death at age 85 is one.
% Following are the orginal survival probabilites I was emailed by Selahattin Imrohoroglu in file 'psurv.dat'
Params.sj=[1.0000000, 0.99851000, 0.99844000, 0.99838000, 0.99832000, 0.99826000, 0.99820000, 0.99816000, 0.99815000, 0.99819000,...
    0.99826000, 0.99834000,0.99840000, 0.99843000, 0.99841000, 0.99835000, 0.99828000, 0.99818000, 0.99807000, 0.99794000,...
    0.99778000, 0.99759000, 0.99737000, 0.99712000, 0.99684000, 0.99653000, 0.99619000, 0.99580000, 0.99535000, 0.99481000,...
    0.99419000, 0.99350000, 0.99278000, 0.99209000, 0.99148000, 0.99088000, 0.99021000, 0.98942000, 0.98851000, 0.98746000,...
    0.98625000, 0.98495000, 0.98350000, 0.98178000, 0.97974000, 0.97743000, 0.97489000, 0.97226000, 0.96965000, 0.96715000,...
    0.96466000, 0.96200000, 0.95907000, 0.95590000, 0.95246000, 0.94872000, 0.94460000, 0.94017000, 0.93555000, 0.93077000, ...
    0.92570000, 0.92030000, 0.91431000, 0.90742000, 0.89948000];
Params.sj(end)=0;


%Age profile of productivity (based on lifetime profile of earnings, Hansen 1993), Epsilon_j. 
% The following line is the contents of 'comboeff.dat' file I received from Selahattin Imrohoroglu (see a few lines above)
Params.epsilon_j=[0.36928407;  0.42120465;  0.49398951; 0.56677436; 0.63955922; 0.71234407; 0.78512892; 0.81551713; 0.84590532; 0.87629354; 0.90668176; 0.93706995; 0.96379826; 0.99052656; 1.0172549; 1.0439831; 1.0707114; 1.0819067; 1.0931018; 1.1042970; 1.1154923; 1.1266874; 1.1435058; 1.1603241; 1.1771424; 1.1939607; 1.2107789; 1.2015913; 1.1924037; 1.1832160; 1.1740283; 1.1648407; 1.1609988; 1.1571568; 1.1533149; 1.1494730; 1.1456310; 1.1202085; 1.0947859; 1.0693634; 1.0439408; 1.0185183; 0.99309576; 0.96767321];
Params.epsilon_j=[Params.epsilon_j; zeros(Params.J-Params.Jr+1,1)]; % Fill in the zeros for retirement age
Params.epsilon_j=Params.epsilon_j/mean(Params.epsilon_j(1:Params.Jr)); % IIJ1995 pg 90, "normalized to average unity over the 44 working years".
% IIJ1995 bans retirees from working. One way to do this is just to set epsilon_j=0 for retirees. 
% Rather than do this via epsilon_j it is done by I_j which in practice is multiplied by epsilon_j.
Params.I_j=[ones(Params.Jr-1,1);zeros(Params.J-Params.Jr+1,1)];
% Note that I_j is redundant using these epsilon_j numbers from IIJ1995 as these are anyway zero during retirement.

% The following is the social security benefits divided by the wage
Params.SSdivw=Params.b*sum(Params.h*Params.epsilon_j(1:(Params.Jr-1)))/(Params.Jr-1);

% The following is lump-sum transfers which are needed to calculate welfare benefits of reforms.
Params.LumpSum=0;  % They are by definition zero whenever not calculating the welfare benefits of reforms

%% Grids and exogenous shocks
k_grid=18*(linspace(0,1,n_k).^3)'; % Grid on capital

z_grid=[0; 1]; % Indicates 'employed'
pi_z=[0.06, 0.94; 0.06, 0.94]; % Note that IIJ1995 reports this on pg 92 with the 'reversed' definition of which state is which.
% This choice implies that the expected duration of unemployment is 1/(1-0.06)=1.0638 periods.
% Notice this actually makes the exogenous shock iid (rows of the markov transition matrix are identical)
statdist_z=pi_z(1,:)'; % follows from iid

%% Beyond the baseline model IIJ1995 have two 'extensions', medical shocks and (deterministic) productivity growth

Params.agej=1:1:Params.J; % We will need this for the medical shocks

% Medical Shocks (IIJ1995, pg 109)
Params.MedicalShock=0; % 0 is the baseline model with no medical shocks. (1 and 2 are the two different specifications of Medical shock considered by IIJ1995, see their pages 110-1.
if Params.MedicalShock>0
   vfoptions.ExogShockFnParamNames={'agej','Jr','MedicalShock'}; % Notice that I_j, which indicates working age, does NOT contain the required information to know when to switch from the shock process being about employment status to being about medical shocks is actually occuring (it just indicates if it has or not), hence we need to use agej and Jr instead.
   simoptions.ExogShockFnParamNames=vfoptions.ExogShockFnParamNames;
   vfoptions.ExogShockFn=@(agej,Jr,MedicalShock) ImrohorogluImrohorogluJoines1995_ExogShockFn(agej,Jr,MedicalShock);
   simoptions.ExogShockFn=vfoptions.ExogShockFn;
end

% (Deterministic labour-augmenting) Productivity Growth (IIJ1995, pg 108)
% Changes the model in two ways: the discount factor and the age-earnings profile
Params.g=0;
% The following captures how productivity growth changes the discount factor
Params.gdiscount=(1+Params.g)^(1-Params.gamma);
% % Modify the age-earnings profile
% Params.epsilon_j=Params.epsilon_j*((1+g).^((1:1:Params.J)-1));
% Params.epsilon_j=Params.epsilon_j/mean(Params.epsilon_j(1:Params.Jr)); % Renormalize to one over working ages (this is necessary due to how some of the replacement rates for various benefits are calculated)
% % Note: Adding productivity growth in this way means it will also have to
% % be accounted for when computing some statistics from the agent
% % distribution in equilibrium (e.g., anything which contains a benefit that
% % is proportional to lifetime earnings; see IIJ pg. 109).

% Note: When there is productivity growth all model output is being detrended to measure it in terms 
% of the 'first year'. You would need to unnormalize it [multiply it by mean((1+g).^((1:1:Params.Jr)-1)) ]
% if you wanted to do something like a transition path, time series or panel data simulation, or
% life-cycle profiles.
% This is exactly the same detrending you do to solve the Solow growth model (switching things into 
% 'per capita per productivity unit' terms), just that we are doing it in a more complex model.
Params.workinglifeincome=mean(Params.epsilon_j(1:Params.Jr)'.*((1+Params.g).^((1:1:Params.Jr)-1))); % Note that when g=0 this is one. We need this for calculating pensions which are defined as a specific fraction of this.

%% Get into format for VFI Toolkit
d_grid=[];
a_grid=k_grid;
% z_grid

n_d=0;
n_a=n_k;
% n_z

%% Calculate the population distribution across the ages (population at time t is made stationary by 
% Population at time t is made stationary by  dividing it by \Pi_{i=1}^{t} (1+n_{i}) 
% (the product of all growth rates up till current period)
% The VFI Toolkit needs to give it a name so that it can be automatically used when calculating model outputs.
Params.mewj=ones(Params.J,1);
for jj=2:Params.J
    Params.mewj(jj)=Params.mewj(jj-1)*(1/(1+Params.n))*Params.sj(jj-1);
end
Params.mewj=Params.mewj/sum(Params.mewj); % normalize to measure one
Params.mewj=Params.mewj'; % Age weights must be a row vector.

AgeWeightsParamNames={'mewj'}; % Many finite horizon models apply different weights to different 'ages'; eg., due to survival rates or population growth rates.

%% General equilibrium, and initial parameter values for the general eqm parameters

GEPriceParamNames={'r','tau_u','tau_s','Tr_beq'};
Params.r=0.06; % interest rate on assets
Params.tau_u=0.019; % set to balance unemployment benefits expenditure
Params.tau_s=0.122; % set to balance social security benefits expenditure
Params.Tr_beq=0.2; % Accidental bequests (IIJ1995 calls this T)

%% Now, create the return function 
DiscountFactorParamNames={'beta','sj','gdiscount'}; % gdiscount is 1 in the baseline, it is needed for the extension to include deterministic productivity growth
 
ReturnFn=@(aprime,a,z,r,tau_u, tau_s,gamma, h,zeta, epsilon_j,I_j,alpha,delta, A,SSdivw, Tr_beq, MedicalShock,workinglifeincome,g,agej,LumpSum) ImrohorogluImrohorogluJoines1995_ReturnFn(aprime,a,z,r,tau_u, tau_s,gamma, h,zeta, epsilon_j,I_j,alpha,delta, A,SSdivw, Tr_beq, MedicalShock,workinglifeincome,g,agej,LumpSum)
% Note: MedicalShock,workinglifeincome,g,agej are only required for the extensions of the model. None of these would be needed for just the baseline model.

%% Now solve the value function iteration problem, just to check that things are working before we go to General Equilbrium

disp('Test ValueFnIter')
vfoptions=struct(); % Use the defaults
[V, Policy]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j, d_grid, a_grid, z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, [],vfoptions);

% max(max(max(max(Policy))))<n_a % Double check that never try to leave top of asset grid.
% sum(sum(sum(sum(Policy==n_a))))

%% Initial distribution of agents at birth (j=1)
jequaloneDist=zeros(n_a,n_z); 
jequaloneDist(1,:)=statdist_z; % All agents born with zero assets and with based on stationary distribution of the exogenous process on labour productivity units (I1998, pg 326)

%% Test
disp('Test StationaryDist')
simoptions=struct(); % Just use the defaults
StationaryDist=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightsParamNames,Policy,n_d,n_a,n_z,N_j,pi_z,Params,simoptions);


%% Set up the General Equilibrium conditions (on assets/interest rate, assuming a representative firm with Cobb-Douglas production function)

% Function to Evaluate (to get aggregates, or life-cycle profiles, etc.)
FnsToEvaluate.K = @(aprime,a,z) a; % Aggregate assets K
FnsToEvaluate.N = @(aprime,a,z,I_j,h,epsilon_j) I_j*h*epsilon_j*z; % Aggregate effective labour supply (in efficiency units), I1998 calls this N
FnsToEvaluate.Tr = @(aprime,a,z,sj,n) (1-sj)*aprime/(1+n); % Tr, accidental bequest transfers % The divided by (1+n) is due to population growth and because this is Tr_{t+1}
FnsToEvaluate.C = @(aprime,a,z,r,tau_u, tau_s,h,zeta,epsilon_j,I_j,alpha,delta, A,SSdivw, Tr_beq,workinglifeincome,g,agej,MedicalShock,LumpSum) ImrohorogluImrohorogluJoines1995_ConsumptionFn(aprime,a,z, r,tau_u, tau_s,h,zeta,epsilon_j,I_j,alpha,delta, A,SSdivw, Tr_beq,workinglifeincome,g,agej,MedicalShock,LumpSum); % Consumption
FnsToEvaluate.LaborIncome = @(aprime,a,z,r,A,alpha,delta,h,I_j,epsilon_j) ((1-alpha)*(A^(1/(1-alpha)))*((r+delta)/alpha)^(alpha/(alpha-1)))*I_j*h*epsilon_j*z; % Labour income (this is also the tax base for various taxes)
FnsToEvaluate.PotentialLaborIncome = @(aprime,a,z,r,A,alpha,delta,h,I_j) ((1-alpha)*(A^(1/(1-alpha)))*((r+delta)/alpha)^(alpha/(alpha-1)))*I_j*h*(1-z); % Potential Labour income (in worked 100 percent of time) of the unemployed (the unemployment benefits are set as a fraction, zeta, of this)
FnsToEvaluate.SSbenefits = @(aprime,a,z,SSdivw,r,A,alpha,delta,I_j) ((1-alpha)*(A^(1/(1-alpha)))*((r+delta)/alpha)^(alpha/(alpha-1)))*SSdivw*(1-I_j); % Total social security benefits: w*SSdivw*(1-I_j)

% General Equilibrium Equations
% Recall that GEPriceParamNames={'r','tau_u','tau_s','Tr_beq'}; In following lines p is the vector of these and so, e.g., p(2) is G.
GeneralEqmEqns.capitalmarket = @(N,K,r,A,alpha,delta) r-(A*(alpha)*(K^(alpha-1))*(N^(1-alpha))-delta); % Rate of return on assets is related to Marginal Product of Capital
GeneralEqmEqns.UnempGenerosity = @(tau_u,LaborIncome,zeta,PotentialLaborIncome) tau_u*LaborIncome-zeta*PotentialLaborIncome; % Equation (16)
GeneralEqmEqns.SSbalance = @(tau_s,LaborIncome,SSbenefits) tau_s*LaborIncome-SSbenefits; % Equation (15)
GeneralEqmEqns.Bequests = @(Tr_beq,Tr) Tr_beq-Tr; % Equation (17): Accidental bequests (adjusted for population growth) are equal to transfers received (this is essentially eqn 14)
% Note that the second of these general eqm conditions, about tau_u and zeta is in fact 
% unnecessary with exogenous labor (as w would cancel from both terms and the rest is exogenous).

%% Test
disp('Test AggVars')
AggVars=EvalFnOnAgentDist_AggVars_FHorz_Case1(StationaryDist, Policy, FnsToEvaluate, Params, [], n_d, n_a, n_z,N_j, d_grid, a_grid, z_grid,[],simoptions);

%% Calculate the general equilibrium
heteroagentoptions.verbose=1; % Give info on how the General eqm conditions are going
[p_eqm,p_eqm_index, GeneralEqmEqnsValues]=HeteroAgentStationaryEqm_Case1_FHorz(jequaloneDist,AgeWeightsParamNames,n_d, n_a, n_z, N_j, 0, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, [], [], [], GEPriceParamNames,heteroagentoptions,simoptions,vfoptions);
Params.r=p_eqm.r;
Params.tau_u=p_eqm.tau_u;
Params.tau_s=p_eqm.tau_s;
Params.Tr_beq=p_eqm.Tr_beq;

%% Calculate some model statistics of the kind reported in IIJ1995
[V, Policy]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j, d_grid, a_grid, z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, [],vfoptions);
StationaryDist=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightsParamNames,Policy,n_d,n_a,n_z,N_j,pi_z,Params,simoptions);

AggVars=EvalFnOnAgentDist_AggVars_FHorz_Case1(StationaryDist, Policy, FnsToEvaluate, Params, [], n_d, n_a, n_z,N_j, d_grid, a_grid, z_grid,[],simoptions);
% Aggregate labor is actually exogenous and equal to N=0.3491, but I
% anyway compute it just to show how (can also be used as a double check)
% N=E[z*h*epsilon_j*I_j] (the mewj and 0.94 of IIJ1995 on pg 90 are implicit in the expectation)
[AggVars.N.Mean,0.3491]

FnsToEvaluate2.K = @(aprime,a,z) a; % Aggregate assets K
FnsToEvaluate2.N = @(aprime,a,z,I_j,h,epsilon_j) I_j*h*epsilon_j*z; % Aggregate effective labour supply (in efficiency units), I1998 calls this N
FnsToEvaluate2.C = @(aprime,a,z,r,tau_u, tau_s,h,zeta,epsilon_j,I_j,alpha,delta, A,SSdivw, Tr_beq,workinglifeincome,g,agej,MedicalShock,LumpSum) ImrohorogluImrohorogluJoines1995_ConsumptionFn(aprime,a,z, r,tau_u, tau_s,h,zeta,epsilon_j,I_j,alpha,delta, A,SSdivw, Tr_beq,workinglifeincome,g,agej,MedicalShock,LumpSum); % Consumption
FnsToEvaluate2.Income = @(aprime,a,z,r,h,zeta, epsilon_j,I_j,alpha,delta, A,SSdivw, Tr_beq,workinglifeincome,g,agej) ImrohorogluImrohorogluJoines1995_IncomeFn(aprime,a,z,r,h,zeta, epsilon_j,I_j,alpha,delta, A,SSdivw, Tr_beq,workinglifeincome,g,agej); % Income
FnsToEvaluate2.Utility = @(aprime,a,z,r,tau_u, tau_s,gamma, h,zeta, epsilon_j,I_j,alpha,delta, A,SSdivw, Tr_beq, MedicalShock,workinglifeincome,g,agej,LumpSum) ImrohorogluImrohorogluJoines1995_ReturnFn(aprime,a,z,r,tau_u, tau_s,gamma, h,zeta, epsilon_j,I_j,alpha,delta, A,SSdivw, Tr_beq, MedicalShock,workinglifeincome,g,agej,LumpSum);

AggVars2=EvalFnOnAgentDist_AggVars_FHorz_Case1(StationaryDist, Policy, FnsToEvaluate2, Params, [], n_d, n_a, n_z,N_j, d_grid, a_grid, z_grid,[],simoptions);

K=AggVars2.K.Mean;
N=AggVars2.N.Mean;
C=AggVars2.C.Mean;
Q=Params.A*(K^(Params.alpha))*(N^(1-Params.alpha)); % Cobb-Douglas production fn

Params.w=(1-Params.alpha)*(Params.A^(1/(1-Params.alpha)))*((Params.r+Params.delta)/Params.alpha)^(Params.alpha/(Params.alpha-1));

fprintf('Social Security replacement rate: %8.2f (theta in IIJ1995 notation, b in my notation) \n ', Params.b);
fprintf('(Social Security) Tax rate (tau_s): %8.3f \n', Params.tau_s);
fprintf('Wage rate (w): %8.3f \n', Params.w);
fprintf('Return to capital (r): %8.3f \n ', Params.r);
fprintf('Average Consumption (C): %8.3f \n', C);
fprintf('Capital Stock (K): %8.3f \n', K);
fprintf('Average Income (Q): %8.3f \n', Q);
fprintf('Average Utility (value fn): %8.3f \n ', sum(sum(sum(V.*StationaryDist))) );
fprintf('K/Q: %8.3f \n ', K/Q);

%% Some life-cycle profiles for income, consumption, and assets
LifeCycleProfiles=LifeCycleProfiles_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate2,[],Params,0,n_a,n_z,N_j,0,a_grid,z_grid,simoptions);
% (this creates much more than just the 'age conditional mean' profiles that we use here)
% (Note: when productivity growth is non-zero then you would need to correct some of these)
% I have assumed income includes capital income, unemployment benefits and
% the transfer of accidental bequests (in addition to earnings and social
% security benefits). Note that IIJ1995 do define a concept of 'disposable
% income' which is after-tax earnings plus unemployment benefits plus
% social security benefits.

%% To create Figures 6 and 7 you would also need the value of assets on the grid
ValuesOnGrid=EvalFnOnAgentDist_ValuesOnGrid_FHorz_Case1(Policy, FnsToEvaluate2, Params, [], n_d, n_a, n_z, N_j, d_grid, a_grid, z_grid,[],simoptions);

%% Welfare benefits (Table 3 reports these, the explanation of calculation is on pg 96-97 of IIJ1995)
UtilityOnGrid=ValuesOnGrid.Utility(:,:,:);
discountongrid=shiftdim(cumprod(Params.beta*Params.sj),-1);
AgeConditionalStationaryDist=StationaryDist./sum(sum(StationaryDist,1),2);
Omega1=sum(sum(sum(discountongrid.*AgeConditionalStationaryDist.*UtilityOnGrid)));

% Start by calculating Q0
GEPriceParamNames_Omega0={'r','Tr_beq'}; % 'tau_u' & 'tau_s' have been removed
GeneralEqmEqns_Omega0.Tr=GeneralEqmEqns.Bequests; % Just the two general eqm eqns
GeneralEqmEqns_Omega0.capitalmarket=GeneralEqmEqns.capitalmarket;
% Start by calculating Q0
Params0=Params;
Params0.b=0;
Params0.SSdivw=Params0.b*sum(Params0.h*Params0.epsilon_j(1:(Params0.Jr-1)))/(Params0.Jr-1);
Params0.tau_s=0;
[p_eqm0,p_eqm_index, GeneralEqmEqnsValues]=HeteroAgentStationaryEqm_Case1_FHorz(jequaloneDist,AgeWeightsParamNames,n_d, n_a, n_z, N_j, 0, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns_Omega0, Params0, DiscountFactorParamNames, [], [], [], GEPriceParamNames_Omega0, heteroagentoptions, simoptions,vfoptions);
Params0.r=p_eqm0.r;
Params0.Tr_beq=p_eqm0.Tr_beq;
[V0, Policy0]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j, d_grid, a_grid, z_grid, pi_z, ReturnFn, Params0, DiscountFactorParamNames, [],vfoptions);
StationaryDist0=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightsParamNames,Policy0,n_d,n_a,n_z,N_j,pi_z,Params0,simoptions);
AggVars0=EvalFnOnAgentDist_AggVars_FHorz_Case1(StationaryDist0, Policy0, FnsToEvaluate2, Params0, [], n_d, n_a, n_z,N_j, d_grid, a_grid, z_grid,[],simoptions);
K0=AggVars0.K.Mean;
N0=AggVars0.N.Mean;
Q0=Params0.A*(K0^(Params0.alpha))*(N0^(1-Params0.alpha)); % Cobb-Douglas production fn

% Initial guess for Lstar, the LumpSum
Params0.LumpSum=0; % This is actually already its value, but I put it here to emphasise that this is the 'initial guess' for the welfare calculation.

absOmega0minusOmega1=@(GEandLumpSum) IIJ1995_absOmega0minusOmega1(GEandLumpSum,Omega1,Params0,jequaloneDist,AgeWeightsParamNames,n_d, n_a, n_z, N_j, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, FnsToEvaluate2, GeneralEqmEqns_Omega0, DiscountFactorParamNames, GEPriceParamNames_Omega0,heteroagentoptions,simoptions,vfoptions);
% Find the Lstar that minimizes this
minoptions = optimset('TolX',10^(-6),'TolFun',10^(-6)); % Note that 10^(-6) for TolX is an accuracy to more digits than are reported (for kappa in Table 3)
[GEandLstar,minval]=fminsearch(absOmega0minusOmega1,[Params0.r, Params0.Tr_beq,Params0.LumpSum],minoptions);
    
Lstar=GEandLstar(3);
% IIJ1995 report kappa
kappa=Lstar/Q0;

