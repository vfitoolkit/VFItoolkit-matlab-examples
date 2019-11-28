% Example based on Aiyagari (1994).

% These codes set up and solve the Aiyagari (1994) model for a given
% parametrization. After solving the model they then show how some of the
% vfitoolkit commands to easily calculate things like the Gini coefficient
% for income, and how to plot the distribution of asset holdings.

Parallel=2 % 2 for GPU, 1 for parallel CPU, 0 for single CPU.

%% Set some basic variables

n_k=2^7;%2^9;
n_s=11; %21;
n_p=0; % Normally you will want n_p=0, setting a non-zero value here activates the use of a grid on prices.

%Parameters
Params.beta=0.96; %Model period is one-sixth of a year
Params.alpha=0.36;
Params.delta=0.08;
Params.mu=3;
Params.sigma=0.2;
Params.rho=0.6;

Params.q=3; %Footnote 33 of Aiyagari(1993WP, pg 25) implicitly says that he uses q=3

% Params has been created as a structure. You can create the individual
% parameters from the structure by running the following command
CreateIndividualParams(Params)

%% Some Toolkit options (most of these are anyway just being set to toolkit defaults)
tauchenoptions.parallel=Parallel;

mcmomentsoptions.T=10^4;
mcmomentsoptions.Tolerance=10^(-9);
mcmomentsoptions.parallel=tauchenoptions.parallel;

vfoptions.lowmemory=0;
vfoptions.parallel=Parallel;

simoptions.burnin=10^4;
simoptions.simperiods=10^5; % For an accurate solution you will either need simperiod=10^5 and iterate=1, or simperiod=10^6 (iterate=0).
simoptions.iterate=1;
simoptions.parallel=Parallel; %Use GPU

heteroagentoptions.verbose=1;

%% Set up the exogenous shock process
%Create markov process for the exogenous labour productivity, l.
% q=3; %Footnote 33 of Aiyagari(1993WP, pg 25) implicitly says that he uses q=3
[s_grid, pi_s]=TauchenMethod(0,(Params.sigma^2)*(1-Params.rho^2),Params.rho,n_s,Params.q,tauchenoptions); %[states, transmatrix]=TauchenMethod_Param(mew,sigmasq,rho,znum,q,Parallel,Verbose), transmatix is (z,zprime)

%[s_grid, pi_s]=TauchenMethod_Param(0,(sigma^2)*(1-rho^2),rho,7,3); %This is the process for ln(l). Aiyagari uses 7 states, 
[s_mean,s_variance,s_corr,~]=MarkovChainMoments(s_grid,pi_s,mcmomentsoptions);
s_grid=exp(s_grid);
%Get some info on the markov process
[Expectation_l,~,~,~]=MarkovChainMoments(s_grid,pi_s,mcmomentsoptions); %Since l is exogenous, this will be it's eqm value 
%Note: Aiyagari (1994) actually then normalizes l by dividing it by
%Expectation_l (so that the resulting process has expectaion equal to 1
%(see Aiyagari (1993WP), footnote 33 pg 25-26).
%The following three lines do just this.
s_grid=s_grid./Expectation_l;
[Expectation_l,~,~,~]=MarkovChainMoments(s_grid,pi_s,mcmomentsoptions);

%% Grids

%In the absence of idiosyncratic risk, the steady state equilibrium is given by
r_ss=1/Params.beta-1;
K_ss=((r_ss+Params.delta)/Params.alpha)^(1/(Params.alpha-1)); %The steady state capital in the absence of aggregate uncertainty.

% Set grid for asset holdings
%Aiyagari uses 25 points, but with a piecewise-linear approx. see Aiyagari (1993WP, pg 28).
%His grid is not directly on k, but implicitly his k grid runs from zero up
%to k_max, where k_max is given by f(k_max,1)=delta*k_max
%k_max=delta^(1/(alpha-1));
%Doing this k_max is slightly less than 10*K_ss. But if I use this as the upper limit on
%the grid the results are wrong (in that increasing this to 15 or 20*K_ss
%gives different results). It might be that Aiyagari gets away with this
%due to his use of piecewise-linear approximation (ie. that policy fn is
%almost linear in this region anyway).
nk1=floor(n_k/3); nk2=floor(n_k/3); nk3=n_k-nk1-nk2;
k_grid=sort([linspace(0,K_ss,nk1),linspace(K_ss+0.0001,3*K_ss,nk2),linspace(3*K_ss+0.0001,15*K_ss,nk3)]');

%Bring model into the notational conventions used by the toolkit
d_grid=0; %There is no d variable
a_grid=k_grid;
%pi_s;
%s_grid

n_d=0;
n_a=n_k;
%n_s

%Create descriptions of SS values as functions of d_grid, a_grid, s_grid &
%pi_s (used to calculate the integral across the SS dist fn of whatever
%functions you define here)
FnsToEvaluateParamNames(1).Names={};
FnsToEvaluateFn_1 = @(aprime_val,a_val,s_val) a_val; %We just want the aggregate assets (which is this periods state)
FnsToEvaluate={FnsToEvaluateFn_1};

%Now define the functions for the General Equilibrium conditions
    %Should be written as LHS of general eqm eqn minus RHS, so that 
    %the closer the value given by the function is to zero, the closer 
    %the general eqm condition is to holding.
%Note: length(AggVars) is as for SSvaluesFn and length(p) is number_p_vars
GeneralEqmEqnParamNames(1).Names={'alpha','delta'};
GeneralEqmEqn_1 = @(AggVars,p,alpha,delta) p-(alpha*(AggVars^(alpha-1))*(Expectation_l^(1-alpha))-delta); %The requirement that the interest rate corresponds to the agg capital level
GeneralEqmEqns={GeneralEqmEqn_1};

disp('sizes')
n_a
n_s
n_p

%%
DiscountFactorParamNames={'beta'};

ReturnFn=@(aprime_val, a_val, s_val,alpha,delta,mu,r) Aiyagari1994_ReturnFn(aprime_val, a_val, s_val,alpha,delta,mu,r);
ReturnFnParamNames={'alpha','delta','mu','r'}; %It is important that these are in same order as they appear in 'Aiyagari1994_ReturnFn'

%%

%Use the toolkit to find the equilibrium price index
GEPriceParamNames={'r'};
%Set initial value for interest rates (Aiyagari proves that with idiosyncratic
%uncertainty, the eqm interest rate is limited above by it's steady state value
%without idiosyncratic uncertainty, that is that r<r_ss).
Params.r=0.04;

%%
V0=ones(n_a,n_s,'gpuArray'); %(a,s)

disp('Calculating price vector corresponding to the stationary eqm')
[p_eqm,~,GeneralEqmCondn]=HeteroAgentStationaryEqm_Case1(V0, n_d, n_a, n_s, n_p, pi_s, d_grid, a_grid, s_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, ReturnFnParamNames, FnsToEvaluateParamNames, GeneralEqmEqnParamNames, GEPriceParamNames,heteroagentoptions, simoptions, vfoptions);

p_eqm

%% Now that we have the GE, let's calculate a bunch of related objects
% Equilibrium wage
Params.w=(1-Params.alpha)*((p_eqm+Params.delta)/Params.alpha)^(Params.alpha/(Params.alpha-1));

disp('Calculating various equilibrium objects')
Params.r=p_eqm;
[~,Policy]=ValueFnIter_Case1(V0, n_d,n_a,n_s,d_grid,a_grid,s_grid, pi_s, ReturnFn, Params, DiscountFactorParamNames, ReturnFnParamNames, vfoptions);

% PolicyValues=PolicyInd2Val_Case1(Policy,n_d,n_a,n_s,d_grid,a_grid, Parallel);

StationaryDist=StationaryDist_Case1(Policy,n_d,n_a,n_s,pi_s, simoptions);

AggVars=EvalFnOnAgentDist_AggVars_Case1(StationaryDist, Policy, FnsToEvaluate,Params, FnsToEvaluateParamNames,n_d, n_a, n_s, d_grid, a_grid,s_grid,Parallel)

% save ./SavedOutput/Aiyagari1994SSObjects.mat p_eqm Policy StationaryDist

% Calculate savings rate:
% We know production is Y=K^{\alpha}L^{1-\alpha}, and that L=1
% (exogeneous). Thus Y=K^{\alpha}.
% In equilibrium K is constant, so aggregate savings is just depreciation, which
% equals delta*K. The agg savings rate is thus delta*K/Y.
% So agg savings rate is given by s=delta*K/(K^{\alpha})=delta*K^{1-\alpha}
aggsavingsrate=Params.delta*AggVars^(1-Params.alpha);

% Calculate Lorenz curves, Gini coefficients, and Pareto tail coefficients
FnsToEvaluateParamNames(1).Names={'w'};
FnsToEvaluate_Earnings = @(aprime_val,a_val,s_val,param) w*s_val;
FnsToEvaluateParamNames(2).Names={'r','w'};
FnsToEvaluate_Income = @(aprime_val,a_val,s_val,r,w) w*s_val+(1+r)*a_val;
FnsToEvaluateParamNames(3).Names={};
FnsToEvaluate_Wealth = @(aprime_val,a_val,s_val) a_val;
FnsToEvaluateFnIneq={FnsToEvaluate_Earnings, FnsToEvaluate_Income, FnsToEvaluate_Wealth};
StationaryDist_LorenzCurves=EvalFnOnAgentDist_LorenzCurve_Case1(StationaryDist, Policy, FnsToEvaluateFnIneq, Params,FnsToEvaluateParamNames, n_d, n_a, n_s, d_grid, a_grid, s_grid, Parallel);

% 3.5 The Distributions of Earnings and Wealth
%  Gini for Earnings
EarningsGini=Gini_from_LorenzCurve(StationaryDist_LorenzCurves(1,:));
IncomeGini=Gini_from_LorenzCurve(StationaryDist_LorenzCurves(2,:));
WealthGini=Gini_from_LorenzCurve(StationaryDist_LorenzCurves(3,:));

% Calculate inverted Pareto coeff, b, from the top income shares as b=1/[log(S1%/S0.1%)/log(10)] (formula taken from Excel download of WTID database)
% No longer used: Calculate Pareto coeff from Gini as alpha=(1+1/G)/2; ( http://en.wikipedia.org/wiki/Pareto_distribution#Lorenz_curve_and_Gini_coefficient)
% Recalculte Lorenz curves, now with 1000 points
StationaryDist_LorenzCurves=EvalFnOnAgentDist_LorenzCurve_Case1(StationaryDist, Policy, FnsToEvaluateFnIneq, Params,FnsToEvaluateParamNames, n_d, n_a, n_s, d_grid, a_grid, s_grid, Parallel,1000);
EarningsParetoCoeff=1/((log(StationaryDist_LorenzCurves(1,990))/log(StationaryDist_LorenzCurves(1,999)))/log(10)); %(1+1/EarningsGini)/2;
IncomeParetoCoeff=1/((log(StationaryDist_LorenzCurves(2,990))/log(StationaryDist_LorenzCurves(2,999)))/log(10)); %(1+1/IncomeGini)/2;
WealthParetoCoeff=1/((log(StationaryDist_LorenzCurves(3,990))/log(StationaryDist_LorenzCurves(3,999)))/log(10)); %(1+1/WealthGini)/2;


%% Display some output about the solution

%plot(cumsum(sum(StationaryDist,2))) %Plot the asset cdf

fprintf('For parameter values sigma=%.2f, mu=%.2f, rho=%.2f \n', [Params.sigma,Params.mu,Params.rho])
fprintf('The table 1 elements are sigma=%.4f, rho=%.4f \n',[sqrt(s_variance), s_corr])

fprintf('The equilibrium value of the interest rate is r=%.4f \n', p_eqm*100)
fprintf('The equilibrium value of the aggregate savings rate is s=%.4f \n', aggsavingsrate)
%fprintf('Time required to find the eqm was %.4f seconds \n',findeqmtime)

