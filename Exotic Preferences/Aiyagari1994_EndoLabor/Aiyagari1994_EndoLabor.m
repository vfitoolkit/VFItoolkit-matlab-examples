% Example based on Aiyagari (1994) with Endogenous Labor
%
% Add the decision variable, l, which is the labor supply.
% Involves:
% 1. Creating l (n_l, l_grid, relevant parameters for utility function)
% 2. Modifying the return function
% 3. Modifying the 'FnsToEvaluate' to include an expression for L, and to
% include l in the arguments of the existing FnsToEvaluate.
%
% Note that the codes contain a variable called nonSeperableUtility which
% can be used to switch between seperable and non-seperable preferences,
% both of which are described in the docs. (Only difference is in the
% return fn.) 
% [Note that when nonSeperableUtility=1, then chi has a
% different interpretation, gamma_c is used as gamma, and gamma_l is
% ignored.]

% The changes can be found on lines:
% 29,83. n_l added
% 38-40, 44: parameters
% 75,78. l_grid
% 92-4, 103-5. FnsToEvaluate & GeneralEqmEqns
% 115-6, Return fn
% 121,125: GEPriceParamNames
% 137, 149, 156, 158, 160

%% Set some basic variables

n_l=51;
n_k=2^9;
n_z=21;
n_p=0; % Normally you will want n_p=0, setting a non-zero value here activates the use of a grid on prices.

%Parameters
Params.beta=0.96; %Model period is one-sixth of a year
Params.alpha=0.36;
Params.delta=0.08;
Params.gamma_c=3; % With non-seperable preferences this is used as gamma
Params.gamma_l=2; % With non-seperable preferences this is not used
Params.chi=0.9;
Params.sigma=0.2;
Params.rho=0.6;

Params.nonSeperableUtility=0; % 0 means seperable, 1 means non-seperable utility (difference in codes is inside return function)

Params.tauchen_q=3; %Footnote 33 of Aiyagari(1993WP, pg 25) implicitly says that he uses q=3

% Params has been created as a structure. You can create the individual
% parameters from the structure by running the following command
CreateIndividualParams(Params)


%% Set up the exogenous shock process
% Create markov process for the exogenous labour productivity, l.
Tauchen_q=3; % Footnote 33 of Aiyagari(1993WP, pg 25) implicitly says that he uses q=3
[z_grid,pi_z]=discretizeAR1_Tauchen(0,Params.rho,sqrt((1-Params.rho^2)*Params.sigma^2),n_z,Tauchen_q);
% Note: sigma is standard deviations of s, input needs to be standard deviation of the innovations
% Because s is AR(1), the variance of the innovations is (1-rho^2)*sigma^2
z_grid=exp(z_grid);

% Note: In the exogenous labor model we normalize E[z]=1, but we won't do it for the endogenous labor model.


%% Grids

%In the absence of idiosyncratic risk, the steady state equilibrium is given by
r_ss=1/Params.beta-1;
K_ss=((r_ss+Params.delta)/Params.alpha)^(1/(Params.alpha-1)); %The steady state capital in the absence of aggregate uncertainty.

% Set grid for asset holdings
k_grid=15*K_ss*(linspace(0,1,n_k).^3)'; % The ^3 means most points are near zero, which is where the derivative of the value fn changes most.

% Set grid for endogenous labor
l_grid=linspace(0,1,n_l)';

%Bring model into the notational conventions used by the toolkit
d_grid=l_grid; %There is no d variable
a_grid=k_grid;
%pi_z;
%z_grid

n_d=n_l;
n_a=n_k;
%n_z

% Create functions to be evaluated
FnsToEvaluate.K = @(d,aprime,a,z) a; % We just want the aggregate assets (which is this periods state)
FnsToEvaluate.N = @(d,aprime,a,z) d*z; % We just want the aggregate effective labor supply

% Now define the functions for the General Equilibrium conditions
    % Should be written as LHS of general eqm eqn minus RHS, so that the closer the value given by the function is to 
    % zero, the closer the general eqm condition is to holding.
GeneralEqmEqns.CapitalMarket = @(r,K,N,alpha,delta) r-(alpha*(K^(alpha-1))*(N^(1-alpha))-delta); %The requirement that the interest rate equals the marginal product of capital
GeneralEqmEqns.LaborMarket = @(w,K,N,alpha) w-(1-alpha)*(K^alpha)*(N^(-alpha)); %The requirement that the wage equals the marginal product of labor
% Inputs can be any parameter, price, or aggregate of the FnsToEvaluate

%%
DiscountFactorParamNames={'beta'};

ReturnFn=@(d_val,aprime_val, a_val, z_val,gamma_c,gamma_l,chi,r,w,nonSeperableUtility) Aiyagari1994_EndoLabor_ReturnFn(d_val,aprime_val, a_val, z_val,gamma_c,gamma_l,chi,r,w,nonSeperableUtility);

%%

%Use the toolkit to find the equilibrium price index
GEPriceParamNames={'r','w'};
%Set initial guesses for the general eqm variables
Params.r=0.04;
Params.w=(1-Params.alpha)*((Params.r+Params.delta)/Params.alpha)^(Params.alpha/(Params.alpha-1));

%% Calculate the general equilibrium 

% % We will just use the default options for vfoptions and simoptions
vfoptions=struct();
simoptions=struct();
disp('Calculating price vector corresponding to the stationary eqm')
heteroagentoptions.verbose=1;
[p_eqm,~,GeneralEqmCondn]=HeteroAgentStationaryEqm_Case1(n_d, n_a, n_z, n_p, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, [], [], [], GEPriceParamNames,heteroagentoptions, simoptions, vfoptions);

p_eqm

%% Now that we have the GE, let's calculate a bunch of related objects
disp('Calculating various equilibrium objects')
Params.r=p_eqm.r;
Params.w=p_eqm.w;
[~,Policy]=ValueFnIter_Case1(n_d,n_a,n_z,d_grid,a_grid,z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, [], vfoptions);

% PolicyValues=PolicyInd2Val_Case1(Policy,n_d,n_a,n_z,d_grid,a_grid);

StationaryDist=StationaryDist_Case1(Policy,n_d,n_a,n_z,pi_z, simoptions);

AggVars=EvalFnOnAgentDist_AggVars_Case1(StationaryDist, Policy, FnsToEvaluate,Params, [],n_d, n_a, n_z, d_grid, a_grid,z_grid);

% Calculate savings rate:
% We know production is Y=K^{\alpha}L^{1-\alpha}.
% In equilibrium K is constant, so aggregate savings is just depreciation, which equals delta*K. 
% The agg savings rate is thus delta*K/Y.
% So agg savings rate is given by s=delta*K/(K^{\alpha}L^{1-\alpha})=delta*(K/L)^{1-\alpha}
aggsavingsrate=Params.delta*(AggVars.K.Mean/AggVars.N.Mean)^(1-Params.alpha);

% Calculate Lorenz curves, Gini coefficients, and Pareto tail coefficients
FnsToEvaluate_Ineq.Earnings = @(d,aprime,a,z,w) w*d*z;
FnsToEvaluate_Ineq.Income = @(d,aprime,a,z,r,w) w*d*z+(1+r)*a;
FnsToEvaluate_Ineq.Wealth = @(d,aprime,a,z) a;
LorenzCurves=EvalFnOnAgentDist_LorenzCurve_Case1(StationaryDist, Policy, FnsToEvaluate_Ineq, Params,[], n_d, n_a, n_z, d_grid, a_grid, z_grid);

% 3.5 The Distributions of Earnings and Wealth
%  Gini for Earnings
EarningsGini=Gini_from_LorenzCurve(LorenzCurves.Earnings);
IncomeGini=Gini_from_LorenzCurve(LorenzCurves.Income);
WealthGini=Gini_from_LorenzCurve(LorenzCurves.Wealth);

% Calculate inverted Pareto coeff, b, from the top income shares as b=1/[log(S1%/S0.1%)/log(10)] (forgammala taken from Excel download of WTID database)
% No longer used: Calculate Pareto coeff from Gini as alpha=(1+1/G)/2; ( http://en.wikipedia.org/wiki/Pareto_distribution#Lorenz_curve_and_Gini_coefficient)
% Recalculte Lorenz curves, now with 1000 points
LorenzCurves=EvalFnOnAgentDist_LorenzCurve_Case1(StationaryDist, Policy, FnsToEvaluate_Ineq, Params,[], n_d, n_a, n_z, d_grid, a_grid, z_grid, [],1000);
EarningsParetoCoeff=1/((log(LorenzCurves.Earnings(990))/log(LorenzCurves.Earnings(999)))/log(10));
IncomeParetoCoeff=1/((log(LorenzCurves.Income(990))/log(LorenzCurves.Income(999)))/log(10));
WealthParetoCoeff=1/((log(LorenzCurves.Wealth(990))/log(LorenzCurves.Wealth(999)))/log(10));

%% Display some output about the solution

%plot(cumsum(sum(StationaryDist,2))) %Plot the asset cdf

fprintf('For parameter values sigma=%.2f, gamma_c=%.2f, gamma_l=%.2f, rho=%.2f \n', [Params.sigma,Params.gamma_c,Params.gamma_l,Params.rho])

fprintf('The equilibrium value of the interest rate is r=%.4f \n', p_eqm.r*100)
fprintf('The equilibrium value of the aggregate savings rate is s=%.4f \n', aggsavingsrate)
%fprintf('Time required to find the eqm was %.4f seconds \n',findeqmtime)

