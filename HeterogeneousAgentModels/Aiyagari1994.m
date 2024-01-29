% Example based on Aiyagari (1994).
%
% These codes set up and solve the Aiyagari (1994) model for a given
% parametrization. After solving the model they then show how some of the
% vfitoolkit commands to easily calculate things like the Gini coefficient
% for income, and how to plot the distribution of asset holdings.
%
% VFI Toolkit automatically detects hardware (GPU? Number of CPUs?) and
% sets defaults accordingly. It will run without a GPU, but slowly. It is
% indended for use with GPU.

%% Set some basic variables

% VFI Toolkit thinks of there as being:
% k: an endogenous state variable (assets)
% z: an exogenous state variable (exogenous labor supply)

% Size of the grids
n_k=2^9;%2^9;
n_z=11; %21;

% Parameters
Params.beta=0.96; %Model period is one-sixth of a year
Params.alpha=0.36;
Params.delta=0.08;
Params.mu=3;
Params.sigma=0.2;
Params.rho=0.6;

%% Set up the exogenous shock process
% Create markov process for the exogenous labour productivity, l.
Tauchen_q=3; % Footnote 33 of Aiyagari(1993WP, pg 25) implicitly says that he uses q=3
[z_grid,pi_z]=discretizeAR1_Tauchen(0,Params.rho,sqrt((1-Params.rho^2)*Params.sigma^2),n_z,Tauchen_q);
% Note: sigma is standard deviations of s, input needs to be standard deviation of the innovations
% Because s is AR(1), the variance of the innovations is (1-rho^2)*sigma^2

[z_mean,z_variance,z_corr,~]=MarkovChainMoments(z_grid,pi_z);
z_grid=exp(z_grid);
% Get some info on the markov process
[Expectation_l,~,~,~]=MarkovChainMoments(z_grid,pi_z); %Since l is exogenous, this will be it's eqm value 
% Note: Aiyagari (1994) actually then normalizes l by dividing it by Expectation_l (so that the resulting process has expectation equal to 1)
z_grid=z_grid./Expectation_l;
[Expectation_l,~,~,~]=MarkovChainMoments(z_grid,pi_z);
% If you look at Expectation_l you will see it is now equal to 1
Params.Expectation_l=Expectation_l;

%% Grids

% In the absence of idiosyncratic risk, the steady state equilibrium is given by
r_ss=1/Params.beta-1;
K_ss=((r_ss+Params.delta)/Params.alpha)^(1/(Params.alpha-1)); %The steady state capital in the absence of aggregate uncertainty.

% Set grid for asset holdings
k_grid=10*K_ss*(linspace(0,1,n_k).^3)'; % linspace ^3 puts more points near zero, where the curvature of value and policy functions is higher and where model spends more time

% Bring model into the notational conventions used by the toolkit
d_grid=0; %There is no d variable
a_grid=k_grid;
% pi_z;
% z_grid

n_d=0;
n_a=n_k;
% n_z

% Create functions to be evaluated
FnsToEvaluate.K = @(aprime,a,s) a; %We just want the aggregate assets (which is this periods state)

% Now define the functions for the General Equilibrium conditions
    % Should be written as LHS of general eqm eqn minus RHS, so that the closer the value given by the function is to 
    % zero, the closer the general eqm condition is to holding.
GeneralEqmEqns.CapitalMarket = @(r,K,alpha,delta,Expectation_l) r-(alpha*(K^(alpha-1))*(Expectation_l^(1-alpha))-delta); %The requirement that the interest rate corresponds to the agg capital level
% Inputs can be any parameter, price, or aggregate of the FnsToEvaluate

fprintf('Grid sizes are: %i points for assets, and %i points for exogenous shock \n', n_a,n_z)

%%
DiscountFactorParamNames={'beta'};

ReturnFn=@(aprime, a, s, alpha,delta,mu,r) Aiyagari1994_ReturnFn(aprime, a, s,alpha,delta,mu,r);
% The first inputs must be: next period endogenous state, endogenous state, exogenous state. Followed by any parameters

%%

% Use the toolkit to find the equilibrium price index
GEPriceParamNames={'r'};
% Set initial value for interest rates (Aiyagari proves that with idiosyncratic
% uncertainty, the eqm interest rate is limited above by it's steady state value
% without idiosyncratic uncertainty, that is that r<r_ss).
Params.r=0.038;


%%
% Solve for the stationary general equilbirium
vfoptions=struct(); % Use default options for solving the value function (and policy fn)
simoptions=struct(); % Use default options for solving for stationary distribution
heteroagentoptions.verbose=1; % verbose means that you want it to give you feedback on what is going on

fprintf('Calculating price vector corresponding to the stationary general eqm \n')
[p_eqm,~,GeneralEqmCondn]=HeteroAgentStationaryEqm_Case1(n_d, n_a, n_z, 0, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, [], [], [], GEPriceParamNames,heteroagentoptions, simoptions, vfoptions);

p_eqm % The equilibrium values of the GE prices

%% Now that we have the GE, let's calculate a bunch of related objects
Params.r=p_eqm.r; % Put the equilibrium interest rate into Params so we can use it to calculate things based on equilibrium parameters

% Equilibrium wage
Params.w=(1-Params.alpha)*((Params.r+Params.delta)/Params.alpha)^(Params.alpha/(Params.alpha-1));

fprintf('Calculating various equilibrium objects \n')
[V,Policy]=ValueFnIter_Case1(n_d,n_a,n_z,d_grid,a_grid,z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, [], vfoptions);
% V is value function
% Policy is policy function (but as an index of k_grid, not the actual values)

% PolicyValues=PolicyInd2Val_Case1(Policy,n_d,n_a,n_s,d_grid,a_grid); % This will give you the policy in terms of values rather than index

StationaryDist=StationaryDist_Case1(Policy,n_d,n_a,n_z,pi_z, simoptions);

AggVars=EvalFnOnAgentDist_AggVars_Case1(StationaryDist, Policy, FnsToEvaluate,Params, [],n_d, n_a, n_z, d_grid, a_grid,z_grid);

% AggVars contains the aggregate values of the 'FnsToEvaluate' (in this model aggregates are equal to the mean expectation over the agent distribution)
% Currently the only FnsToEvaluate is assets, so we get aggregate capital stock
AggVars.K.Mean

% Calculate savings rate:
% We know production is Y=K^{\alpha}L^{1-\alpha}, and that L=1
% (exogeneous). Thus Y=K^{\alpha}.
% In equilibrium K is constant, so aggregate savings is just depreciation, which
% equals delta*K. The agg savings rate is thus delta*K/Y.
% So agg savings rate is given by s=delta*K/(K^{\alpha})=delta*K^{1-\alpha}
aggsavingsrate=Params.delta*(AggVars.K.Mean)^(1-Params.alpha);

% Calculate Lorenz curves, Gini coefficients, and Pareto tail coefficients
% Declare the function we want to evaluate
FnsToEvaluateIneq.Earnings = @(aprime,a,z,w) w*z;
FnsToEvaluateIneq.Income = @(aprime,a,z,r,w) w*z+(1+r)*a;
FnsToEvaluateIneq.Wealth = @(aprime,a,s) a;

% By default, lorenz curves are calculated with 100 points, we will want to
% use 1000 points, so
simoptions.npoints=1000; % number of points to use for lorenz curves (default is 100) 
AllStats=EvalFnOnAgentDist_AllStats_Case1(StationaryDist, Policy, FnsToEvaluateIneq,Params, [],n_d, n_a, n_z, d_grid, a_grid,z_grid,simoptions);

% 3.5 The Distributions of Earnings and Wealth
%  Gini for Earnings
fprintf('The Gini coeff for Earnings is: %8.4f \n ',AllStats.Earnings.Gini)
fprintf('The Gini coeff for Income is: %8.4f \n ',AllStats.Income.Gini)
fprintf('The Gini coeff for Wealth is: %8.4f \n ',AllStats.Wealth.Gini)

% Calculate inverted Pareto coeff, b, from the top income shares as b=1/[log(S1%/S0.1%)/log(10)] (formula taken from Excel download of WTID database)
% No longer used: Calculate Pareto coeff from Gini as alpha=(1+1/G)/2; ( http://en.wikipedia.org/wiki/Pareto_distribution#Lorenz_curve_and_Gini_coefficient)
% Recalculate Lorenz curves, now with 1000 points
EarningsParetoCoeff=1/((log(AllStats.Earnings.LorenzCurve(990))/log(AllStats.Earnings.LorenzCurve(999)))/log(10)); %(1+1/EarningsGini)/2;
IncomeParetoCoeff=1/((log(AllStats.Income.LorenzCurve(990))/log(AllStats.Income.LorenzCurve(999)))/log(10)); %(1+1/IncomeGini)/2;
WealthParetoCoeff=1/((log(AllStats.Wealth.LorenzCurve(990))/log(AllStats.Wealth.LorenzCurve(999)))/log(10)); %(1+1/WealthGini)/2;


%% Display some output about the solution

% plot(cumsum(sum(StationaryDist,2))) %Plot the asset cdf

fprintf('For parameter values sigma=%.2f, mu=%.2f, rho=%.2f \n', [Params.sigma,Params.mu,Params.rho])
fprintf('The table 1 elements are sigma=%.4f, rho=%.4f \n',[sqrt(z_variance), z_corr])

fprintf('The equilibrium value of the interest rate is r=%.4f \n', p_eqm.r*100)
fprintf('The equilibrium value of the aggregate savings rate is s=%.4f \n', aggsavingsrate)

