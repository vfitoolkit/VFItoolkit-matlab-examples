% Example based on Aiyagari model with Quasi-Hyperbolic discouting.
%
%% Use Quasi-Hyperbolic discounting. I have put all the lines that relate to Quasi-hyperbolic discounting together to simplify reading.
% There are essentially three parts to using Quasi-Hyperbolic discounting.
% 1. Use vfoptions to state that you are using Quasi-Hyperbolic discounting.
% 2. (optional) Use vfoptions to state whether you want 'naive' or 'sophisticated' Quasi-Hyperbolic discounting. (Optional as naive by default)
% 3. Set the appropriate preference parameters.
% 4. Minor adjustment to 'discount factors'.

% 1. Use vfoptions to state that you are using Quasi-Hyperbolic discounting.
% The quasi-hyperbolic preferences are controlled using,
vfoptions.exoticpreferences='QuasiHyperbolic' % Use the quasi-hyperbolic preferences
% To turn off, either don't delare vfoptions.exoticpreferences, or set vfoptions.exoticpreferences='None'
% If you turn this off you would just solve the same model but with standard exponential discounting (parameter beta0 removed from discount factor by code below). 
% This is intended as purely illustrative. A serious comparison of the preference types would require you to recalibrate the model.

% 2. (optional) Use vfoptions to state whether you want 'naive' or 'sophisticated' Quasi-Hyperbolic discounting. (Optional as naive by default)
% Also need to choose which of the two quasi-hyperbolic solutions, naive or sophisticated, to use.
% vfoptions.quasi_hyperbolic='Naive' % This is the default, alternative is 'Sophisticated'.
vfoptions.quasi_hyperbolic='Naive' % This is the default, alternative is 'Sophisticated'.

% WARNING: Sophisticated Quasi-Hyperbolic discounting in infinite horizon appears to
% require insanely large number of grid points for accurate convergence of
% V, so it is imposing a maximum number of iterations and stopping after
% than rather than when V converges. (Note that this is not a problem for
% Quasi-Hyperbolic discounting in finite horizon, only infinite horizon)

% 3. Set the appropriate preference parameters.
Params.beta0=0.85; % The quasi-hyperbolic discounting parameter controlling 'additional' discounting between 'today and tomorrow'
Params.beta=0.96; % The quasi-hyperbolic discounting parameter controlling discounting between any two periods
% Quasi-hyperbolic discounting is sometimes refered to as 'alpha-beta'
% discounting. In this case alpha would refer to beta0*beta which is the discount factor between the current period and next period (and beta is beta :).
% Note that setting beta0=1 would give standard exponential discounting.

% 4. Minor adjustment to 'discount factors'.
% For quasi-hyperbolic preferences the last element of DiscountFactorParamNames must be the 'today-to-tomorrow' additional discount factor.
DiscountFactorParamNames={'beta','beta0'};
% Note that when using Quasi-hyperbolic discounting beta just acts like a
% standard discounting parameter, only beta0 needs to be treated specially.

% If not using quasi-hyperbolic discounting then disable beta0 by setting
if strcmp(vfoptions.exoticpreferences,'None')
    DiscountFactorParamNames={'beta','sj','gdiscount'};
end

% That is all. Every other line of code is unchanged!! Quasi-Hyperbolic discounting really is that easy ;)

%% Set some basic variables

n_k=2^10;
n_z=21;
n_p=0; % Normally you will want n_p=0, setting a non-zero value here activates the use of a grid on prices.

%Parameters
Params.beta=0.96; %Model period is one-sixth of a year
Params.alpha=0.36;
Params.delta=0.08;
Params.gamma=3;
Params.sigma=0.2;
Params.rho=0.6;

% I HAVE CHANGED THIS TO 2
Params.q=2; %Footnote 33 of Aiyagari(1993WP, pg 25) implicitly says that he uses q=3

%% Set up the exogenous shock process
% Create markov process for the exogenous labour productivity, l.
Tauchen_q=3; % Footnote 33 of Aiyagari(1993WP, pg 25) implicitly says that he uses q=3
[z_grid,pi_z]=discretizeAR1_Tauchen(0,Params.rho,sqrt((1-Params.rho^2)*Params.sigma^2),n_z,Tauchen_q);
% Note: sigma is standard deviations of s, input needs to be standard deviation of the innovations
% Because s is AR(1), the variance of the innovations is (1-rho^2)*sigma^2

[z_mean,z_variance,z_corr,~]=MarkovChainMoments(z_grid,pi_z);
z_grid=exp(z_grid);
% Get some info on the markov process
[Expectation_z,~,~,~]=MarkovChainMoments(z_grid,pi_z); %Since l is exogenous, this will be it's eqm value 
% Note: Aiyagari (1994) actually then normalizes l by dividing it by Expectation_z (so that the resulting process has expectation equal to 1)
z_grid=z_grid./Expectation_z;
[Expectation_z,~,~,~]=MarkovChainMoments(z_grid,pi_z);
% If you look at Expectation_z you will see it is now equal to 1
Params.Expectation_z=Expectation_z; % Need to put it in Params so it can be used as part of the general eqm eqn

%% Grids

%In the absence of idiosyncratic risk, the steady state equilibrium is given by
r_ss=1/Params.beta-1;
K_ss=((r_ss+Params.delta)/Params.alpha)^(1/(Params.alpha-1)); %The steady state capital in the absence of aggregate uncertainty.

% Set grid for asset holdings
% I MADE THIS 25 INSTEAD OF 15
k_grid=10*K_ss*(linspace(0,1,n_k).^3)'; % The ^3 means most points are near zero, which is where the derivative of the value fn changes most.

%Bring model into the notational conventions used by the toolkit
d_grid=0; %There is no d variable
a_grid=k_grid;
%pi_z;
%z_grid

n_d=0;
n_a=n_k;
%n_z

% Create functions to be evaluated
FnsToEvaluate.K = @(aprime,a,z) a; %We just want the aggregate assets (which is this periods state)

% Now define the functions for the General Equilibrium conditions
    % Should be written as LHS of general eqm eqn minus RHS, so that the closer the value given by the function is to 
    % zero, the closer the general eqm condition is to holding.
GeneralEqmEqns.CapitalMarket = @(r,K,alpha,delta,Expectation_z) r-(alpha*(K^(alpha-1))*(Expectation_z^(1-alpha))-delta); %The requirement that the interest rate corresponds to the agg capital level
% Inputs can be any parameter, price, or aggregate of the FnsToEvaluate


%%
% DiscountFactorParamNames was set above as part of implementing Quasi-Hyperbolic discounting.
ReturnFn=@(aprime_val, a_val, z_val,alpha,delta,gamma,r) Aiyagari1994_ReturnFn(aprime_val, a_val, z_val,alpha,delta,gamma,r);

%%

%Use the toolkit to find the equilibrium price index
GEPriceParamNames={'r'};
%Set initial value for interest rates (Aiyagari proves that with idiosyncratic
%uncertainty, the eqm interest rate is limited above by it's steady state value
%without idiosyncratic uncertainty, that is that r<r_ss).
Params.r=0.04;

%% Calculate the general equilibrium 

% We will just use the default options for simoptions
simoptions=struct();
disp('Calculating price vector corresponding to the stationary eqm')
heteroagentoptions.verbose=1;
[p_eqm,~,GeneralEqmCondn]=HeteroAgentStationaryEqm_Case1(n_d, n_a, n_z, n_p, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, [], [], [], GEPriceParamNames,heteroagentoptions, simoptions, vfoptions);

p_eqm

%% Now that we have the GE, let's calculate a bunch of related objects
% Equilibrium wage
Params.r=p_eqm.r;
Params.w=(1-Params.alpha)*((Params.r+Params.delta)/Params.alpha)^(Params.alpha/(Params.alpha-1));

disp('Calculating various equilibrium objects')
[~,Policy]=ValueFnIter_Case1(n_d,n_a,n_z,d_grid,a_grid,z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, [], vfoptions);

% PolicyValues=PolicyInd2Val_Case1(Policy,n_d,n_a,n_z,d_grid,a_grid);

simoptions.verbose=1;
StationaryDist=StationaryDist_Case1(Policy,n_d,n_a,n_z,pi_z, simoptions);

figure(2)
plot(cumsum(sum(StationaryDist,2)))

AggVars=EvalFnOnAgentDist_AggVars_Case1(StationaryDist, Policy, FnsToEvaluate,Params, [],n_d, n_a, n_z, d_grid, a_grid,z_grid);

% Calculate savings rate:
% We know production is Y=K^{\alpha}L^{1-\alpha}, and that L=1
% (exogeneous). Thus Y=K^{\alpha}.
% In equilibrium K is constant, so aggregate savings is just depreciation, which
% equals delta*K. The agg savings rate is thus delta*K/Y.
% So agg savings rate is given by s=delta*K/(K^{\alpha})=delta*K^{1-\alpha}
aggsavingsrate=Params.delta*AggVars.K.Mean^(1-Params.alpha);

% Calculate Lorenz curves, Gini coefficients, and Pareto tail coefficients
FnsToEvaluate_Ineq.Earnings = @(aprime,a,z,w) w*z;
FnsToEvaluate_Ineq.Income = @(aprime,a,z,r,w) w*z+(1+r)*a;
FnsToEvaluate_Ineq.Wealth = @(aprime,a,z) a;
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

fprintf('For parameter values sigma=%.2f, gamma=%.2f, rho=%.2f \n', [Params.sigma,Params.gamma,Params.rho])
fprintf('The table 1 elements are sigma=%.4f, rho=%.4f \n',[sqrt(z_variance), z_corr])

fprintf('The equilibrium value of the interest rate is r=%.4f \n', p_eqm.r*100)
fprintf('The equilibrium value of the aggregate savings rate is s=%.4f \n', aggsavingsrate)
%fprintf('Time required to find the eqm was %.4f seconds \n',findeqmtime)

