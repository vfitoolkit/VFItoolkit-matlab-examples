% Example of computing a general eqm transition path for the model of Aiyagari (1994).
%
% These codes set up and solve the Aiyagari (1994) model for a given
% parametrization. They then show how to solve for the general equilibrium
% transition path in reposonse to a 'surprise' one off change in the
% parameter alpha (the capital share in production) from 0.36 to 0.40
%
% Transition Path commands require a GPU.

%% Set some basic variables

n_k=512;
n_z=21; 

%Parameters
Params.beta=0.96;
Params.alpha=0.36; % This is actually redundant as declare this below when looking at initial and final eqm
Params.delta=0.08;
Params.mu=3;
Params.sigma=0.2;
Params.rho=0.6;

%% Since this example is intended to show working of transtion paths, make it verbose (print output)
transpathoptions.verbose=1;

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

%%
% Create functions to be evaluated
FnsToEvaluate.K = @(aprime,a,s) a; % We just want the aggregate assets (which is this periods state)

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

%% Compute the initial general equilibrium
Params.alpha=0.36;

% Solve for the stationary general equilbirium
vfoptions=struct(); % Use default options for solving the value function (and policy fn)
vfoptions.gridinterplayer=0; % 0=pure discretization (deafult),1=linear interpolation
vfoptions.ngridinterp=15;
simoptions=struct(); % Use default options for solving for stationary distribution
simoptions.gridinterplayer=vfoptions.gridinterplayer;
simoptions.ngridinterp=vfoptions.ngridinterp;
heteroagentoptions.verbose=1; % verbose means that you want it to give you feedback on what is going on

fprintf('Calculating price vector corresponding to the stationary general eqm \n')
[p_eqm_init,GeneralEqmCondn_init]=HeteroAgentStationaryEqm_Case1(n_d, n_a, n_z, 0, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, [], [], [], GEPriceParamNames,heteroagentoptions, simoptions, vfoptions);

p_eqm_init % The equilibrium values of the GE prices
GeneralEqmCondn_init % The value of the general equilibrium equation, at the the general eqm that we found
% Note: GeneralEqmCondn_init should be essentially zero (that is what it means to be a general eqm)

% For the transition path we will need the initial agents distribution
Params.r=p_eqm_init.r;
[~,Policy_init]=ValueFnIter_Case1(n_d,n_a,n_z,d_grid,a_grid,z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, [], vfoptions);

StationaryDist_init=StationaryDist_Case1(Policy_init,n_d,n_a,n_z,pi_z,simoptions);

% Following line is just a check
AggVars_init=EvalFnOnAgentDist_AggVars_Case1(StationaryDist_init, Policy_init, FnsToEvaluate, Params, [], n_d, n_a, n_z, d_grid, a_grid, z_grid,simoptions);

%% Compute the final general equilbrium
Params.r=0.038; % Initial guess
Params.alpha=0.4;

% Note: if the change in parameters affected pi_z this would need to be recalculated here.

disp('Calculating price vector corresponding to the final stationary eqm')
[p_eqm_final,GeneralEqmCondn_final]=HeteroAgentStationaryEqm_Case1(n_d, n_a, n_z, 0, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, [], [], [], GEPriceParamNames,heteroagentoptions, simoptions, vfoptions);

p_eqm_final % The equilibrium values of the GE prices
% Note: GeneralEqmCondn_final will be essentially zero, it is the value of the general equilibrium equation

% For the transition path we will need the final value function
Params.r=p_eqm_final.r;
[V_final,Policy_final]=ValueFnIter_Case1(n_d,n_a,n_z,d_grid,a_grid,z_grid, pi_z, ReturnFn,Params, DiscountFactorParamNames,[],vfoptions);

StationaryDist_final=StationaryDist_Case1(Policy_final,n_d,n_a,n_z,pi_z);
AggVars_final=EvalFnOnAgentDist_AggVars_Case1(StationaryDist_final, Policy_final, FnsToEvaluate, Params, [], n_d, n_a, n_z, d_grid, a_grid, z_grid,simoptions);

% surf(k_grid*ones(1,n_s),ones(n_a,1)*s_grid',V_final)

%% Compute the transition path
% For this we need the following extra objects: PricePathOld, PriceParamNames, ParamPath, ParamPathNames, T, V_final, StationaryDist_init
% (already calculated V_final & StationaryDist_init above)

% Number of time periods to allow for the transition (if you set T too low
% it will cause problems, too high just means run-time will be longer).
T=150

% We want to look at a one off unanticipated change of alpha [not necessarily a sensible 'reform', but shows how we can do these transition paths]
% ParamPath & PathParamNames are thus given by
ParamPath.alpha=0.4*ones(T,1); % For each parameter that changes value, ParamPath is matrix of size T-by-1
% (the way ParamPath is set is designed to allow for a series of changes in the parameters)

% We need to give an initial guess for the price path on interest rates
% (this is deliberately not a good guess, so you can see that the transition path can be found)
PricePath0.r=[linspace(p_eqm_init.r, p_eqm_final.r, floor(T/2))'; p_eqm_final.r*ones(T-floor(T/2),1)]; % For each price, PricePath0 is matrix of size T-by-1

% General equilibrium conditions (for the transition path)
TransPathGeneralEqmEqns.CapitalMarket = @(r,K,alpha,delta,Expectation_l) r-(alpha*(K^(alpha-1))*(Expectation_l^(1-alpha))-delta);
% Note: For this model the transition path has the same general equilibrium conditions as the stationary equilibrium, but this will not always be true for more complex models.

transpathoptions.GEnewprice=3;
% Need to explain to transpathoptions how to use the GeneralEqmEqns to
% update the general eqm transition prices (in PricePath).
transpathoptions.GEnewprice3.howtoupdate=... % a row is: GEcondn, price, add, factor
    {'CaptialMarket','r',0,0.2}; % CaptialMarket GE condition will be positive if r is too big, so subtract
% Note: the update is essentially new_price=price+factor*add*GEcondn_value-factor*(1-add)*GEcondn_value
% Notice that this adds factor*GEcondn_value when add=1 and subtracts it what add=0
% A small 'factor' will make the convergence to solution take longer, but too large a value will make it 
% unstable (fail to converge). Technically this is the damping factor in a shooting algorithm.


% Now just run the TransitionPath_Case1 command (all of the other inputs
% are things we had already had to define to be able to solve for the
% initial and final equilibria)
transpathoptions.weightscheme=1;
transpathoptions.verbose=1;

[PricePath,GeneralEqmCondnPath]=TransitionPath_Case1(PricePath0, ParamPath, T, V_final, StationaryDist_init, n_d, n_a, n_z, pi_z, d_grid,a_grid,z_grid, ReturnFn, FnsToEvaluate, TransPathGeneralEqmEqns, Params, DiscountFactorParamNames, transpathoptions,vfoptions,simoptions);

figure(1)
plot(0:1:T, [p_eqm_init.r,PricePath.r])
title('interest rate path for transtion')

%% Look at results
[VPath,PolicyPath]=ValueFnOnTransPath_Case1(PricePath, ParamPath, T, V_final, Policy_final, Params, n_d, n_a, n_z, pi_z, d_grid, a_grid,z_grid, DiscountFactorParamNames, ReturnFn, transpathoptions, vfoptions);

AgentDistPath=AgentDistOnTransPath_Case1(StationaryDist_init, PolicyPath,n_d,n_a,n_z,pi_z,T,simoptions);

AggVarsPath=EvalFnOnTransPath_AggVars_Case1(FnsToEvaluate,AgentDistPath,PolicyPath,PricePath,ParamPath, Params, T, n_d, n_a, n_z, d_grid, a_grid,z_grid,simoptions);

figure(2)
plot(0:1:T, [AggVars_init.K.Mean, AggVarsPath.K.Mean])
title('path of aggregate physical capital for transtion')


