% Example based on Aiyagari (1994).
%
% This script sets up and solves the Aiyagari (1994) incomplete-markets model
% for a given parameterization. After solving the model, it illustrates how
% to use selected VFI Toolkit commands to compute inequality statistics
% (e.g. Gini coefficients) and to plot the stationary distribution of assets.
%
% The VFI Toolkit automatically detects available hardware (GPU, number of
% CPU cores) and sets defaults accordingly. The code will run without a GPU,
% but substantially more slowly. It is intended primarily for GPU use.

%% Basic setup

% VFI Toolkit notation:
%   a (or k): endogenous state variable (assets)
%   z       : exogenous state variable (labor productivity)

% Grid sizes
n_k = 601;
n_z = 11;

%% Parameters
Params.beta  = 0.96;   % Discount factor, note that the model period is one-sixth of a year
Params.alpha = 0.36;   % Capital elasticity
Params.delta = 0.08;   % Depreciation rate
Params.mu    = 3;      % Coefficient of relative risk aversion in utility
Params.sigma = 0.2;    % Standard deviation of log productivity
Params.rho   = 0.6;    % Autocorrelation of log productivity

%% Exogenous shock process

% Discretize AR(1) labor productivity using Tauchen (1986)
Tauchen_q = 3; % As in footnote 33 of Aiyagari (1993 WP)

% Tauchen requires the std. dev. of innovations.
% Since z_t is AR(1), Var(epsilon) = (1 - rho^2) * sigma^2
[z_grid, pi_z] = discretizeAR1_Tauchen(0,Params.rho,sqrt((1-Params.rho^2)*Params.sigma^2),n_z,Tauchen_q);
% Note: Nowadays Tauchen is outdated and you would be better of using discretizeAR1_FarmerToda(), 
% but following what Aiyagari (1994) originally did and using Tauchen here.

% Moments of the discretized process (in logs)
[z_mean, z_variance, z_corr, ~] = MarkovChainMoments(z_grid, pi_z);

% Convert from logs to levels
z_grid = exp(z_grid);

% Normalize labor productivity so that E[z] = 1 (taking exp(z_grid) makes this slightly off)
[Expectation_l, ~, ~, ~] = MarkovChainMoments(z_grid, pi_z);
z_grid = z_grid/Expectation_l;
[Expectation_l, ~, ~, ~] = MarkovChainMoments(z_grid, pi_z);

Params.Expectation_l = Expectation_l;

%% Asset grid

% Deterministic steady state (no idiosyncratic risk)
r_ss = 1/Params.beta - 1;
K_ss = ((r_ss + Params.delta) / Params.alpha)^(1/(Params.alpha - 1));

% Asset grid (more points near zero)
k_min  = 0;
k_max  = 10 * K_ss;
k_grid = k_min + (k_max - k_min) * (linspace(0, 1, n_k).^3)'; % ^3 is putting more points near 0, which is where value fn has higher curvature, this is more accurate

% Toolkit notation
n_d = 0; % No decision choice variable (i.e. no labor supply)
n_a = n_k;

d_grid = []; % will be ignored, as n_d=0
a_grid = k_grid;

%% Discount factor and Return function

DiscountFactorParamNames = {'beta'};

% Return function
% First inputs must be: a', a, z
% Everything after this is interpreted as a parameter
ReturnFn = @(aprime, a, z, alpha, delta, mu, r) ...
    Aiyagari1994_ReturnFn(aprime, a, z, alpha, delta, mu, r);
% You should take a look at the contents of Aiyagari1994_ReturnFn(), this
% is where most of the household problem is written out.

%% Functions to evaluate on the distribution

% Aggregate capital
FnsToEvaluate.K = @(aprime, a, z) a;
% Note that we called this K, and so we can use this name in the GeneralEqmEqns below

%% General Equilibrium
% In the Aiyagari (1994) model, general equilibrium is about finding the
% interest rate r to satisfy the capital market clearing condition.

GEPriceParamNames = {'r'};

% Note: with idiosyncratic risk, r is bounded above by the deterministic steady state
Params.r = 0.038; % Initial guess

%% General equilibrium condition

% Capital market clearing:
% r = alpha * K^{alpha-1} * L^{1-alpha} - delta
GeneralEqmEqns.CapitalMarket = @(r, K, alpha, delta, Expectation_l) ...
    r - (alpha * K^(alpha - 1) * Expectation_l^(1 - alpha) - delta);
% GeneralEqmEqns should be written to evaluate to zero in general equilibrium.

fprintf('Grid sizes: %d asset points, %d productivity states \n', n_a, n_z);

%% Solve for the stationary general equilibrium

vfoptions  = struct();  % Default options for value function iteration
simoptions = struct();  % Default options for stationary distribution

% If you run out of memory, you can enable low-memory mode (but will be slower):
% vfoptions.lowmemory = 1;

% If you want to solve the VFI using linear interpolation (GPU only), uncomment the following:
% vfoptions.gridinterplayer  = 1;
% vfoptions.ngridinterp      = 15;   % number of extra points between grid points
% simoptions.gridinterplayer = vfoptions.gridinterplayer;
% simoptions.ngridinterp     = vfoptions.ngridinterp;

heteroagentoptions.verbose = 1;  % Display progress information

fprintf('Calculating price vector corresponding to the stationary general equilibrium\n');

[p_eqm, GeneralEqmCondn] = HeteroAgentStationaryEqm_Case1(n_d, n_a, n_z, 0, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, [], [], [], GEPriceParamNames, heteroagentoptions, simoptions, vfoptions);

fprintf(' \n');
fprintf('The equilibrium values of the GE prices: %2.4f \n', p_eqm.r);

%% Equilibrium objects

% Now that we found the equilibrium interest rate, we have to put it into Params so we can use it
Params.r = p_eqm.r;

% Wage from firm first-order conditions
Params.w = (1-Params.alpha)*((Params.r + Params.delta) / Params.alpha)^(Params.alpha / (Params.alpha - 1));

fprintf('Calculating various equilibrium objects \n');

[V, Policy] = ValueFnIter_Case1(n_d, n_a, n_z, d_grid, a_grid, z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, [], vfoptions);

% Convert policy from indices to values
PolicyValues = PolicyInd2Val_Case1(Policy, n_d, n_a, n_z, d_grid, a_grid, vfoptions);

% Stationary distribution
StationaryDist = StationaryDist_Case1(Policy, n_d, n_a, n_z, pi_z, simoptions);

% Aggregate variables
AggVars = EvalFnOnAgentDist_AggVars_Case1(StationaryDist, Policy, FnsToEvaluate, Params, [], n_d, n_a, n_z, d_grid, a_grid, z_grid, simoptions);

fprintf('Aggregate capital: %f \n', AggVars.K.Mean);
% Note that AggVars uses the name 'K' we gave to the FnsThEvaluate

% Aggregate savings rate:
% Y = K^alpha (since L = 1), savings = delta*K
% s = delta*K / Y = delta*K^{1-alpha}
aggsavingsrate = Params.delta * AggVars.K.Mean^(1 - Params.alpha);

%% Distributional statistics

FnsToEvaluate2.Earnings    = @(aprime, a, z, w) w * z;
FnsToEvaluate2.Income      = @(aprime, a, z, r, w) w * z + (1 + r) * a;
FnsToEvaluate2.Wealth      = @(aprime, a, z) a;
FnsToEvaluate2.Consumption = @(aprime, a, z, alpha, delta, r) ...
    Aiyagari1994_ConsumptionFn(aprime, a, z, alpha, delta, r);

% Values on the grid
ValuesOnGrid = EvalFnOnAgentDist_ValuesOnGrid_Case1(Policy, FnsToEvaluate2, Params, [], n_d, n_a, n_z, d_grid, a_grid, z_grid, simoptions);

% Lorenz curves and inequality statistics (only on GPU)
if gpuDeviceCount>0
    simoptions.npoints = 1000;
    AllStats = EvalFnOnAgentDist_AllStats_Case1(StationaryDist, Policy, FnsToEvaluate2, Params,[],n_d, n_a, n_z, d_grid, a_grid, z_grid, simoptions);
    fprintf(' \n');
    fprintf('Distributions of Earnings and Wealth \n');
    fprintf('Gini (Earnings): %8.4f \n', AllStats.Earnings.Gini);
    fprintf('Gini (Income):   %8.4f \n', AllStats.Income.Gini);
    fprintf('Gini (Wealth):   %8.4f \n', AllStats.Wealth.Gini);
end

%% Output
fprintf('For parameter values sigma=%.2f, mu=%.2f, rho=%.2f \n', ...
    Params.sigma, Params.mu, Params.rho);
fprintf('Table 1 moments: sigma=%.4f, rho=%.4f \n', ...
    sqrt(z_variance), z_corr);

fprintf('Equilibrium interest rate r = %.4f%% \n', 100 * p_eqm.r);
fprintf('Aggregate savings rate s = %.4f \n', aggsavingsrate);

%% Plots

% Policy function for next-period assets
figure
plot(a_grid, a_grid, '--', 'LineWidth', 2)
hold on
plot(a_grid, PolicyValues(1,:,1),   'LineWidth', 2)
plot(a_grid, PolicyValues(1,:,n_z), 'LineWidth', 2)
legend('45 line', 'z_{1}', 'z_{nz}')
grid on
title('Policy function a''(a,z)')
ylabel('Next-period assets')
xlabel('Current assets')

% Policy function for consumption
figure
plot(a_grid, ValuesOnGrid.Consumption(:,1),   'LineWidth', 2)
hold on
plot(a_grid, ValuesOnGrid.Consumption(:,n_z), 'LineWidth', 2)
legend('z_{1}', 'z_{nz}')
grid on
title('Policy function c(a,z)')
ylabel('Consumption')
xlabel('Current assets')

% Asset cumulative distribution function
figure
plot(cumsum(sum(StationaryDist, 2)), 'LineWidth', 2)
grid on
title('Cumulative distribution of assets')
ylabel('Probability')
xlabel('Current assets')
