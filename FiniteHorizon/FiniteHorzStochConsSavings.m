% This code implements a simple ten-period stochastic consumption savings
% problem to demonstrate the finite value function command of the VFI
% Toolkit.

% Main strength of the finite value function command is that you simply
% declare parameters that depend on age as row vectors, and the VFI Toolkit
% automatically realizes that these are age dependent.


% Agents live ten periods. Income consists of a deterministic function of
% age Wj plus an AR(1) shock Wz.

% Age
% One endogenous variable: assets
% One stochastic exogenous variables: income shock

J=10; % Ages 1 to 10 inclusive.

% Grid sizes to use
n_a=501;
n_z=15; % income
N_j=J; % Number of periods in finite horizon


%% Declare the model parameters

Params.gamma=3;
% Gamma plays three roles:      
        % gamma is coefficient of relative risk aversion
        % 1/gamma is the intertemporal elasticity of substitution of consumption
        % gamma+1 is the coefficient of relative prudence (Kimball, 1990)

Params.beta=0.96; % rate at which agents discount future
Params.r=0.03; % interest rate on assets

% Declare the age dependent parameters. This is a simple matter of creating
% the parameter as a row vector of length J (the VFI Toolkit automatically
% detects which parameters depend on age and which do not).
Params.Wj=[1,2,3,5,7,8,8,5,4,4]; % deterministic income depends on age

% Stochastic Wz: use Tauchen method to discretize the AR(1) process log(Wz):
Params.Wz_rho=0.7;
Params.Wz_sigmasqepsilon=0.05;
Params.Wz_sigmasqu=Params.Wz_sigmasqepsilon./(1-Params.Wz_rho.^2);
Params.q=3; % For tauchen method
[z_grid, pi_z]=TauchenMethod(0,Params.Wz_sigmasqu, Params.Wz_rho, n_z, Params.q);

%% Grids
maxa=150;
a_grid=linspace(0,maxa,n_a)'; % Could probably do better by adding more grid points near zero

%% Now, create the return function 
DiscountFactorParamNames={'beta'};

ReturnFn=@(aprime,a,Wz,gamma,r,Wj) FiniteHorzStochConsSavings_ReturnFn(aprime,a,Wz,gamma,r,Wj)
ReturnFnParamNames={'gamma','r','Wj'}; %It is important that these are in same order as they appear in 'FiniteHorzStochConsSavings_ReturnFn'

%% Now solve the value function iteration problem

vfoptions.verbose=0;
tic;
[V, Policy]=ValueFnIter_Case1_FHorz(0,n_a,n_z,N_j, 0, a_grid, z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, ReturnFnParamNames,vfoptions);
toc

% max(max(max(max(Policy))))<n_a % Double check that never try to leave top of asset grid.

% V is assets-by-Wz-by-age


%% Plot of how Value function in asset holdings and age 
% (three subplots for minimum, maximum, and average values of stochastic component of income)

subplot(2,2,1)
surf(a_grid*ones(1,N_j),ones(n_a,1)*(1:1:N_j),reshape(V(:,1,:),[n_a,N_j]))
subplot(2,2,2)
surf(a_grid*ones(1,N_j),ones(n_a,1)*(1:1:N_j),reshape(V(:,n_z,:),[n_a,N_j]))
subplot(2,2,3)
surf(a_grid*ones(1,N_j),ones(n_a,1)*(1:1:N_j),reshape(V(:,ceil(n_z/2),:),[n_a,N_j]))






