% This code implements a simple ten-period stochastic consumption savings
% problem to demonstrate the finite value function command of the VFI
% Toolkit.
% This 'extended' version then goes on to simulate panel data and
% life-cycle profiles for earnings and assets. ('extension' begins on line 84)

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
[z_grid,pi_z]=discretizeAR1_FarmerToda(0,Params.Wz_rho,Params.Wz_sigmasqepsilon,n_z); % Farmer-Toda discretizes an AR(1) (is better than Tauchen, & better than Rouwenhorst for rho<0.99)

%% Grids
maxa=150;
a_grid=linspace(0,maxa,n_a)'; % Could probably do better by adding more grid points near zero

%% Now, create the return function 
DiscountFactorParamNames={'beta'};

ReturnFn=@(aprime,a,Wz,gamma,r,Wj) FiniteHorzStochConsSavings_ReturnFn(aprime,a,Wz,gamma,r,Wj)

%% Now solve the value function iteration problem

vfoptions.verbose=0;
tic;
[V, Policy]=ValueFnIter_Case1_FHorz(0,n_a,n_z,N_j, 0, a_grid, z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, [],vfoptions);
toc

% max(max(max(max(Policy))))<n_a % Double check that never try to leave top of asset grid.

% V is assets-by-Wz-by-age


%% Plot of how Value function in asset holdings and age 
% (three subplots for minimum, maximum, and average values of stochastic component of income)

figure(1)
subplot(2,2,1)
surf(a_grid*ones(1,N_j),ones(n_a,1)*(1:1:N_j),reshape(V(:,1,:),[n_a,N_j]))
subplot(2,2,2)
surf(a_grid*ones(1,N_j),ones(n_a,1)*(1:1:N_j),reshape(V(:,n_z,:),[n_a,N_j]))
subplot(2,2,3)
surf(a_grid*ones(1,N_j),ones(n_a,1)*(1:1:N_j),reshape(V(:,ceil(n_z/2),:),[n_a,N_j]))


%% Create a simulated panel data set from this model.
% Panel will contain 1000 life-times (J=10 period) simulations for earnings, assets and age.
simoptions.numbersims=1000;

% Define the ages (lets consider each period j to be 5 years, and this
% model to cover ages 21-70).
Params.age=21:5:70; % Have set age to be starting value of each five-year 'age bin'.

% 1. Define the 'functions to evaluate' for earnings and assets.
FnsToEvaluate.Assets = @(aprime,a,z) a; 
FnsToEvaluate.Earnings= @(aprime,a,z,Wj) Wj+exp(z);
FnsToEvaluate.age = @(aprime,a,z,age) age;
% As with return function the first inputs must be (any decision variables), next period endogenous
% state, this period endogenous state (any exogenous shocks). After that come any parameters.

% 2. Define the initial (probability) distribution from which households in simulation are 'drawn/born'. Also called the initial conditions.
InitialDist=zeros([n_a,n_z]);
InitialDist(1,ceil(n_z/2))=1; % All agents born with zero assets (a_grid(1)) and the 'average' shock (middle value of z_grid).

% 3. Simulate Panel Data
% Same variables as we used for the life-cycle profiles.
SimPanelValues=SimPanelValues_FHorz_Case1(InitialDist,Policy, FnsToEvaluate,[],Params,0,n_a,n_z,N_j,0,a_grid,z_grid,pi_z, simoptions);

% So for example one simulated life-time looks like
SimPanelValues.Assets(:,1)
SimPanelValues.Earnings(:,1)
SimPanelValues.age(:,1)
% The columns index the model period j (which in this model represents age).

%% Plot the mean life-cycle profiles for earnings and assets.

% We will create them from the 

% 1. Define the 'functions to evalute' for earnings and assets.
% No need to do this again, as we already have done so for the panel data.

simoptions.lifecyclepercentiles=4; % Just mean and median, no percentiles. (By default is 20, so also gives min and ventiles, the later includes max by definition as the 20th ventile.)
LifeCycleProfiles=SimLifeCycleProfiles_FHorz_Case1(InitialDist,Policy, FnsToEvaluate,Params,[],0,n_a,n_z,N_j,0,a_grid,z_grid,pi_z,simoptions);
% There is also 'LifeCycleProfiles_FHorz_Case1()' which uses the StationaryDist, rather than as a simulation from InitialDist

% Figure: Assets
figure(2)
plot(Params.age,LifeCycleProfiles.Assets.Mean,Params.age,LifeCycleProfiles.Assets.Median)
title('Life-cycle Profile of Assets')
legend('Mean', 'Median')
xlabel('Age')
% Figure: Earnings
figure(3)
plot(Params.age,LifeCycleProfiles.Earnings.Mean,Params.age,LifeCycleProfiles.Earnings.Median)
title('Life-cycle Profile of Earnings')
legend('Mean', 'Median')
xlabel('Age')

% Notice that since agents are not allowed to borrow (as minimum value of
% a_grid is zero) they simply consume all their income for first few
% periods and keep zero assets, then during the high incomes in mid-late
% life they save some of their income to consume during retirement when
% their income is lower and they run down the assets they saved.


