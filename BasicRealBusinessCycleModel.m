% Example using a variant of the Basic RBC model (following Aruoba, Fernandez-Villaverde, & Rubio-Ramirez, 2006)
%
% Note: On GPU this will take fractions of a second, on CPUs you will want to make the grids smaller (n_d,n_a,n_z) or it will take a while.

tic;
UseAlternativeParams=0;

%% Set up

% Number of grid points
n_d=50;
n_a=250;
n_z=21;  % Note, for the figures to correctly use z=0 this must be an odd number (make it 39)

% Discount factor
Params.beta = 0.9896;

% Parameter values
Params.alpha = 0.4; % alpha
Params.theta=0.357;
Params.rho = 0.95; % rho
Params.tau=2;
Params.delta = 0.0196; % delta
Params.sigma_epsilon=0.007;

if UseAlternativeParams==1
    Params.tau=50;
    Params.sigma_epsilon=0.035;
end

%% Compute the steady state (just use these when picking grids)
varphi=(1/Params.alpha*(1/Params.beta-1+Params.delta))^(1/(1-Params.alpha));
Psi=Params.theta/(1-Params.theta)*(1-Params.alpha)*varphi^(-Params.alpha);
Omega=varphi^(1-Params.alpha)-Params.delta;
K_ss=Psi/(Omega+varphi*Psi);


%% Create grids (it is very important that each grid is defined as a column vector)
Tauchen_q=3; % Parameter for the Tauchen method
[z_grid,pi_z]=discretizeAR1_Tauchen(0,Params.rho,Params.sigma_epsilon,n_z,Tauchen_q);

a_grid=linspace(0.01,5*K_ss,n_a)';
d_grid=linspace(0,1,n_d)';

%% Now, create the return function 
DiscountFactorParamNames={'beta'};

ReturnFn=@(d_val,aprime_val, a_val, s_val,alpha,delta,theta,tau) BasicRealBusinessCycleModel_ReturnFn(d_val,aprime_val, a_val, s_val,alpha,delta,theta,tau);

%% Solve
vfoptions.gridinterplayer=1; % for next period capital, interpolate between grid points...
vfoptions.ngridinterp=20; % ...with 20 more evenly spaced grid points

% Do the value function iteration. Returns both the value function itself, and the optimal policy function.
[V,Policy]=ValueFnIter_Case1(n_d,n_a,n_z,d_grid,a_grid,z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames,[],vfoptions);
time1=toc;

%% Report some output
fprintf('Run time for value function iteration: %8.2f seconds \n', time1)

%% Graph of the value function.
figure(1)
surf(a_grid*ones(1,n_z), ones(n_a,1)*z_grid', V)
title('Value function')

%% Policy contains indexes, and because we are using vfoptions.gridinterplayer=1 is not trivial to interpret
% So instead we can just convert it to PolicyVals which contains the values (by evaluating Policy based on the grids)
PolicyVals=PolicyInd2Val_Case1(Policy,n_d,n_a,n_z,d_grid,a_grid,vfoptions);

%% Graph of the policy function

figure(2)
plot(a_grid,a_grid,'--')
hold on
plot(a_grid,PolicyVals(2,:,1))
hold on
plot(a_grid,PolicyVals(2,:,n_z))
legend('45 line','Given z(1)','Given z(n_z)')
title('Policy function for next period capital')
xlabel('Current period capital')

figure(3)
plot(a_grid,PolicyVals(1,:,1))
hold on
plot(a_grid,PolicyVals(1,:,n_z))
legend('Given z(1)','Given z(n_z)')
title('Policy function for labor')
xlabel('Current period capital')

