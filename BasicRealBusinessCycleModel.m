%%% Example using a variant of the Basic RBC model (following Aruoba, Fernandez-Villaverde, & Rubio-Ramirez, 2006)

tic;
UseAlternativeParams=1;

%% Set up
% Skip setting any options and just use the defaults. There is thus no need
% to pass these as inputs to the relevant commands (as was done in StochasticNeoClassicalGrowthModel example).
% tauchenoptions.parallel=2;
% vfoptions.parallel=2;

%Aruoba, Fernandez-Villaverde, & Rubio-Ramirez (2006) use 40 points for z
%and 25000 points for a (they use functional maximization so have no d
%grid). [Their runtime was about one week.]
n_d=50;
n_a=250;
n_z=21;  %Note, for the figures to correctly use z=0 this must be an odd number (make it 39)
q=3; %Parameter for the Tauchen method

%Discounting rate
Params.beta = 0.9896;

%Parameter values
Params.alpha = 0.4; % alpha
Params.theta=0.357;
Params.rho = 0.95; % rho
Params.tau=2;
Params.delta = 0.0196; % delta
Params.sigmasq_epsilon=0.007^2;

if UseAlternativeParams==1
    Params.tau=50;
    Params.sigmasq_epsilon=0.035^2;
end

% Params has been created as a structure. You can create the individual
% parameters from the structure by running the following command
CreateIndividualParams(Params)


%% Compute the steady state (just use these when picking grids)
varphi=(1/alpha*(1/beta-1+delta))^(1/(1-alpha));
Psi=theta/(1-theta)*(1-alpha)*varphi^(-alpha);
Omega=varphi^(1-alpha)-delta;
K_ss=Psi/(Omega+varphi*Psi);


%% Create grids (it is very important that each grid is defined as a column vector)
[z_grid, pi_z]=TauchenMethod(0,sigmasq_epsilon,rho,n_z,q); %[states, transmatrix]=TauchenMethod_Param(mew,sigmasq,rho,znum,q,[tauchenoptions]), transmatix is (z,zprime)

a_grid=linspace(0.01,5*K_ss,n_a)';
d_grid=linspace(0,1,n_d)';

%% Now, create the return function 
DiscountFactorParamNames={'beta'};

ReturnFn=@(d_val,aprime_val, a_val, s_val,alpha,delta,theta,tau) BasicRealBusinessCycleModel_ReturnFn(d_val,aprime_val, a_val, s_val,alpha,delta,theta,tau);
ReturnFnParams={'alpha','delta','theta','tau'}; %It is important that these are in same order as they appear in 'BasicRealBusinessCycleModel_ReturnFn'

%% Solve
%Do the value function iteration. Returns both the value function itself, and the optimal policy function.
V0=zeros(n_a,n_z);
[V,Policy]=ValueFnIter_Case1(V0, n_d,n_a,n_z,d_grid,a_grid,z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, ReturnFnParams);
time1=toc;

%% Report some output
fprintf('Run time for value function iteration: %8.2f seconds \n', time1)

%% Graph of the value function.
surf(a_grid*ones(1,n_z), ones(n_a,1)*z_grid', V)


