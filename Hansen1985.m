%%% Example using the model of Hansen (1985)

% l_z=1; %tech shock
% l_a=1; %assets
% l_d=1; %decide whether to work

n_z=21;
n_a=101;
n_d=21;

%% Some Toolkit options
Parallel=2; % Use GPU

tauchenoptions.parallel=Parallel;
mcmomentsoptions.parallel=tauchenoptions.parallel;

vfoptions.parallel=Parallel;

simoptions.parallel=Parallel;

%% Setup

%Discounting rate
beta = 0.96;

%Parameter values
alpha = 0.33; % alpha
gamma=-2*log(1-0.53)/0.53;
rho = 0.95; % rho
delta = 0.10; % delta
sigmasq_epsilon=0.09;
%Step 1: Compute the steady state
K_ss=((alpha*beta)/(1-beta*(1-delta)))^(1/(1-alpha));

%Create grids (it is very important that each grid is defined as a column
%vector)
q=3;
[z_grid, pi_z]=TauchenMethod(0,sigmasq_epsilon,rho,n_z,q,tauchenoptions);
a_grid=linspace(0,2*K_ss,n_a)';
d_grid=linspace(0,1,n_d)';

n_z=length(z_grid);
n_a=length(a_grid);
n_d=length(d_grid);

ReturnFn=@(d_val, aprime_val, a_val, z_val, alpha, delta, gamma) Hansen1985_ReturnFn(d_val, aprime_val, a_val, z_val, alpha, delta, gamma)
ReturnFnParams=[alpha, delta, gamma];

%% Solve
V0=ones([n_a,n_z],'gpuArray'); %(a,z)
[V,Policy]=ValueFnIter_Case1(V0, n_d,n_a,n_z,d_grid,a_grid,z_grid, pi_z, beta, ReturnFn,vfoptions,ReturnFnParams);


SteadyStateDist=SteadyState_Case1_Simulation(Policy,n_d,n_a,n_z,pi_z, simoptions);
SteadyStateDist=SteadyState_Case1(SteadyStateDist,Policy,n_d,n_a,n_z,pi_z,simoptions);


%% Generate output
% NOT YET DONE












