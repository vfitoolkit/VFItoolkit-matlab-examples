% Imrohoroglu (1989) - Cost of Business Cycles with Indivisibilities and Liquidity Constraints

%Choose one of the following
%EconomyEnvironment pairs
%1: Economy=1, Environment=A
%1: Economy=2, Environment=A
%1: Economy=1, Environment=B
%1: Economy=2, Environment=B
%1: Economy=1, Environment=C
%1: Economy=2, Environment=C
EconomyEnvironment=1

%% Some Toolkit options
Parallel=2; % Use GPU

tauchenoptions.parallel=Parallel;

mcmomentsoptions.parallel=tauchenoptions.parallel;

vfoptions.parallel=Parallel;

simoptions.parallel=Parallel;

%% Setup

% num_z_vars=1; %This is 'n' in Imrohoroglu
% num_s_vars=1; %This is 'i' in Imrohoroglu
% num_a_vars=1; %This is 'a' in Imrohoroglu
% num_d_vars=0; %There is none in Imrohorglu (you only choose aprime, which Case 1 codes assume is chosen anyway)

%Discounting rate
beta = 0.995;

%Parameters
sigma=6.2; %Lucas (1987) uses 6.2; also compares for value of 1.5
theta=0.25;
r_b=8; %Interest rate on borrowing
r_l=0; %Interest rate on lending (in environment A, this is the only r)
y=1; %Normalization (THIS IS NOT FROM PAPER, I JUST GUESSED)

E1_z_grid=[1,2]';
E2_z_grid=[1]';
s_grid=[y,theta*y]';
E1_pi_sz=[0.9141,0.0234,0.0587, 0.0038; 0.5625, 0.3750, 0.0269, 0.0356; 0.0608, 0.0016, 0.8813, 0.0563; 0.0375, 0.0250, 0.4031, 0.5344];
E2_pi_s=[0.9565,0.0435;0.5,0.5];
E2_pi_z=[1]';
EA_a_grid=linspace(0,8,301)';
EB_a_grid=linspace(-8,8,601)';


if EconomyEnvironment==1 || EconomyEnvironment==3 || EconomyEnvironment==5
    z_grid=E1_z_grid;
    pi_sz=E1_pi_sz;
elseif EconomyEnvironment==2 || EconomyEnvironment==4 || EconomyEnvironment==6
    z_grid=E2_z_grid;
    pi_sz=kron(E2_pi_s,E2_pi_z);
end

if EconomyEnvironment==1 || EconomyEnvironment==4
    a_grid=EA_a_grid;
elseif EconomyEnvironment==2 || EconomyEnvironment==5
    a_grid=EB_a_grid;
elseif EconomyEnvironment==3 || EconomyEnvironment==6
    disp('EconomyEnvironment 3 & 6 not yet implemented')
end

%Note that from the point of view of the value function, there is no
%difference between s & z, so we combine them.
sz_grid=[s_grid;z_grid];

n_s=length(s_grid);
n_z=length(z_grid);
n_sz=[n_s,n_z];
n_a=length(a_grid);

n_d=0; % No d variable
d_grid=0;

ReturnFn=@(aprime_val, a_val, s_val,sigma,theta,y,r_l,r_b) Imrohoroglu1989_ReturnFn(aprime_val, a_val, s_val,sigma,theta,y,r_l,r_b);
ReturnFnParams=[sigma,theta,y,r_l,r_b];

%% Solve the model
%Do the value function iteration. Returns both the value function itself,
%and the optimal policy function.
V0=ones(n_a,n_s,n_z,'gpuArray');
[V,Policy]=ValueFnIter_Case1(V0, n_d, n_a, n_sz, d_grid, a_grid, sz_grid, pi_sz, beta, ReturnFn, vfoptions, ReturnFnParams);

SteadyStateDist=ones([n_a,n_s,n_z],'gpuArray')./prod([n_a,n_sz]);
SteadyStateDist=SteadyState_Case1(SteadyStateDist,Policy,n_d,n_a,n_sz,pi_sz, simoptions);

%% Generate some output following what is reported in Imrohoroglu (1989)
% NOT YET DONE





