% Example using Epstein-Zin preferences
% Based on Caldara, Fernandez-Villaverde, Rubio-Ramirez, & Yao (2012)
% Stochastic Variance and Recursive (Epstein-Zin) Parameters
% Henceforth CFVRR2012
%
% There are a few things I do not replicate. Rf (in Tables 2/4) and
% Welfare Cost of Business Cycles in (Tables 3/5) I skip as paper does not
% seem to describe the exact calculation of these. I skip Euler eqn errors for
% discrete VFI these will be anyway be large (I say this based on past
% experience).
%
% Outputs: Figure 1 relates to the grid on z; Figure 2 shows policy/decision rules
% (compare to CFVRRY2012 Figures 1&2); Figure 3 shows densities (actually
% historgrams; compare to CFVRRY2012 Figures 3&4); Table is printed out
% (compare to CFVRRY2012 Tables 2&4).
%
% Discrete VFI is easily able to get everything but the investment
% volatility. It is likely that if the policy fn where interpolated during
% the time series simulations then you would get this too.
%
% Only place I 'cheat' relative to CFVRRY2012 is using a grid for capital which concentrates
% points near the steady state. However this is anyway a standard
% practice when creating grids for capital (CFVRRY2012 use a purely uniform
% grid).

ParamCombination=2; %A number from 1 to 8; determines which values of psi & gamma to use.
% ParamCombination=2 gives 'baseline' case of CFVRRY2012 (gamma=5, psi=0.5)
% ParamCombination=4 gives 'extreme' case of CFVRRY2012 (gamma=40, psi=0.5)
ExtremeCase=0; %To use the 'extreme' case of CFVRRY2012 you also must set this to 1 (default=0) to change values of sigmabar and eta

vfoptions.exoticpreferences=2; %Use Epstein-Zin preferences 
% If you set this zero you would just solve the same model but with
% standard recursive vonNeumann-Morgerstern expected utility preferences
% (parameters gamma and psi would be irrelevant to model). This is intended
% as purely illustrative. A serious comparison of the preference types would 
% require you to recalibrate the model.

%% Grid sizes
n_k=200;
n_l=31;
n_z=25;     %25 points for productivity (z)
n_sigma=5;  %5 points for volatility (sigma)
sprintf('Grid sizes n_k=%i, n_l=%i, n_z=%i, n_sigma=%i',n_k,n_l,n_z,n_sigma)
vfoptions.lowmemory=0; % If you want to make the grids larger you will need to use vfoptions.lowmemory=1

%% Declare parameters
ParamCombinationSub=ind2sub_homemade([2,4],ParamCombination);

%discount factor
Params.beta=0.991;
%utility
psiMatrix=[0.5,1.5];
Params.psi=psiMatrix(ParamCombinationSub(1));
Params.upsilon=0.357; %Paper contains a typo. Table 1 is correct, but then first paragraph of "Calibration" section accidently write theta=0.357 when it should say upsilon=0.357.
gammaMatrix=[2,5,10,40];
Params.gamma=gammaMatrix(ParamCombinationSub(2));
Params.theta=(1-Params.gamma)/(1-1/Params.psi); %Eqn at top of pg 3 (190)
%budget constraint
Params.zeta=0.3;
Params.delta=0.0196;
%z
Params.lambda=0.95;
%sigma AR(1)
Params.eta=0.06;
Params.rho=0.9;
Params.sigmabar=log(0.007); %Paper contains a typo. Table 1 reports value of exp(sigmabar) (for sigmabar defined in eqns on pg 190), not the value of log(sigmabar) which is what the table describes it as.
                       %Making this correction means that the std dev of
                       %the shocks to z [the exp(sigma)*omega term] is
                       %0.007, which would be in line with other papers.
if ExtremeCase==1
    Params.sigmabar=log(0.021);
    Params.eta=0.1;
end
 

%% Create the return function
if vfoptions.exoticpreferences==2 % This if statement is only needed so that you can easily check what happens if not using Epstein-Zin preferences
    DiscountFactorParamNames={'beta','gamma','psi'}; % The 'Epstein-Zin parameters'
elseif vfoptions.exoticpreferences==0
    DiscountFactorParamNames={'beta'};
end

ReturnFn =@(l,kprime,k,z,zeta,delta,upsilon) EpsteinZinPreferences_FmatrixFn(l,kprime,k,z,zeta, delta, upsilon);
ReturnFnParamNames={'zeta','delta','upsilon'};


%% Create grids
% Create grids using CFVRRY2012 notation of l & k, then convert into VFI toolkit notation of d & a

Omega=(1/Params.zeta*(1/Params.beta-1+Params.delta))^(1/(Params.zeta-1)); % pg 205 of CFVRRY2012
Phi=(Params.upsilon/(1-Params.upsilon))*(1-Params.zeta)*(Omega^Params.zeta); % pg 205 of CFVRRY2012
k_ss=Phi*Omega/(Omega^Params.zeta-Params.delta*Omega+Phi); % pg 205 of CFVRRY2012

k_grid=sort([linspace(3/4*k_ss,1.5*k_ss,ceil(n_k/3)),linspace(k_ss/2,2*k_ss,floor(n_k/3)),linspace(0.01,35,floor(n_k/3))])'; % Use 0 to 35 so at graphs of policy functions can have correct x-axis
% CFVRRY2012: mid pg 9 says a uniform grid (I use two overlaid uniform grids as we know most action is around the steady-state), top of pg 10 says 3000 points.

% Figure out the index of k_ss in k_grid (needed for drawing some graphs later on; not needed to solve model)
[~,k_ss_index]=min(abs(k_grid-k_ss));


%My methods also require a grid for labour supply l
l_grid=linspace(0.1,0.65,n_l)'; % In theory you might prefer to use 0 to 1, but a quick investigation shows anything outside this range would be wasted  (investigation shows=solve with 0 to 1, look at actual policy function)
if vfoptions.lowmemory==1
    l_grid=linspace(0,1,n_l)';
end

d_grid=l_grid;
a_grid=k_grid;
n_d=length(d_grid); %kprime by l
n_a=length(a_grid); % k


%% These lines (if uncommented) create the z_grid and sigma_grid following CFVRRY2012
% This implements my understanding of Section 3.3 (mid pg 196 of CFVRR2012)
[sigma_grid, pi_sigma]=TauchenMethod((1-Params.rho)*Params.sigmabar,Params.eta^2,Params.rho,n_sigma,1);
sigma_grid=exp(sigma_grid);

[z_grid,pi_z]=TauchenMethod(0,sigma_grid(1)^2,Params.lambda,n_z,1);
s_grid=z_grid;
pi_s=kron(pi_sigma(1,:),pi_z);
for sigma_c=2:n_sigma
    [z_grid,pi_z]=TauchenMethod(0,sigma_grid(sigma_c)^2,Params.lambda,n_z,1);
    % The following two should be predeclared, but am feeling lazy
    s_grid=[s_grid; z_grid];
    pi_s=[pi_s; kron(pi_sigma(sigma_c,:),pi_z)];
end

n_s=length(s_grid); %z by sigma

% Just out of interest
figure(1)
plot(sort(s_grid))

%% Now, do the value function iteration
vfoptions.howards=20;
vfoptions.verbose=1;
V0=zeros(n_a,n_s,'gpuArray');
[V, Policy]=ValueFnIter_Case1(V0, n_d,n_a,n_s,d_grid,a_grid,s_grid, pi_s, DiscountFactorParamNames, ReturnFn, vfoptions,Params,ReturnFnParamNames);

%% Generate output based on the solution

midpointofs=ceil(n_s/2);
midpointsofz=ceil(n_z/2)+0:n_z:n_s-1;

%% 1. Plot the Decision Rules (Figure 1/2)
PolicyMatrix_l=l_grid(shiftdim(Policy(1,:,:),1));
PolicyMatrix_kprime=k_grid(shiftdim(Policy(2,:,:),1));
% Decision rule for consumption
c_fn=@(l,kprime,k,z,zeta,delta) exp(z)*(k^zeta)*(l^(1-zeta))+(1-delta)*k-kprime;
PolicyMatrix_c=arrayfun(c_fn, PolicyMatrix_l, PolicyMatrix_kprime, k_grid, s_grid',Params.zeta,Params.delta); % Value

figure(2)
subplot(2,3,1), plot(k_grid,PolicyMatrix_c(:,midpointofs))
subplot(2,3,2), plot(k_grid,PolicyMatrix_l(:,midpointofs))
subplot(2,3,3), plot(k_grid,PolicyMatrix_kprime(:,midpointofs))
subplot(2,3,4), plot(sigma_grid,PolicyMatrix_c(k_ss_index,midpointsofz))
subplot(2,3,5), plot(sigma_grid,PolicyMatrix_l(k_ss_index,midpointsofz))
subplot(2,3,6), plot(sigma_grid,PolicyMatrix_kprime(k_ss_index,midpointsofz))


%% Set simulation options for toolkit and do simulations that are used for 2 & 3
% simoptions are set following Section 5.2 of CFVRRY2012 (bottom pg. 198).
% They start simulations from the deterministic steady-state, I just use
% VFI toolkit default (mid-point of grid; you could overrule this using simoptions.seedpoint).
simoptions.burnin=10^3;
simoptions.simperiods=10^4;
simoptions.iterate=0;

TimeSeriesFn_c=@(l,kprime,k,z) exp(z)*(k^Params.zeta)*(l^(1-Params.zeta))+(1-Params.delta)*k-kprime;
TimeSeriesFn_l=@(l,kprime,k,z) l;
TimeSeriesFn_k=@(l,kprime,k,z) k;
TimeSeriesFn_y=@(l,kprime,k,z) exp(z)*(k^Params.zeta)*(l^(1-Params.zeta));
TimeSeriesFn_i=@(l,kprime,k,z) kprime-(1-Params.delta)*k;
TimeSeriesFn_Rk=@(l,kprime,k,z) Params.zeta*exp(z)*(k^(Params.zeta-1))*(l^(1-Params.zeta))-Params.delta;
TimeSeriesFn={TimeSeriesFn_c,TimeSeriesFn_l,TimeSeriesFn_k,TimeSeriesFn_y,TimeSeriesFn_i,TimeSeriesFn_Rk};

TimeSeriesIndexes=SimTimeSeriesIndexes_Case1(Policy,n_d,n_a,n_s,pi_s,simoptions);
TimeSeries=TimeSeries_Case1(TimeSeriesIndexes, Policy, TimeSeriesFn, n_d, n_a, n_s, d_grid, a_grid, s_grid,simoptions);

%% 2. Plot the Densities (Figure 3/4)
% Not clear in CFVRRY2012 exactly how the density plots are created. I use
% the matlab histogram function 'histcounts'.

nbins=30;

[densityfn_c, edges, bins]=histcounts(TimeSeries(1,:),nbins);
xaxis_c=edges(1:end-1)+(edges(2:end)-edges(1:end-1))/2;
[densityfn_l, edges, bins]=histcounts(TimeSeries(2,:),nbins);
xaxis_l=edges(1:end-1)+(edges(2:end)-edges(1:end-1))/2;
[densityfn_k, edges, bins]=histcounts(TimeSeries(3,:),nbins);
xaxis_k=edges(1:end-1)+(edges(2:end)-edges(1:end-1))/2;
[densityfn_y, edges, bins]=histcounts(TimeSeries(4,:),nbins);
xaxis_y=edges(1:end-1)+(edges(2:end)-edges(1:end-1))/2;
% [densityfn_Rf, edges, bins]=histcounts(TimeSeries(2,:),100);
% xaxis_Rf=edges(1:end-1)+(edges(2:end)-edges(1:end-1))/2;
[densityfn_Rk, edges, bins]=histcounts(TimeSeries(6,:),nbins);
xaxis_Rk=edges(1:end-1)+(edges(2:end)-edges(1:end-1))/2;

% [densityfn_c,xaxis_c] = ksdensity(gather(TimeSeries(1,:)));
% [densityfn_l,xaxis_l] = ksdensity(gather(TimeSeries(2,:)));
% [densityfn_k,xaxis_k] = ksdensity(gather(TimeSeries(3,:)));
% [densityfn_y,xaxis_y] = ksdensity(gather(TimeSeries(4,:)));
% % [densityfn_Rf,xaxis_Rf] = ksdensity(gather(TimeSeries(1,:)));
% [densityfn_Rk,xaxis_Rk] = ksdensity(gather(TimeSeries(6,:)));

figure(3)
subplot(3,2,1), plot(xaxis_c,densityfn_c)
subplot(3,2,2), plot(xaxis_l,densityfn_l)
subplot(3,2,3), plot(xaxis_k,densityfn_k)
subplot(3,2,4), plot(xaxis_y,densityfn_y)
% subplot(3,2,5), plot(xaxis_Rf,densityfn_Rf)
subplot(3,2,6), plot(xaxis_Rk,densityfn_Rk)


%% 3. Numbers Relating to Business Cycles (Table 2/5) 

BusCycStats_means=mean(TimeSeries,2);
BusCycStats_var=100*std([log(TimeSeries(1:5,:));TimeSeries(6,:)],[],2).^2;

Table=[BusCycStats_means(1), BusCycStats_means(4), BusCycStats_means(5), nan, 100*BusCycStats_means(6); ...
    BusCycStats_var(1), BusCycStats_var(4), BusCycStats_var(5),nan, BusCycStats_var(6)]

%% Save all of the parts that relate to output
save ./SavedOutput/CaldaraFVillaverdeRRamirezYao2012.mat k_grid sigma_grid V Policy PolicyMatrix_c PolicyMatrix_l PolicyMatrix_kprime midpointofs midpointsofz k_ss_index TimeSeries Table















