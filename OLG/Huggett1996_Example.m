% Replication of Huggett (1996) - Wealth Distribution in Life Cycle Economies

% I follow notation of Huggett (1996) for all parameters except the number
% of period (and retirement age) which I denote by J (JR) rather than N.

% One endogenous variable: assets
% One stochastic exogenous variables: income shock
% Age

Params.J=79; % Ages 20 to 98 inclusive.

% Grid sizes to use
n_a=2501;
% n_z=18; % income (these 18 points are hardcoded into z_grid and pi_z, done this way due to how Huggett sets them up)
N_j=Params.J; % Number of periods in finite horizon

% Huggett solves more than one economy (parametrization). This example just does his baseline.

%% Declare the model parameters
% Note that r, w, and G will be determined in General equilbrium, so these are really just initial guesses.

% Preference
Params.beta=1.011; % rate at which agents discount future (is greater than one due to the 'stochastic probability of death')
Params.sigma=1.5; % Risk-aversion

Params.bc_equalsminusw=0; % An indicator variable for the first of the two possible values for the borrowing constraint, 0 and -w

% Demographics
Params.JR=46; % Retirement at age 65 (Note: never actually used directly as is implicit in the deterministic earnings profile and the retirement benefits
Params.n=0.012; % Population growth rate of 1%

% Tax rates
% Params.tau % Determined based on r: Params.tau=0.195/(1-Params.delta*K/Y) % Note that K/Y can be calculated from r, see below
Params.tau=0.195/(1-0.06*3); % Just give an initial guess for tau here
Params.theta=0.1;
% Accidental bequests
% Params.T % Determined in GE
% Retirement benefits: I set them equal to b*bvec (in Huggetts notation this is just b which is itself a function of retirement status, 
% I seperate it into a scalar and an indicator on whether or not you are retired).
Params.bvec=[zeros(1,Params.JR-1),ones(1,1+Params.J-Params.JR)]; % Set it up as an age dependent vector that is zero before retirement age
% Note: bvec is really just an indicator of retirement

% Production function
Params.A=0.895944;
Params.alpha=0.36; % Capital-share in Cobb-Douglas production function (I currently live in New Zealand where this is the actual capital-share of GDP ;)
Params.delta=0.06; % Depreciation rate

% Survival probabilities, based on Jordan, C., Life contingencies, 2nd ed. (Society of Actuaries).
% Huggett (1996) dates it as the 1975 print, the following link is to the 1991 imprint. But since it is still the 2nd edition (which first appeared 1967)
% Seems to be the right numbers. (the pg 342 of the pdf, 346 of book, appears to be the closest thing)
% https://vdocuments.mx/download/life-contingencies-chester-wallace-jordanpdf
Params.dj=[0.00159, 0.00169, 0.00174, 0.00172, 0.00165, 0.00156, 0.00149, 0.00145, 0.00145, 0.00149,...
    0.00156, 0.00163, 0.00171, 0.00181, 0.00193, 0.00207, 0.00225, 0.00246, 0.00270, 0.00299,...
    0.00332, 0.00368, 0.00409, 0.00454, 0.00504, 0.00558, 0.00617, 0.00686, 0.00766, 0.00865,...
    0.00955, 0.01058, 0.01162, 0.01264, 0.01368, 0.01475, 0.01593, 0.01730, 0.01891, 0.02074,...
    0.02271, 0.02476, 0.02690, 0.02912, 0.03143, 0.03389, 0.03652, 0.03930, 0.04225, 0.04538,...
    0.04871, 0.05230, 0.05623, 0.06060, 0.06542, 0.07066, 0.07636, 0.08271, 0.08986, 0.09788,...
    0.10732, 0.11799, 0.12895, 0.13920, 0.14861, 0.16039, 0.17303, 0.18665, 0.20194, 0.21877,...
    0.23601, 0.25289, 0.26973, 0.28612, 0.30128, 0.31416, 0.32915, 0.34450, 0.36018]; % conditional probability of death
Params.sj=1-Params.dj; % Conditional survival probabilities. Act as a discount rate.

% Declare the age dependent parameters. This is a simple matter of creating
% the parameter as a row vector of length J (the VFI Toolkit automatically
% detects which parameters depend on age and which do not).
% Params.ybarj=log(0.5289)+log([linspace(0.3,1.25,12),linspace(1.3,1.5,6),linspace(1.5,1.5,11),linspace(1.49,1.05,11),linspace(1,0.1,9),linspace(0.08,0,5),zeros(1,Params.J-54)]); % deterministic income depends on age
Corbae_deterministicincome = [0.0911 0.1573 0.2268 0.2752 0.3218 0.3669 0.4114 0.4559 0.4859 0.5164 0.5474 0.5786 0.6097 0.6311 0.6517 0.6711 0.6893 0.7060 0.7213 0.7355 0.7489 0.7619 0.7747 0.7783 0.7825 0.7874 0.7931 0.7994 0.7923 0.7850 0.7771 0.7679 0.7567 0.7351 0.7105 0.6822 0.6500 0.6138 0.5675 0.5183 0.4672 0.3935 0.3239 0.2596 0.1955 0.1408 0.0959 0.0604 0.0459 0.0342 0.0246 0.0165 0.0091 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]; % Figure 1 of Huggett (1996) shows: plot(19+(1:1:Params.J),exp(Params.ybarj))
Params.ybarj=log(Corbae_deterministicincome);
% The log(0.5289)+ come from a 2014 Homework Handout by Dean Corbae who
% states that he got them from Mark Huggett (I got the Corbae handout files from a friend).
% The only thing that is odd is the normalization by 0.5289. I can't find
% mention of this in Huggett (1996) but according to the Corbae Homework
% it is about getting the ratio of 'work' L to population N correct; to
% 0.5289 as in data.

% Stochastic income, y:
Params.gamma=0.96;
Params.sigmasqepsilon=0.045;
Params.sigmasqy=Params.sigmasqepsilon./(1-Params.gamma.^2); 
% Initial distribution of income
Params.sigmasqy1=0.38;

% Huggett use 17 states equally spaced between +-4*sigmasqy1, with an 18th
% state at 6*sigmasqy1, "The transition probabilities between states are
% calculated by integrating the area under the normal distribution conditional on
% the current value of the state."
n_z=18;
z_grid=[linspace(-4*sqrt(Params.sigmasqy1),4*sqrt(Params.sigmasqy1),17),6*sqrt(Params.sigmasqy1)]';
pi_z=nan(18,18);
% Following lines implement the transition matrix, they are largely just a copy of some code from the TauchenMethod() command.
sigma=sqrt(Params.sigmasqepsilon); %stddev of e
for ii=1:length(z_grid)
    pi_z(ii,1)=normcdf(z_grid(1)+(z_grid(2)-z_grid(1))/2-Params.gamma*z_grid(ii),0,sigma);
    for jj=2:(length(z_grid)-1)
        pi_z(ii,jj)=normcdf(z_grid(jj)+(z_grid(jj+1)-z_grid(jj))/2-Params.gamma*z_grid(ii),0,sigma)-normcdf(z_grid(jj)-(z_grid(jj)-z_grid(jj-1))/2-Params.gamma*z_grid(ii),0,sigma);
    end
    pi_z(ii,end)=1-normcdf(z_grid(end)-(z_grid(end)-z_grid(end-1))/2-Params.gamma*z_grid(ii),0,sigma);
end
if gpuDeviceCount>0 % This is just so code can be run with or without gpu
    z_grid=gpuArray(z_grid);
    pi_z=gpuArray(pi_z);
end
% Double check: cumsum(pi_z,2) shows all each row adding to one.

%% General eqm variables: give some initial values
GEPriceParamNames={'r','b','T'};
Params.r=0.06; % interest rate on assets
Params.b=1.2; % Benefits level for retirees
Params.T=0.2; % lumpsum transfers made out of the accidental bequests
% I originally had b=0.8, T=0.6; switched to these after solving for the GE as I know they are
% closer to the true values, which helps make things run a bit faster.

% Following are coded into the return function to get the values of w and tau given the value of r
% KdivL=((Params.r+Params.delta)/(Params.alpha*Params.A))^(1/(Params.alpha-1));
% KdivY=(KdivL^(1-Params.alpha))/Params.A;
% Params.w=Params.A*(1-Params.alpha)*(KdivL^Params.alpha); % wage rate (per effective labour unit)
% Params.tau=0.195/(1-Params.delta*KdivY);

% Note that G is not part of the GEPriceParamNames, this is because it is
% effectively just a residual of the model and plays no part in the actual
% computation. It doesn't effect anyones behaviour, so we don't need it to
% solve the model, and once we solve the model if we want to know what is it
% we can just calculate it based on the governments budget balance equation.


%% Grids
% w is dec in r, so on the assumption that r will not be any lower than zero, we know that w will not be higher than:
maxw=Params.A*(1-Params.alpha)*((((0+Params.delta)/(Params.alpha*Params.A))^(1/(Params.alpha-1)))^Params.alpha);
maxa=250; % Note: If anyone actually had this many assets in first 5 periods they would hold on to them all and this maxa would be too low, but noone will actually ever get near it in first few periods, and by the time they do they would not want to hold such high level of assets.
a_grid=[linspace(-maxw,15,ceil(n_a/2))'; linspace(15,maxa,n_a-ceil(n_a/2))']; % Chose 15 as this is anyway greater than the mean
a_grid(ceil(n_a/2))=0; %There are two 15 values, replace one with 0 
a_grid=sort(a_grid); % Double-check: length(unique(a_grid))==n_a
% Note that the borrowing constraint is instead enforced inside the return
% function. This means some grid points are wasted, but was a bit cleaner.

%% Now, create the return function 
DiscountFactorParamNames={'beta','sj'};
 
ReturnFn=@(aprime,a,z,sigma,r,ybarj,theta,b,bvec,T,delta,alpha,A,bc_equalsminusw) Huggett1996_ReturnFn(aprime,a,z,sigma,r,ybarj,theta,b,bvec,T,delta,alpha,A,bc_equalsminusw)
% For the return function the first inputs must be (any decision variables), next period endogenous
% state, this period endogenous state (any exogenous shocks). After that come any parameters.

%% Now solve the value function iteration problem, just to check that things are working before we go to General Equilbrium

disp('Test ValueFnIter')
% vfoptions.verbose=0;
% vfoptions.policy_forceintegertype=2; % Policy was not being treated as integers (one of the elements was 10^(-15) different from an integer)
vfoptions=struct(); % Just use the defaults
tic;
[V, Policy]=ValueFnIter_Case1_FHorz(0,n_a,n_z,N_j, 0, a_grid, z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, [],vfoptions);
toc

% max(max(max(max(Policy))))<n_a % Double check that never try to leave top of asset grid.
% sum(sum(sum(sum(Policy==n_a))))

%% Initial distribution of agents at birth (j=1)
if gpuDeviceCount>0 % This is just so code can be run with or without gpu
    jequaloneDist=zeros(n_a,n_z,'gpuArray');
else
    jequaloneDist=zeros(n_a,n_z);
end
[trash,zeroassets_index]=min(abs(a_grid));
% Place them on the existing z_grid based on Params.sigmasqy1 (the variance
% of earnings at age 1) under assumption of normal distribution.
% Following lines implement this, they are largely just a copy of some code from the TauchenMethod() command.
sigma=sqrt(Params.sigmasqy1); %stddev of e
for ii=1:length(z_grid)
    jequaloneDist(zeroassets_index,1)=normcdf(z_grid(1)+(z_grid(2)-z_grid(1))/2,0,sigma);
    for jj=2:(length(z_grid)-1)
        jequaloneDist(zeroassets_index,jj)=normcdf(z_grid(jj)+(z_grid(jj+1)-z_grid(jj))/2,0,sigma)-normcdf(z_grid(jj)-(z_grid(jj)-z_grid(jj-1))/2,0,sigma);
    end
    jequaloneDist(zeroassets_index,end)=1-normcdf(z_grid(end)-(z_grid(end)-z_grid(end-1))/2,0,sigma);
end

%% Agents age distribution
% Almost all OLG models include some kind of population growth, and perhaps
% some other things that create a weighting of different ages that needs to
% be used to calculate the stationary distribution and aggregate variable.
Params.mewj=ones(1,Params.J);
for jj=2:length(Params.mewj)
    Params.mewj(jj)=Params.sj(jj-1)*Params.mewj(jj-1)/(1+Params.n);
end
Params.mewj=Params.mewj./sum(Params.mewj);
AgeWeightsParamNames={'mewj'}; % Many finite horizon models apply different weights to different 'ages'; eg., due to survival rates or population growth rates.

Params.fractionretired=sum(Params.mewj.*Params.bvec); % Note: bvec is really just an indicator of retirement
%% Test
disp('Test StationaryDist')
simoptions.parallel=3; % Sparse matrix on cpu (results moved back to gpu)
% simoptions=struct(); % Just use the defaults
tic;
StationaryDist=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightsParamNames,Policy,0,n_a,n_z,N_j,pi_z,Params,simoptions);
toc

%% Set up the General Equilibrium conditions (on assets/interest rate, assuming a representative firm with Cobb-Douglas production function)

% Functions to evaluate
FnsToEvaluate.K = @(aprime,a,z) a; % Aggregate assets (which is this periods state)
FnsToEvaluate.L = @(aprime,a,z,ybarj) exp(z+ybarj); % Aggregate labour supply (in efficiency units)
FnsToEvaluate.Beq = @(aprime,a,z,sj,r,tau) (1-sj)*aprime*(1+r*(1-tau)); % Total accidental bequests
% Note that the aggregate labour supply is actually entirely exogenous and so I could just precompute it, but am feeling lazy.

% General Equilibrium Equations
% Recall that GEPriceParamNames={'r','b','T'};
GeneralEqmEqns.capitalmarket = @(r,L,K,A,alpha,delta) r-(A*alpha*(K^(alpha-1))*(L^(1-alpha))-delta); % Rate of return on assets is related to Marginal Product of Capital
GeneralEqmEqns.SSbalance = @(b,K,L,theta,fractionretired,alpha,A) b*fractionretired-theta*(A*(1-alpha)*(K^(alpha))*(L^(-alpha)))*L; % Retirement benefits equal Payroll tax revenue: b*fractionretired-theta*w*L
GeneralEqmEqns.Bequests = @(T,Beq,n) T-Beq/(1+n); % Lump-sum transfers equal Accidental bequests 

%% Test
disp('Test AggVars')
tic;
AggVars=EvalFnOnAgentDist_AggVars_FHorz_Case1(StationaryDist, Policy, FnsToEvaluate, Params, [], 0, n_a, n_z,N_j, 0, a_grid, z_grid);
toc

%% Solve for the General Equilibrium
% Use the toolkit to find the equilibrium price index. There are two ways
% to do this. In what follows I use the 'search' approach to calculate the
% (initial) General equilibrium. But the commented-out lines that follow it
% show how to set up the grid approach.

% Without p_grid, just searching. Use n_p=0. (Setting the actual algorithm
% used to 'search' can be done with heteroagentoptions.fminalgo)
heteroagentoptions.verbose=1;
[p_eqm,p_eqm_index, GeneralEqmEqnsValues]=HeteroAgentStationaryEqm_Case1_FHorz(jequaloneDist,AgeWeightsParamNames,0, n_a, n_z, N_j, 0, pi_z, 0, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, [], [], [], GEPriceParamNames,heteroagentoptions,simoptions,vfoptions);
Params.r=p_eqm.r;
Params.b=p_eqm.b;
Params.T=p_eqm.T;
% save ./SavedOutput/Huggett1996.mat Params

% Using p_grid. This can be helpful if you want to, e.g., look for
% possibility of multiple equilibria.
% % GEPriceParamNames={'r','b','T'}; % Already delared above.
% r_grid=linspace(0.5,2,21)'*Params.r;
% b_grid=linspace(0.5,2,21)'*Params.b;
% T_grid=linspace(0.5,2,21)'*Params.T;
% p_grid=[r_grid,b_grid, T_grid];
% 
% disp('Calculating price vector corresponding to the stationary eqm')
% n_p=[length(r_grid),length(b_grid),length(T_grid)];
% heteroagentoptions.pgrid=p_grid;
% heteroagentoptions.verbose=1;
% [p_eqm,p_eqm_index, GeneralEqmEqnsValues]=HeteroAgentStationaryEqm_Case1_FHorz(jequaloneDist,AgeWeights,0, n_a, n_z, N_j, n_p, pi_z, 0, a_grid, z_grid, ReturnFn, SSvaluesFn, GeneralEqmEqns, Params, DiscountFactorParamNames, [], [], [], GEPriceParamNames, heteroagentoptions, simoptions, vfoptions);
% Params.r=p_eqm(1);
% Params.b=p_eqm(2);
% Params.T=p_eqm(3);
% save ./SavedOutput/Huggett1996grid.mat Params

%% Compute a few things about the equilibrium
[V, Policy]=ValueFnIter_Case1_FHorz(0,n_a,n_z,N_j, 0, a_grid, z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, [],vfoptions);
StationaryDist=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightsParamNames,Policy,0,n_a,n_z,N_j,pi_z,Params,simoptions);

% Start with some basics (these were previously being done inside return function):
% Rearranging that r=MPK-delta gives the following eqn (MPK is marginal product of capital)
KdivL=((Params.r+Params.delta)/(Params.alpha*Params.A))^(1/(Params.alpha-1));
% From K/Y, substituting the production function gives
KdivY=(KdivL^(1-Params.alpha))/Params.A;
% We know w=MPL (MPL is marginal product of labour)
Params.w=Params.A*(1-Params.alpha)*(KdivL^Params.alpha); % wage rate (per effective labour unit)
% Huggett (1996) calibrates tau to the following (see pg 478 for explanation)
Params.tau=0.195/(1-Params.delta*KdivY);

% Aggregate wealth transfers.
AggregateWealthTransers=zeros(1,N_j);
for jj=1:Params.J
    for ii=1:jj
        AggregateWealthTransers(jj)=AggregateWealthTransers(jj)+Params.T*(1+Params.r*(1-Params.tau))^(ii-1);
    end
end
AggregateWealthTransers=sum(Params.mewj.*AggregateWealthTransers);
% Total wealth
FnsToEvaluate.TotalWealth = @(aprime,a,z) a;
AggVars=EvalFnOnAgentDist_AggVars_FHorz_Case1(StationaryDist, Policy, FnsToEvaluate, Params, [], 0, n_a, n_z,N_j, 0, a_grid, z_grid);
% Transfer Wealth Ratio
TransferWealthRatio=AggregateWealthTransers/AggVars.TotalWealth.Mean;


% Calculate fraction of population with zero or negative wealth
FnsToEvaluate.ZeroOrNegAssets = @(aprime,a,z) (a<=0); % Indicator for ZeroOrNegAssets
AggVars=EvalFnOnAgentDist_AggVars_FHorz_Case1(StationaryDist, Policy, FnsToEvaluate, Params, [], 0, n_a, n_z,N_j, 0, a_grid, z_grid);
FractionWithZeroOrNegAssets=100*AggVars.ZeroOrNegAssets.Mean;

% Calculate wealth lorenz curve (and thus all percentile shares) and also
% the that for earnings (in text at bottom of pg 480, top of pg 481, there
% are a bunch of descriptions of model earnings, conditional on working age)
FnsToEvaluate.Wealth = @(aprime,a,z) a; % Notice that wealth is just the same as aggregate assets. I am going to evaluate it again anyway but this is total overkill.
AllStats=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDist, Policy, FnsToEvaluate, Params, 0, n_a, n_z,N_j, 0, a_grid, z_grid);
TopWealthShares=100*(1-AllStats.Wealth.LorenzCurve([80,95,99],1)); % Need the 20,5, and 1 top shares for Tables of Huggett (1996)
% Calculate the wealth gini
WealthGini=Gini_from_LorenzCurve(AllStats.Wealth.LorenzCurve);

AgeConditionalStats=LifeCycleProfiles_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate,Params,[],0,n_a,n_z,N_j,0,a_grid,z_grid);

%% Draw figures from Huggett (1996)

% Fig 1
figure(1);
plot(19+(1:1:Params.J),exp(Params.ybarj)./sum(exp(Params.ybarj).*Params.mewj))
title({'Earnings Profile (ratio to overall mean)'})
xlabel('Age')
ylabel('Earnings')
% Fig 2
figure(2);
plot(19+(1:1:71),AgeConditionalStats.Wealth.Mean(1:end-8),19+(1:1:71),AgeConditionalStats.Wealth.QuantileCutoffs(2,1:end-8),19+(1:1:71),AgeConditionalStats.Wealth.QuantileCutoffs(5,1:end-8),19+(1:1:71),AgeConditionalStats.Wealth.QuantileCutoffs(10,1:end-8))
legend('Mean','10th','25th','Median')
title({'Wealth Profiles: Uncertain Lifetimes'})
xlabel('Age')
ylabel('Wealth')
% Fig 4(ish)
figure(3);
plot(19+(6:1:71),AgeConditionalStats.Wealth.Gini(6:(end-8)))
title('Wealth Gini at different ages')
xlabel('Age')
ylabel('Wealth Gini')


%% Give some of the outputs that Huggett (1996) reports

fprintf('The capital-output ratio, K/Y is: %8.2f \n', KdivY);

fprintf('The transfer-wealth ratio is: %8.2f \n', TransferWealthRatio);

fprintf('The Wealth Gini is: %8.2f \n', WealthGini);

fprintf('The percentage share of total wealth held by top 1 percent of wealth-holders is: %8.2f \n', TopWealthShares(3));
fprintf('The percentage share of total wealth held by top 5 percent of wealth-holders is: %8.2f \n', TopWealthShares(2));
fprintf('The percentage share of total wealth held by top 20 percent of wealth-holders is: %8.2f \n', TopWealthShares(1));

fprintf('The percentage of population with zero or negative wealth is: %8.2f \n', FractionWithZeroOrNegAssets);

%% A double-check, calculate some statistics about earnings (in text at bottom of pg 480, top of pg 481, there
% are a bunch of descriptions of model earnings, conditional on working age). Since earnings inequality is completely exogenous (earnings depends
% endogenously on w, but this will disappear from all the relative inequality measures).

FnsToEvaluate.Earnings = @(aprime,a,z,w,ybarj) w*exp(z+ybarj); % Earnings

% options.agegroupings=1:1:N_j; % for each age, this is anyway the default
AgeConditionalStats=LifeCycleProfiles_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate,Params,[],0,n_a,n_z,N_j,0,a_grid,z_grid);
simoptions.agegroupings=[1,Params.JR]; % Working age, Retired (not actually interested in the numbers for retired)
AllEmployedStats=LifeCycleProfiles_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate,Params,[],0,n_a,n_z,N_j,0,a_grid,z_grid,simoptions);

fprintf('The Earnings Gini for 20 year olds is: %8.2f \n', AgeConditionalStats.Earnings.Gini(1));
fprintf('The Earnings Gini for 65 year olds is: %8.2f \n', AgeConditionalStats.Earnings.Gini(Params.JR));
fprintf('The Earnings Gini for Working Age Population is: %8.2f \n', AllEmployedStats.Earnings.Gini(1));
fprintf('The percentage share of earnings going to top 1 percent of earners is: %8.2f \n', 100*(1-AllEmployedStats.Earnings.LorenzCurve(99,1)));
fprintf('The percentage share of earnings going to top 5 percent of earners is: %8.2f \n', 100*(1-AllEmployedStats.Earnings.LorenzCurve(95,1)));
fprintf('The percentage share of earnings going to top 10 percent of earners is: %8.2f \n', 100*(1-AllEmployedStats.Earnings.LorenzCurve(90,1)));
fprintf('The percentage share of earnings going to top 20 percent of earners is: %8.2f \n', 100*(1-AllEmployedStats.Earnings.LorenzCurve(80,1)));

% %% Dean Corbae's Homework 5 inputs:
% % According to Dean Corbae's HW5 files he got the original numbers from
% % Huggett and they should be
% Corbae_z_grid = [0.0849 0.1156 0.1573 0.2141 0.2915 0.3967 0.5399 0.7348 1.0000 1.3610 1.8523 2.5210 3.4311 4.6697 6.3555 8.6499 11.7725 40.3927]; % F06Hzvec.txt
% % Compare this to exp(z_grid)
% Corbae_jinitaldist=[0.0001 0.0005 0.0024 0.0092 0.0278 0.0656 0.1210 0.1747 0.1974 0.1747 0.1210 0.0656 0.0278 0.0092 0.0024 0.0005 0.0001 0.0000]; % F06Hz1probvec.txt
% % Compare this to jequaloneDist(zeroassets_index,:)
% Corbae_deterministicincome = [0.0911 0.1573 0.2268 0.2752 0.3218 0.3669 0.4114 0.4559 0.4859 0.5164 0.5474 0.5786 0.6097 0.6311 0.6517 0.6711 0.6893 0.7060 0.7213 0.7355 0.7489 0.7619 0.7747 0.7783 0.7825 0.7874 0.7931 0.7994 0.7923 0.7850 0.7771 0.7679 0.7567 0.7351 0.7105 0.6822 0.6500 0.6138 0.5675 0.5183 0.4672 0.3935 0.3239 0.2596 0.1955 0.1408 0.0959 0.0604 0.0459 0.0342 0.0246 0.0165 0.0091 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]; % F06Hyvec.txt
% % "Finally, assume that the normalized (by the population) labor efficiency hours are given by L = Lt/Nt = 0.5289"
% % Compare this to exp(Params.ybarj)
% % 0.5289*exp(Params.ybarj)
% plot(1:1:79, Corbae_deterministicincome, 1:1:79, exp(Params.ybarj))
% 
% Corbae_pi_z=[0.6031555	0.3536218	0.042449579	0.00076694961	0.0000019141335	0.000000000626149	0.000000000000026	0	0	0	0	0	0	0	0	0	0	0;...
%     0.128545	0.496843	0.3364697	0.037517324	0.00062712817	0.0000014441507	0.000000000435219	0.000000000000017	0	0	0	0	0	0	0	0	0	0;...
%     0.0057332767	0.1354109	0.5060558	0.3192306	0.033057541	0.0005111752	0.0000010860242	0.000000000301517	0.000000000000011	0	0	0	0	0	0	0	0	0;...
%     0.000043726395	0.0067113182	0.1478038	0.5139832	0.3020029	0.029039346	0.00041533652	0.00000081404778	0.000000000208198	0.000000000000007	0	0	0	0	0	0	0	0;...
%     0.000000052497327	0.000055510805	0.0078789536	0.1608557	0.5205598	0.2848815	0.025431748	0.00033639232	0.00000060819696	0.000000000143287	0.000000000000004	0	0	0	0	0	0	0;...
%     0.000000000009551	0.000000072109074	0.000070316441	0.009221036	0.1745456	0.5257324	0.2679545	0.022204313	0.00027158641	0.00000045291998	0.000000000098289	0.000000000000003	0	0	0	0	0	0;...
%     0	0.000000000014203	0.000000098749716	0.000088785499	0.010758306	0.1888452	0.5294583	0.2513036	0.019327182	0.00021856703	0.00000033618835	0.000000000067199	0.000000000000002	0	0	0	0	0;...
%     0	0	0.000000000021049	0.0000001347899	0.00011174595	0.012513041	0.2037184	0.531706	0.2350039	0.016771287	0.00017533687	0.00000024872674	0.000000000045792	0.000000000000001	0	0	0	0;...
%     0	0	0.000000000000001	0.000000000031094	0.00000018338343	0.00014019408	0.014509136	0.2191215	0.5324572	0.2191228	0.014508678	0.00014020711	0.0000001834176	0.000000000031101	0.000000000000001	0	0	0;...
%     0	0	0	0.000000000000001	0.000000000045781	0.00000024868015	0.00017532188	0.016771846	0.2350028	0.5317057	0.2037196	0.012512634	0.00011175699	0.00000013481558	0.000000000021054	0	0	0;...
%     0	0	0	0	0.000000000000002	0.000000000067183	0.00000033612619	0.00021855129	0.019327847	0.2513023	0.529458	0.1888461	0.010757899	0.00008879503	0.000000098768027	0.000000000014205	0	0;...
%     0	0	0	0	0	0.000000000000003	0.000000000098266	0.00000045283821	0.00027156886	0.022205003	0.2679527	0.5257323	0.1745464	0.0092207119	0.000070324822	0.00000007212293	0.000000000009553	0;...
%     0	0	0	0	0	0	0.000000000000004	0.000000000143254	0.00000060808856	0.0003363732	0.025432453	0.28488	0.5205595	0.1608566	0.0078786584	0.000055517812	0.000000052500724	0;...
%     0	0	0	0	0	0	0	0.000000000000007	0.00000000020815	0.00000081390453	0.0004153155	0.029040076	0.3020017	0.5139828	0.1478043	0.006711069	0.000043726046	0.000000000000021;...
%     0	0	0	0	0	0	0	0	0.000000000000011	0.000000000301446	0.0000010858275	0.0005111529	0.033058293	0.3192292	0.5060557	0.1354114	0.0057331258	0.000000000362423;...
%     0	0	0	0	0	0	0	0	0	0.000000000000017	0.000000000435124	0.0000014439051	0.00062711048	0.037518173	0.3364691	0.4968424	0.1285409	0.00000094089796;...
%     0	0	0	0	0	0	0	0	0	0	0.000000000000026	0.000000000626009	0.0000019138029	0.00076692866	0.042450383	0.3536209	0.6027857	0.00037465594;...
%     0	0	0	0	0	0	0	0	0	0	0	0	0	0	0.000000000000151	0.000000002594984	0.013603883	0.9863952];
% % Compare this to pi_z
% 
% % I have exact same grid on z
% % My pi_z is same to order of 10^-5
% % I have copied the deterministic income spline from Corbae's HW5 numbers.
