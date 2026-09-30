% Example just solving the baseline model of Bruggemann (2021) - Higher Taxes at the Top: The Role of Entrepreneurs
%
% assets: assets
% l: labor supply
% e: entrepreneur/worker
% age: young/old
% eta: labor productivity (B2021 calls this epsilon)
% theta: entrepreneurial ability
%
% Just solves the baseline stationary general equilibrium and a few of the Table 3 statistics


%% Model size
% Two decision variables (you can change n_asset and n_l, the rest are hardcoded by the grids below)
n_l=150; % labor supply
n_entre=2; % 0=worker, 1=entrepreneur
% We only need n_l for the workers (as entrepreneur has exogenous fixed
% labor supply) so we use a joint-grid for the decision variables.
% One endogenous state
n_asset=501; % assets
% Three markov exogenous states
n_age=2; % age (1=young & 2=old)
n_eta=6; % eta (labor productivity)
n_theta=4; % theta (entrepreneurial ability)
% But importantly, we only care about eta & theta when young, we do not need them when old.
% We therefore set this up as a joint-grid, so there is one single 'old' state rather than
% n_eta*n_theta copies of it that would all behave identically.

%% Parameters

% Preferences
Params.beta=0.9; % discount factor
Params.sigma1=1.5; % CRRA risk aversion
Params.sigma2=1.7; % Inverse of Frisch elasticity of labor supply
Params.ell=3; % household time endowment (not really used, see creation of l_grid below)
Params.xi=0.716; % weight on disutility of labor supply

% Stochastic ageing
Params.pi_y=0.978; % probability of remaining young
Params.pi_o=0.911; % probability of remaining old

% Entrepreneur's production fn
Params.gamma=0.359; % relative importance of capital
Params.upsilon=0.864; % (decreasing) returns to scale [Lucas span-of-control parameter]
% Collateral constraint
Params.lambda=1.5;
% Entrepreneurs labor input
Params.lbar=1;

% Corporate sector production fn
Params.Z=1; % technology level (normalization; B2021 calls this A)
Params.alpha=0.33; % relative importance of capital

% Depreciation rate
Params.delta=0.06; % depreciation rate of capital (in both sectors)

% Taxes
Params.ybar=1; % average income (initial guess, will be determined in eqm)
Params.tau_c=0.110; % consumption tax rate
Params.d=0.2385; % (times ybar) standard deduction (from income, before income taxes)
Params.tau_i_r1_stat=0.1; % Statutory Income tax rate, bracket 1
Params.tau_i_r2_stat=0.15; % Statutory Income tax rate, bracket 2
Params.tau_i_r3_stat=0.25; % Statutory Income tax rate, bracket 3
Params.tau_i_r4_stat=0.28; % Statutory Income tax rate, bracket 4
Params.tau_i_r5_stat=0.33; % Statutory Income tax rate, bracket 5
Params.tau_i_r6_stat=0.35; % Statutory Income tax rate, bracket 6
Params.tau_i_t1=0; % Income tax threshold 1
Params.tau_i_t2=0.214; % (times ybar) % Income tax threshold 2
Params.tau_i_t3=0.868; % (times ybar) % Income tax threshold 3
Params.tau_i_t4=1.753; % (times ybar) % Income tax threshold 4
Params.tau_i_t5=2.672; % (times ybar) % Income tax threshold 5
Params.tau_i_t6=4.771; % (times ybar) % Income tax threshold 6
Params.tau_i_adj=0.669; % linear scaling factor so that income tax raises the 'right' revenue
Params.tau_s=0; % flat-tax on income representing state and local taxes (initial guess)

% Put the adjustment onto the statutory rates to get the ones the model uses
Params.tau_i_r1=Params.tau_i_adj*Params.tau_i_r1_stat; % Income tax rate, bracket 1
Params.tau_i_r2=Params.tau_i_adj*Params.tau_i_r2_stat; % Income tax rate, bracket 2
Params.tau_i_r3=Params.tau_i_adj*Params.tau_i_r3_stat; % Income tax rate, bracket 3
Params.tau_i_r4=Params.tau_i_adj*Params.tau_i_r4_stat; % Income tax rate, bracket 4
Params.tau_i_r5=Params.tau_i_adj*Params.tau_i_r5_stat; % Income tax rate, bracket 5
Params.tau_i_r6=Params.tau_i_adj*Params.tau_i_r6_stat; % Income tax rate, bracket 6
% Note: the effective top marginal tax rate is tau_i_adj times the statutory one. It is worth
% keeping track of which of the two any given number refers to, as the two differ by a third.

% Prices (initial guesses; these are general eqm parameters)
Params.r=0.1; % interest rate
Params.w=0.7; % wage per-productivity-unit-per-unit-of-time

% Government spending
Params.G=0.354; % government spending
Params.pension=0.4; % pension (B2021 calls this b) (initial guess, will be determined in eqm)

% Lump-Sum transfers
Params.lumpsum=0; % 0 in baseline, but B2021 uses it elsewhere to return any additional tax revenues to households

%% Grids

% Labor supply, l
Params.maxl=1.6;
l_grid=linspace(0,Params.maxl,n_l)'; % labor supply
% make sure lbar is a point in the grid [lbar is the fixed labor supply that entrepreneurs must provide]
[~,lbarindex]=min(abs(l_grid-Params.lbar));
l_grid(lbarindex)=Params.lbar;

entre_grid=[0;1]; % 0=worker, 1=entrepreneur

% Set grid for asset holdings
assetmaxfactor=520; % This is the max assets
asset_grid=assetmaxfactor*(linspace(0,1,n_asset).^3)'; % linspace ^3 puts more points near zero, where the curvature of value and policy functions is higher

age_grid=[1;2]; % 1=young, 2=old

pi_age=[Params.pi_y, 1-Params.pi_y;...
    1-Params.pi_o, Params.pi_o]; % transitions between young and old

% B2021 pg 12 "I take the values for the first five levels of the labor ability process from Cagetti & De Nardi (2009)
% ...[and] introduce a high sixth level." [CDN2009 report the grid and transition probabilities in their Appendix A]
eta_grid1to5=[0.2468, 0.4473, 0.7654, 1.3097, 2.3742]';
pi_eta1to5=[0.7376, 0.2473, 0.0150, 0.0002, 0.0000;....
    0.1947, 0.5555, 0.2328, 0.0169, 0.0001;...
    0.0113, 0.2221, 0.5333, 0.2221, 0.0113;...
    0.0001 0.0169 0.2328 0.5555 0.1947;...
    0.0000 0.0002 0.0150 0.2473 0.7376];
% Now add in the sixth point
Params.eta6=26; % value of 6th eta point
eta_grid=[eta_grid1to5; Params.eta6];
Params.prob_eta6=0.00160; % probability of going to 6th eta point [from Appendix, you cannot use rounded version in the main paper]
Params.prob_eta6toeta3=0.071; % from 6th eta point, you either go to 3rd point, or remain in 6th
pi_eta=[pi_eta1to5*(1-Params.prob_eta6), Params.prob_eta6*ones(5,1);...
    0,0,Params.prob_eta6toeta3,0,0,1-Params.prob_eta6toeta3];
pi_eta=pi_eta./sum(pi_eta,2); % renormalize rows (some were 1.0001)

% B2021 Table 2 reports theta grid and transition probabilities
theta_grid=[0,0.682,1.750,2.818]'; % entrepreneurial ability
pi_theta=[0.963, 0.037, 0, 0;...
    0.275, 0.581, 0.144, 0;...
    0, 0.275, 0.581, 0.144;...
    0, 0, 0.304, 0.696]; % B2021 calls this Lambda

% We need the stationary distributions for both eta and theta, as newborns are
% assumed to draw their eta and theta from these
[eta_mean,~,~,eta_statdist]=MarkovChainMoments(eta_grid,pi_eta);
[theta_mean,~,~,theta_statdist]=MarkovChainMoments(theta_grid,pi_theta);

%% Get into form for VFI toolkit

% Joint grid on the decision variables. A cross-product grid would be wasteful here because entrepreneurs have no
% labor supply choice: they must supply lbar. So the joint grid is the n_l worker rows,
% each with its own l and e=0, plus one single entrepreneur row with l=lbar and e=1.
n_d=[n_l+1,1]; % hardcodes n_entre=2
if n_entre~=2
    error('joint grid for decision variable hardcodes n_entre=2')
end
d_grid=[[l_grid; l_grid(lbarindex)],[entre_grid(1)*ones(n_l,1); entre_grid(2)]];

n_a=n_asset;

% Joint grid on the exogenous states, n_eta*n_theta points for young plus one point for old.
N_young=n_eta*n_theta; % number of young exogenous states
n_z=[N_young+1,1,1];
a_grid=asset_grid;
% For z_grid, the first N_young rows are young, the last row is old
z_grid=[age_grid(1)*ones(N_young,1), repmat(eta_grid,n_theta,1), repelem(theta_grid,n_eta,1);
        age_grid(2),0,0];
% note: use zero for eta and theta when old, so you will get issues if you try to do anything with them

% When agents are 'old' it becomes irrelevant what their eta & theta are, so we use a
% joint-grid on z and do not track them. When agents become 'young' the "new born household's
% two abilities are uncorrelated with the abilities of the parent household", so use
% pi_newborn based on the stationary distributions of eta and theta (which are i.i.d).
% When agents are 'young', the transitions are based on pi_eta and pi_theta created above.
pi_etatheta=kron(pi_theta, pi_eta); % in reverse order [young-young]

etatheta_statdist=kron(theta_statdist,eta_statdist); % in reverse order
pi_newborn=etatheta_statdist'; % i.i.d. based on etatheta_statdist

pi_z=[pi_age(1,1)*pi_etatheta, pi_age(1,2)*ones(N_young,1);...  % top left is young-young, top right is young-old,
    pi_age(2,1)*pi_newborn, pi_age(2,2)]; % bottom left is old-young, bottom right is old-old

%% Return fn
DiscountFactorParamNames={'beta'};

ReturnFn=@(l,e,aprime,a,age,eta,theta,r,w,sigma1,sigma2,xi,lbar,lambda,delta,gamma,upsilon,lumpsum,pension,tau_c,d,tau_s,ybar, tau_i_r1,tau_i_r2,tau_i_r3,tau_i_r4,tau_i_r5,tau_i_r6, tau_i_t1, tau_i_t2, tau_i_t3, tau_i_t4, tau_i_t5, tau_i_t6)...
    Bruggemann2021_ReturnFn(l,e,aprime,a,age,eta,theta,r,w,sigma1,sigma2,xi,lbar,lambda,delta,gamma,upsilon,lumpsum,pension,tau_c,d,tau_s,ybar, tau_i_r1,tau_i_r2,tau_i_r3,tau_i_r4,tau_i_r5,tau_i_r6, tau_i_t1, tau_i_t2, tau_i_t3, tau_i_t4, tau_i_t5, tau_i_t6);
% The first inputs must be: decision variables, next period endogenous state, endogenous state, exogenous state. Followed by any parameters

%% Aggregates

% Create functions to be evaluated
FnsToEvaluate.K_noncorp = @(l,e,aprime,a,age,eta,theta, r,w,lambda,delta,gamma,upsilon,lbar)...
    B2021_kFn(l,e,aprime,a,age,eta,theta,r,w,lambda,delta,gamma,upsilon,lbar); % Assets used in non-corporate sector (=entrepreneurs)
FnsToEvaluate.A = @(l,e,aprime,a,age,eta,theta) a; % Total assets of households (workers and entrepreneurs)
FnsToEvaluate.N_noncorp = @(l,e,aprime,a,age,eta,theta,r,w,lambda,delta,gamma,upsilon,lbar)...
    B2021_nFn(l,e,aprime,a,age,eta,theta,r,w,lambda,delta,gamma,upsilon,lbar); % Labor used in non-corporate sector (=entrepreneurs), excluding own labor
FnsToEvaluate.N_lbar = @(l,e,aprime,a,age,eta,theta,lbar) lbar*e*(age==1); % Entrepreneurs own labor supply
FnsToEvaluate.L = @(l,e,aprime,a,age,eta,theta) l*eta*(e==0)+l*(e==1); % Total labor supply
FnsToEvaluate.Y_noncorp =  @(l,e,aprime,a,age,eta,theta,r,w,lambda,delta,gamma,upsilon,lbar)...
    B2021_YnoncorpFn(l,e,aprime,a,age,eta,theta,r,w,lambda,delta,gamma,upsilon,lbar); % output of non-corporate sector
FnsToEvaluate.IncomeTaxRevenue = @(l,e,aprime,a,age,eta,theta,r,w,lambda,delta,gamma,upsilon,lbar,pension,d,tau_s,ybar, tau_i_r1,tau_i_r2,tau_i_r3,tau_i_r4,tau_i_r5,tau_i_r6, tau_i_t1, tau_i_t2, tau_i_t3, tau_i_t4, tau_i_t5, tau_i_t6)...
    B2021_IncomeTaxRevenueFn(l,e,aprime,a,age,eta,theta,r,w,lambda,delta,gamma,upsilon,lbar,pension,d,tau_s,ybar, tau_i_r1,tau_i_r2,tau_i_r3,tau_i_r4,tau_i_r5,tau_i_r6, tau_i_t1, tau_i_t2, tau_i_t3, tau_i_t4, tau_i_t5, tau_i_t6); % Tax Revenue from the Income Tax
FnsToEvaluate.ConsumptionTaxRevenue = @(l,e,aprime,a,age,eta,theta,r,w,lambda,delta,gamma,upsilon,lbar,lumpsum,pension,tau_c,d,tau_s,ybar, tau_i_r1,tau_i_r2,tau_i_r3,tau_i_r4,tau_i_r5,tau_i_r6, tau_i_t1, tau_i_t2, tau_i_t3, tau_i_t4, tau_i_t5, tau_i_t6)...
    tau_c*B2021_ConsumptionFn(l,e,aprime,a,age,eta,theta,r,w,lambda,delta,gamma,upsilon,lbar,lumpsum,pension,tau_c,d,tau_s,ybar, tau_i_r1,tau_i_r2,tau_i_r3,tau_i_r4,tau_i_r5,tau_i_r6, tau_i_t1, tau_i_t2, tau_i_t3, tau_i_t4, tau_i_t5, tau_i_t6); % Tax Revenue from the Consumption Tax
FnsToEvaluate.PensionSpending=@(l,e,aprime,a,age,eta,theta,pension) pension*(age==2);
FnsToEvaluate.Income=@(l,e,aprime,a,age,eta,theta,r,w,lambda,delta,gamma,upsilon,lbar,pension)...
    B2021_IncomeFn(l,e,aprime,a,age,eta,theta,r,w,lambda,delta,gamma,upsilon,lbar,pension); % taxable income, before deductions
FnsToEvaluate.Entrepreneur =  @(l,e,aprime,a,age,eta,theta) (e==1)*(age==1); % fraction of entrepreneurs
FnsToEvaluate.Young =  @(l,e,aprime,a,age,eta,theta) (age==1); % fraction young

%% General equilbrium

% IntermediateEqns. Note that corporate labor nets out BOTH the labor entrepreneurs hire (N_noncorp)
% and their own lbar (N_lbar, which is already an input inside Y_noncorp via their production fn).
heteroagentoptions.intermediateEqns.N_corp=@(L,N_noncorp,N_lbar) L-(N_noncorp+N_lbar);
heteroagentoptions.intermediateEqns.K_corp=@(A,K_noncorp) A-K_noncorp;
heteroagentoptions.intermediateEqns.Y_corp=@(K_corp,N_corp,alpha,Z) Z*(K_corp^alpha)*(N_corp^(1-alpha));
heteroagentoptions.intermediateEqns.Y=@(Y_corp,Y_noncorp) Y_corp+Y_noncorp;

% GE params
% r, w, and tau_s (which balances the government budget) are the three competitive-equilibrium prices.
% Three more are calibration targets:
%   ybar is set to the mean income
%   pension is set to target the pension replacement rate
% Note that pension is calibrated as a ratio, not levels.
GEPriceParamNames={'r','w','tau_s','pension','ybar'}; % No need to actually have w here (because it is Cobb-Douglas prodn fn in corporate sector, and we can therefore derive a relation between w and r that must hold in stationary eqm (but I don't bother here)

GeneralEqmEqns.CapitalMarket = @(r,K_corp,N_corp,alpha,delta,Z) r-(alpha*Z*(K_corp/N_corp)^(alpha-1)-delta); %The requirement that the interest rate corresponds to the marginal product of capital in corporate sector
GeneralEqmEqns.LaborMarket = @(w,K_corp,N_corp,alpha,Z) w-(1-alpha)*Z*(K_corp/N_corp)^(alpha); %The requirement that the wage corresponds to the marginal product of labor in corporate sector
GeneralEqmEqns.GovBudget = @(G,PensionSpending,IncomeTaxRevenue,ConsumptionTaxRevenue,lumpsum) G+PensionSpending+lumpsum-IncomeTaxRevenue-ConsumptionTaxRevenue; %Government runs balanced budget [take advantage of the fact that lumpsum is same for all, so adding up across everyone just gives the lumpsum parameter value; note, lumpsum=0 in baseline]
% Two calibration targets
GeneralEqmEqns.Pensions = @(pension,ybar,pensionreplacementrate) pension/ybar - pensionreplacementrate;
GeneralEqmEqns.AvgIncome = @(ybar,Income) ybar - Income; % get ybar to be the average income

%%
vfoptions.gridinterplayer=1; % set to 0 for debugging
vfoptions.ngridinterp=50;
vfoptions.maxaprimediff=20;
vfoptions.lowmemory=1;

simoptions.gridinterplayer=vfoptions.gridinterplayer;
simoptions.ngridinterp=vfoptions.ngridinterp;

heteroagentoptions.verbose=1; % verbose means that you want it to give you feedback on what is going on
heteroagentoptions.fminalgo=[8,1]; % fast but not really high accuracy, then a higher accuracy
heteroagentoptions.toleranceGEcondns=[1e-4,1e-5]; % high accuracy on final solve
% Note: see the note in Bruggemann2021.m on which MATLAB versions fminalgo=8 (lsqnonlin) needs

% Set initial value for general eqm
% (these are decent initial guesses based on the solution; the comments afterwards show the initial guess used the first time this was run)
Params.r=0.02; % 0.05
Params.w=1.3; % 1.2
Params.tau_s=0.1; % 0.01
Params.pension = 0.1; %0.4
Params.ybar    = 1; % 1


%% Solve the stationary general eqm
GEPriceParamNames={'r','w','tau_s','pension','ybar'};
heteroagentoptions.constrainpositive={'w'};

[p_eqm,GeneralEqmCondn]=HeteroAgentStationaryEqm_InfHorz(n_d, n_a, n_z, 0, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, [], [], [], GEPriceParamNames,heteroagentoptions, simoptions, vfoptions);

p_eqm % The equilibrium values of the GE prices

Params.r=p_eqm.r;
Params.w=p_eqm.w;
Params.tau_s=p_eqm.tau_s;

%% Now that we have the GE, let's calculate a bunch of related objects

[V,Policy]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, [], vfoptions);

PolicyVals=PolicyInd2Val_InfHorz(Policy,n_d,n_a,n_z,d_grid,a_grid, vfoptions); % This will give you the policy in terms of values rather than index

StationaryDist=StationaryDist_InfHorz(Policy,n_d,n_a,n_z,pi_z, simoptions);

%% Check that we don't hit the top of asset grid (this is not a figure from Bruggemann 2021, just something I want to see)
Fig13=figure(13);
assetdist_young_e=cumsum(sum(StationaryDist(:,1:N_young).*shiftdim(PolicyVals(2,:,1:N_young),1),2),1); % recall: joint-grid on z, and e is the second decision variable
assetdist_young_note=cumsum(sum(StationaryDist(:,1:N_young).*shiftdim(1-PolicyVals(2,:,1:N_young),1),2),1);
assetdist_old=cumsum(StationaryDist(:,N_young+1),1);
plot(asset_grid,assetdist_young_note,asset_grid,assetdist_young_e,asset_grid,assetdist_old)
title('cdf of HHs over assets')
legend('worker','entrepreneur','retiree')

%% Calculate various statistics related to the eqm, a few of which are in Table 3
FnsToEvaluate.EntrepreneurHire = @(l,e,aprime,a,age,eta,theta,r,w,lambda,delta,gamma,upsilon,lbar) (B2021_nFn(l,e,aprime,a,age,eta,theta,r,w,lambda,delta,gamma,upsilon,lbar)>0); % n>0 (note: n excludes lbar)
FnsToEvaluate.HoursWorked=@(l,e,aprime,a,age,eta,theta,lbar) l*(1-e)*(age==1)+lbar*e*(age==1); % l for workers and lbar for entrepreneurs

simoptions.conditionalrestrictions.entrepreneurs=@(l,e,aprime,a,age,eta,theta) e*(age==1); % entrepreneurs
simoptions.conditionalrestrictions.workers=@(l,e,aprime,a,age,eta,theta) (1-e)*(age==1); % workers
simoptions.conditionalrestrictions.workerretirees=@(l,e,aprime,a,age,eta,theta) (age==2) + (age==1)*(1-e); % workers and retirees
% Note: B2021's own do-files classify a household as an entrepreneur by kstar>0, which puts the
% retirees in with the workers. That is what the workerretirees restriction is for.

simoptions.npoints=100; % 100 points for the lorenz curve, so we can read off the top 1 percent income share
AllStats=EvalFnOnAgentDist_AllStats_InfHorz(StationaryDist, Policy, FnsToEvaluate,Params, [],n_d, n_a, n_z, d_grid, a_grid,z_grid,simoptions);
simoptions=rmfield(simoptions,'conditionalrestrictions');

% Just give some feedback so we can see nothing looks odd
AllStats.Entrepreneur.Mean
[AllStats.A.Mean, AllStats.L.Mean]
[AllStats.entrepreneurs.A.Mean, AllStats.workers.A.Mean]
[AllStats.A.Mean, AllStats.K_noncorp.Mean, AllStats.A.Mean-AllStats.K_noncorp.Mean] % assets, capital to entrepreneurs, capital to corporate
[AllStats.L.Mean, AllStats.N_noncorp.Mean, AllStats.L.Mean-AllStats.N_lbar.Mean-AllStats.N_noncorp.Mean] % labor supply, labor to entrepreneurs, labor to corporate

Output_corp=Params.Z*((AllStats.A.Mean-AllStats.K_noncorp.Mean)^Params.alpha)*((AllStats.L.Mean-AllStats.N_lbar.Mean-AllStats.N_noncorp.Mean)^(1-Params.alpha));
Y=Output_corp+AllStats.Y_noncorp.Mean;
[Y,AllStats.Y_noncorp.Mean, Output_corp]

%% A handful of the Table 3 statistics
fprintf('\n')
fprintf('Capital-output ratio               %6.2f \n', AllStats.A.Mean/Y)
fprintf('Top 1 percent income share         %6.2f \n', 1-AllStats.Income.LorenzCurve(99))
fprintf('Wealth Gini                        %6.2f \n', AllStats.A.Gini)
fprintf('Fraction of entrepreneurs          %6.2f \n', AllStats.Entrepreneur.Mean)
fprintf('Entrepreneurs income Gini          %6.2f \n', AllStats.entrepreneurs.Income.Gini)
fprintf('Share of entrepreneurs who hire    %6.2f \n', AllStats.entrepreneurs.EntrepreneurHire.Mean)
fprintf('Average working time               %6.2f \n', AllStats.workers.HoursWorked.Mean)
fprintf('Workers income Gini                %6.2f \n', AllStats.workerretirees.Income.Gini)
fprintf('\n')

