% Example based on Imrohoroglu, Imrohoroglu & Joines (2003) - Time-Inconsistent Preferences and Social Security
% Demonstrates use of quasi-hyperbolic preferences.

% NEED TO FIX epsilon_j

% The quasi-hyperbolic preferences are controlled using,
% vfoptions.exoticpreferences='QuasiHyperbolic' % Use the quasi-hyperbolic preferences
% Also need to choose which of the two quasi-hyperbolic solutions, naive or sophisticated, to use.
% vfoptions.quasi_hyperbolic='Naive' % This is the default, alternative is 'Sophisticated'.
% These are set below.

% One decision variable: labor supply
% Two endogenous state variable: assets (k), retirement state (b)
% One stochastic exogenous state variable (z, which is either 'employed' or 'unemployed')
% Age

Params.J=65; % Number of period in life-cycle (represents ages 21 to 85)

% Grid sizes to use
n_l=101; % labour supply
n_k=501; % Assets (IIJ2003 call this a, rather than k); IIJ2003 use 4097 points
n_ebar=201; % Average Indexed Monthly Earnings (determines social security benefits; in model it is actuall annual earnings)
n_b=2; % Retired or not retired
n_z=2; % employed or unemployed
N_j=Params.J; % Number of periods in finite horizon

% IIJ2003 discuss 'Effect 1' and 'Effect 2'. Quasi-hyperbolic preferences
% is 'Effect 2'. Since 'Effect 1' does not change the actual policy
% functions it is observationally equivalent to time-consistend preferences (leads to same policy fn)
% with exponential discounting and so other than philosophically is of little to no practical/scientific interest.
% I therefore do not consider 'Effect 1' here.

%% Set model parameters

%Preference parameters
Params.gamma=2;
% IIJ2003 consider three values for gamma: 1,2,3
Params.beta0= 0.85; % The 'first' quasi-hyperbolic discounting parameter (IIJ2003 call this beta)
% IIJ2003 consider three values for beta0: 0.6, 0.85, 0.9
Params.beta=0.95; % The 'second' quasi-hyperbolic discounting parameter (IIJ2003 call this delta_f)
% Just an initial guess, beta is determined to make a capital-output ratio of 2.52
Params.varphi=0.33; % share of consumption in the utility fn

%Technology parameters
Params.A=1.76; % TFP, set to normalize output to 1 (IIJ2003 call this B)
% Just an initial guess, A is determined to normalize output to 1
Params.alpha=0.31; % Capital-share in Cobb-Douglas production fn (IIJ2003 use alpha as labor-share, so my value is 1-alpha in their notation)
Params.delta=0.044; % Depreciation rate (IIJ2003 call this d)
Params.g=Params.alpha*0.0165; % Rate of growth of producivity (IIJ2003 pg 757, TFP growth is alpha*rho, IIJ2003 pg 764)

% People become economically active (j=1) at age 21, retire at 65, and max age is 85
% Params.J=N_j; %=85-21
Params.Jr=45; % =65-21 % Age at which retirement becomes possible

% Population growth of 1.2% IIJ2003, pg 762
Params.n=0.012;
% Social Security benefits
% IIJ2003 say they follow Huggett & Venture (1999) in defining this.
Params.b=0.1242; % HV1999 set b to 0.1242*output. IIJ2003 normalize to get output equal to one
Params.emax=2.47; % The maximum addition to AIME is set to 2.47 times average earnings
% Just an initial guess. Because labor is endogenous emax is determined in equilibrium.
Params.ebendpt1=0.2; % HV1999 set ebendpt1 to 0.2 times average earnings
% Just an initial guess. Because labor is endogenous ebend1 is determined in equilibrium.
Params.ebendpt2=1.24; % HV1999 set ebendpt2 to 1.24 times average earnings
% Just an initial guess. Because labor is endogenous ebend2 is determined in equilibrium.
Params.brate1=0.9; % The replacement rates for Social Security, see HV1999 pg 376
Params.brate2=0.32;
Params.brate3=0.15;

% Unemployment replacement rate
Params.phi=0.25; % IIJ2003, pg 763

% Taxes
tau_c=0.055; % Consumption tax (IIJ2003, pg 764)
tau_k=0.4; % Tax rate on interest earnings (IIJ2003, pg 764; IIJ2003 call this tau_a)
tau_l=0.2; % Labor income tax
tau_s=0.1; % Social security payroll tax
tau_u=0.05; % Unemployment-benefits tax
% Just an initial guess, tau_u is determined to fund the unemployment benefits


% Government spending
Params.G=0.18; % 0.18 of output, which is normalized to 1 in equilibrium

%Probability of surviving, conditional on being alive (IIJ2003 mention getting them from Faber, 1982; but this document does not exist online)
% Following are the orginal survival probabilites I was emailed by Selahattin Imrohoroglu in file 'psurv.dat'
Params.sj=[1.0000000; 0.99851000; 0.99844000; 0.99838000; 0.99832000; 0.99826000; 0.99820000; 0.99816000; 0.99815000; 0.99819000; 0.99826000; 0.99834000; 0.99840000; 0.99843000; 0.99841000; 0.99835000; 0.99828000; 0.99818000; 0.99807000; 0.99794000; 0.99778000; 0.99759000; 0.99737000; 0.99712000; 0.99684000; 0.99653000; 0.99619000; 0.99580000; 0.99535000; 0.99481000; 0.99419000; 0.99350000; 0.99278000; 0.99209000; 0.99148000; 0.99088000; 0.99021000; 0.98942000; 0.98851000; 0.98746000; 0.98625000; 0.98495000; 0.98350000; 0-98178000; 0.97974000; 0.97743000; 0.97489000; 0.97226000; 0.96965000; 0.96715000; 0.96466000; 0.96200000; 0.95907000; 0.95590000; 0.95246000; 0.94872000; 0.94460000; 0.94017000; 0.93555000; 0.93077000; 0.92570000; 0.92030000; 0.91431000; 0.90742000; 0.89948000];
Params.sj(end)=0;

% NEED TO CORRECT EPSILON_j (IT NEEDS TO BE NON ZERO FOR 'RETIREMENT' AGES
% Deterministic Life-cycle profile of earnings are interpolated from Hansen (1983).
% The following line is the contents of 'comboeff.dat' file I received from Selahattin Imrohoroglu (see a few lines above)
Params.epsilon_j=[0.36928407;  0.42120465;  0.49398951; 0.56677436; 0.63955922; 0.71234407; 0.78512892; 0.81551713; 0.84590532; 0.87629354; 0.90668176; 0.93706995; 0.96379826; 0.99052656; 1.0172549; 1.0439831; 1.0707114; 1.0819067; 1.0931018; 1.1042970; 1.1154923; 1.1266874; 1.1435058; 1.1603241; 1.1771424; 1.1939607; 1.2107789; 1.2015913; 1.1924037; 1.1832160; 1.1740283; 1.1648407; 1.1609988; 1.1571568; 1.1533149; 1.1494730; 1.1456310; 1.1202085; 1.0947859; 1.0693634; 1.0439408; 1.0185183; 0.99309576; 0.96767321];
Params.epsilon_j=[Params.epsilon_j; zeros(Params.J-Params.Jr+1,1)]; % Fill in the zeros for retirement age
Params.epsilon_j=epsilon_j/mean(epsilon_j); % IIJ2003 pg 763, "normalized to attain an average of unity over j=1,2,...,65".


%% Grids and exogenous shocks
l_grid=linspace(0,1,n_l)';
k_grid=20*linspace(0,1,n_k)'; % Phiaprime (see below) requires this to be a linear grid, and assumes the minimum grid value is zero
ebar_grid=linspace(0,2,n_ebar)'; % Phiaprime (see below) requires this to be a linear grid, and assumes the minimum grid value is zero
b_grid=[0;1]; % Indicates 'retired'

z_grid=[0; 1]; % Indicates 'employed'
pi_z=[0.06, 0.94; 0.06, 0.94]; % Note that IIJ2003 do not report this, you can find it in the 1999 working paper on page 16)
% This choice implies that the expected duration of unemployment is 1/(1-0.06)=1.0638 periods.
% Notice this actually makes the exogenous shock iid (rows of the markov transition matrix are identical)
statdist_z=pi_z(1,:)'; % follows from iid

% Following are needed for phiaprime (Case2 mapping from decision variables to endogenous states)
Params.kspacing=abs(k_grid(2)-k_grid(1));
Params.ebarspacing=abs(ebar_grid(2)-ebar_grid(1));

%% Get into format for VFI Toolkit
d_grid=[l_grid; k_grid; b_grid]; % Note that choosing l is implicitly determining ebar
a_grid=[k_grid; ebar_grid; b_grid];
% z_grid
% pi_z

n_d=[n_l,n_k,n_b];
n_a=[n_k,n_ebar,n_b];
% n_z

%% Calculate the population distribution across the ages (population at time
% t is made stationary by dividing it by \Pi_{i=1}^{t} (1+n_{i}) (product of all
% growth rates up till current period))
% I998 uses mewj, but has different notation for n (calls it rho) and sj (calls it psi_j). The VFI Toolkit needs to give
% it a name so that it can be automatically used when calculating model outputs.
Params.mewj=ones(Params.J,1);
for jj=2:Params.J
    Params.mewj(jj)=Params.mewj(jj-1)*(1/(1+Params.n))*Params.sj(jj-1);
end
Params.mewj=Params.mewj/sum(Params.mewj); % normalize to measure one
Params.mewj=Params.mewj'; % Age weights must be a row vector.

AgeWeightsParamNames={'mewj'}; % Many finite horizon models apply different weights to different 'ages'; eg., due to survival rates or population growth rates.

%% General equilibrium, and initial parameter values for the general eqm parameters

GEPriceParamNames={'r','tau_u','tau_s','Tr_beq','ebendpt1','ebendpt2','A','emax','beta'};
Params.r=0.06; % interest rate on assets
Params.tau_u=0.019; % set to balance unemployment benefits expenditure
Params.tau_s=0.122; % set to balance social security benefits expenditure

Params.Tr_beq=0.2; % Accidental bequests (IIJ2003 calls this kappa)

%% Now, create the return function 
% For quasi-hyperbolic preferences the last element of
% DiscountFactorParamNames must be the 'today-to-tomorrow' discount factor.
DiscountFactorParamNames={'beta','sj','gdiscount','beta0'}; % gdiscount is 1 in the baseline, it is needed for the extension to include deterministic productivity growth

ReturnFn=@(l,kprime,bprime,k,ebar,b,z,r,tau_c, tau_k, tau_l,tau_u, tau_s,gamma, epsilon_j,alpha,delta, A,Tr_beq,g,agej,Jr,ebarspacing,ebendpt1,ebendpt2,brate1,brate2,brate3) ImrohorogluImrohorogluJoines2003_ReturnFn(l,kprime,bprime,k,ebar,b,z,r,tau_c, tau_k, tau_l,tau_u, tau_s,gamma, epsilon_j,alpha,delta, A,Tr_beq,g,agej,Jr,ebarspacing,ebendpt1,ebendpt2,brate1,brate2,brate3)
ReturnFnParamNames={'r','tau_c', 'tau_k', 'tau_l','tau_u', 'tau_s','gamma', 'epsilon_j','alpha','delta', 'A','Tr_beq','g','agej','Jr','ebarspacing','ebendpt1','ebendpt2','brate1','brate2','brate3'}; %It is important that these are in same order as they appear in 'ImrohorogluImrohogluJoines2003_ReturnFn'

% Case 2 requires 'phiaprime' which determines next periods endogenous states from this periods decisions.
Case2_Type=12; % aprime=phi(d,a,z)
% For IIJ2003 phiaprime maps from the decision variables (l,kprime,bprime)
% to the next period exogenous states (kprime,ebarprime,bprime)
vfoptions.phiaprimematrix=2;
vfoptions.phiaprimedependsonage=1;
Params.n_a1=n_a(1); Params.n_a2=n_a(2); Params.n_a3=n_a(3);
PhiaprimeParamNames={'r','alpha','delta','A','epsilon_j','agej', 'Jr', 'emax','ebar','kgridspacing','ebargridspacing','n_a1','n_a2','n_a3'};
Phi_aprimeFn=@(l,kprime,bprime,k,ebar,b,z,r,alpha,delta,A,epsilon_j,agej, Jr, emax,ebar,kgridspacing,ebargridspacing,n_a1,n_a2,n_a3) ImrohorogluImrohorogluJoines2003_PhiaprimeFn(l,kprime,bprime,k,ebar,b,z,r,alpha,delta,A,epsilon_j,agej, Jr, emax,ebar,kgridspacing,ebargridspacing,n_a1,n_a2,n_a3)
% Note that Case2 is only required because the way ebar and l interact is
% different for working ages than it is for retirement.


%% Test: Now solve the value function iteration problem, just to check that things are working before we go to General Equilbrium

% First, the naive quasi-hyperbolic solution.
vfoptions.exoticpreferences='QuasiHyperbolic'
vfoptions.quasi_hyperbolic='Naive' % This is the default
tic;
[V_naive, Policy_naive]=ValueFnIter_Case2_FHorz(n_d,n_a,n_z,N_j, d_grid, a_grid, z_grid, pi_z, Phi_aprime, Case2_Type, ReturnFn, Params, DiscountFactorParamNames, ReturnFnParamNames);
toc

% Now, the sophisticated quasi-hyperbolic solution
vfoptions.exoticpreferences='QuasiHyperbolic'
vfoptions.quasi_hyperbolic='Sophisticated'
tic;
[V_sophisticated, Policy_sophisticated]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j, d_grid, a_grid, z_grid, pi_z, Phi_aprime, Case2_Type, ReturnFn, Params, DiscountFactorParamNames, ReturnFnParamNames);
toc

% Now, just the exponential discounting solution (without quasi-hyperbolic preferences)
vfoptions.exoticpreferences=0; % Turn of quasi-hyperbolic preferences (exponential discounting is the default)
DiscountFactorParamNames_Exp={'beta','sj','gdiscount'};
tic;
[V, Policy]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j, d_grid, a_grid, z_grid, pi_z, Phi_aprime, Case2_Type, ReturnFn, Params, DiscountFactorParamNames_Exp, ReturnFnParamNames);
toc
% Switch quasi-hyperbolic preferences back on
vfoptions.exoticpreferences='QuasiHyperbolic'


% max(max(max(max(Policy))))<n_a % Double check that never try to leave top of asset grid.
% sum(sum(sum(sum(Policy==n_a))))

%% Once we have the policy functions simulations, agent distributions, etc. are all computed just as normal.

%% Initial distribution of agents at birth (j=1)
jequaloneDist=zeros(n_k,n_ebar,n_b,n_z); 
jequaloneDist(1,1,1,:)=statdist_z; % All agents born with zero assets, zero ebar, and not retired and with z based on stationary distribution of the exogenous process on labour productivity units (I1998, pg 326)

%% Test: StationaryDist commands
StationaryDist_naive=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightsParamNames,Policy_naive,n_d,n_a,n_z,N_j,pi_z,Params);
StationaryDist_sophisticated=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightsParamNames,Policy_sophisticated,n_d,n_a,n_z,N_j,pi_z,Params);
StationaryDist=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightsParamNames,Policy,n_d,n_a,n_z,N_j,pi_z,Params);

%% Set up all the conditions for general equilibrium
% Steady State Aggregates (important that ordering of Names and Functions is the same)
FnsToEvaluateParamNames=struct();
FnsToEvaluateParamNames(1).Names={};
FnsToEvaluate_K = @(l_val,kprime_val,ebarprime_val,b_val,k_val,ebar_val,a_val,z_val) k_val; % Aggregate assets K
FnsToEvaluateParamNames(2).Names={'epsilon_j'};
FnsToEvaluate_2 = @(l_val,kprime_val,ebarprime_val,b_val,k_val,ebar_val,a_val,z_val,epsilon_j) epsilon_j*l*z_val; % Aggregate effective labour supply (in efficiency units)
FnsToEvaluateParamNames(3).Names={'sj','n'};
FnsToEvaluate_3 = @(l_val,kprime_val,ebarprime_val,b_val,k_val,ebar_val,a_val,z_val,sj,n) (1-sj)*kprime_val/(1+n); % Tr, accidental bequest transfers % The divided by (1+n) is due to population growth and because this is Tr_{t+1}
FnsToEvaluateParamNames(4).Names={'SSdivw','r','A','alpha','delta','I_j'}; % MODIFY TO ALLOW FOR g
% FnsToEvaluate_4 = @(aprime_val,a_val,z_val,SSdivw,r,A,alpha,delta,I_j) ((1-alpha)*(A^(1/(1-alpha)))*((r+delta)/alpha)^(alpha/(alpha-1)))*SSdivw*(1-I_j); % Total social security benefits
FnsToEvaluateParamNames(5).Names={'r','A','alpha','delta','h','I_j','epsilon_j'};
FnsToEvaluate_5 = @(aprime_val,a_val,z_val,r,A,alpha,delta,h,I_j,epsilon_j) ((1-alpha)*(A^(1/(1-alpha)))*((r+delta)/alpha)^(alpha/(alpha-1)))*I_j*h*epsilon_j*z_val; % Labour income (this is also the tax base for various taxes)
FnsToEvaluateParamNames(6).Names={'r','A','alpha','delta','h','I_j','epsilon_j'};
FnsToEvaluate_6 = @(aprime_val,a_val,z_val,r,A,alpha,delta,h,I_j,epsilon_j) ((1-alpha)*(A^(1/(1-alpha)))*((r+delta)/alpha)^(alpha/(alpha-1)))*I_j*epsilon_j*(1-z_val); % Unemployment benefits
FnsToEvaluate={FnsToEvaluate_K,FnsToEvaluate_2,FnsToEvaluate_3,FnsToEvaluate_4,FnsToEvaluate_5,FnsToEvaluate_6};

% General Equilibrium Equations
% Recall that GEPriceParamNames={'r','tau_u','tau_s','Tr_beq','ebendpt1','ebendpt2','A','emax','beta'};
% In following lines GEprices is the vector of these and so, e.g., GEprices(2) is tau_u.
GeneralEqmEqnParamNames=struct();
GeneralEqmEqnParamNames(1).Names={'A','alpha','delta'};
GeneralEqmEqn_1 = @(AggVars,GEprices,A,alpha,delta) GEprices(1)-(A*(alpha)*(AggVars(1)^(alpha-1))*(AggVars(2)^(1-alpha))-delta); % Rate of return on assets is related to Marginal Product of Capital
GeneralEqmEqnParamNames(2).Names={'zeta'};
GeneralEqmEqn_2 = @(AggVars,GEprices,zeta) GEprices(2)*AggVars(5)-zeta*AggVars(6); % Unemployment benefits are budget balances
GeneralEqmEqnParamNames(3).Names={};
GeneralEqmEqn_3 = @(AggVars,GEprices) GEprices(3)*AggVars(5)-AggVars(4); % Social security is budget balances
GeneralEqmEqnParamNames(4).Names={};
GeneralEqmEqn_4 = @(AggVars,GEprices) GEprices(4)-AggVars(3); % Accidental bequests (adjusted for population growth) are equal to transfers received (this is essentially eqn 14)
GeneralEqmEqnParamNames(5).Names={};
GeneralEqmEqn_5 = @(AggVars,GEprices) GEprices(5)-0.2*AggVars(5); % ebendpt1 is 0.2 of average earnings
GeneralEqmEqnParamNames(6).Names={};
GeneralEqmEqn_6 = @(AggVars,GEprices) GEprices(6)-1.24*AggVars(5); % ebendpt2 is 1.24 of average earnings
GeneralEqmEqnParamNames(7).Names={};
GeneralEqmEqn_7 = @(AggVars,GEprices) 1-GEprice(7)*(AggVars(1)^alpha)*(AggVars(2)^(1-alpha)); % Choose A to normalize output to 1 (output is Cobb-Douglas prodn fn)
GeneralEqmEqnParamNames(8).Names={};
GeneralEqmEqn_8 = @(AggVars,GEprices) GEprices(8)-2.47*AggVars(5); % emax is 2.47 times average earnings
GeneralEqmEqnParamNames(9).Names={};
GeneralEqmEqn_9 = @(AggVars,GEprices) 2.52-AggVars(1); % Choose beta to get capital-output ratio of 2.52 (output is already normalized to 1, so just need to target capital stock of 2.52)
GeneralEqmEqns={GeneralEqmEqn_1, GeneralEqmEqn_2, GeneralEqmEqn_3, GeneralEqmEqn_4, GeneralEqmEqn_5, GeneralEqmEqn_6, GeneralEqmEqn_7, GeneralEqmEqn_8, GeneralEqmEqn_9};


%% Computing General Eqm is just as normal (then create life-cycle profile of assets)
FnsToEvaluateParamNames2=struct();
FnsToEvaluateParamNames2(1).Names={};
FnsToEvaluate2={FnsToEvaluate_K}

% Quasi-hyperbolic, naive
vfoptions.exoticpreferences='QuasiHyperbolic'
vfoptions.quasi_hyperbolic='Naive' % This is the default
[p_eqm_naive,p_eqm_index, GeneralEqmEqnsValues]=HeteroAgentStationaryEqm_Case2_FHorz(jequaloneDist,AgeWeightsParamNames,n_d, n_a, n_z, N_j, 0, pi_z, d_grid, a_grid, z_grid, Phi_aprimeKron, Case2_Type, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, ReturnFnParamNames, PhiaprimeParamNames, FnsToEvaluateParamNames, GeneralEqmEqnParamNames, GEPriceParamNames,[], [], vfoptions);
Params.r=p_eqm_naive.r;
Params.tau_u=p_eqm_naive.tau_u;
Params.tau_s=p_eqm_naive.tau_s;
Params.Tr_beq=p_eqm_naive.Tr_beq;
Params.ebendpt1=p_eqm_naive.ebendpt1;
Params.ebendpt2=p_eqm_naive.ebendpt2;
Params.A=p_eqm_naive.A;
Params.emax=p_eqm_naive.emax;
Params.beta=p_eqm_naive.beta;

[V_naive, Policy_naive]=ValueFnIter_Case2_FHorz(n_d,n_a,n_z,N_j, d_grid, a_grid, z_grid, pi_z, Phi_aprime, Case2_Type, ReturnFn, Params, DiscountFactorParamNames, ReturnFnParamNames);
StationaryDist_naive=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightsParamNames,Policy_naive,n_d,n_a,n_z,N_j,pi_z,Params);
AgeConditionalStats_naive=LifeCycleProfiles_FHorz_Case2(StationaryDist_naive,Policy_naive,FnsToEvaluate2,FnsToEvaluateParamNames2,Params,0,n_a,n_z,N_j,0,a_grid,z_grid);

% Quasi-hyperbolic, sophisticated
vfoptions.exoticpreferences='QuasiHyperbolic'
vfoptions.quasi_hyperbolic='Sophisticated' % This is the default
[p_eqm_sophisticated,p_eqm_index, GeneralEqmEqnsValues]=HeteroAgentStationaryEqm_Case2_FHorz(jequaloneDist,AgeWeightsParamNames,n_d, n_a, n_z, N_j, 0, pi_z, d_grid, a_grid, z_grid, Phi_aprimeKron, Case2_Type, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, ReturnFnParamNames, PhiaprimeParamNames, FnsToEvaluateParamNames, GeneralEqmEqnParamNames, GEPriceParamNames,[], [], vfoptions);
Params.r=p_eqm_sophisticated.r;
Params.tau_u=p_eqm_sophisticated.tau_u;
Params.tau_s=p_eqm_sophisticated.tau_s;
Params.Tr_beq=p_eqm_sophisticated.Tr_beq;
Params.ebendpt1=p_eqm_sophisticated.ebendpt1;
Params.ebendpt2=p_eqm_sophisticated.ebendpt2;
Params.A=p_eqm_sophisticated.A;
Params.emax=p_eqm_sophisticated.emax;
Params.beta=p_eqm_sophisticated.beta;

[V_sophisticated, Policy_sophisticated]=ValueFnIter_Case2_FHorz(n_d,n_a,n_z,N_j, d_grid, a_grid, z_grid, pi_z, Phi_aprime, Case2_Type, ReturnFn, Params, DiscountFactorParamNames, ReturnFnParamNames);
StationaryDist_sophisticated=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightsParamNames,Policy_sophisticated,n_d,n_a,n_z,N_j,pi_z,Params);
AgeConditionalStats_sophisticated=LifeCycleProfiles_FHorz_Case2(StationaryDist_sophisticated,Policy_sophisticated,FnsToEvaluate2,FnsToEvaluateParamNames2,Params,0,n_a,n_z,N_j,0,a_grid,z_grid);

% Standard Exponential discounting
vfoptions.exoticpreferences=0
[p_eqm_expdisc,p_eqm_index, GeneralEqmEqnsValues]=HeteroAgentStationaryEqm_Case2_FHorz(jequaloneDist,AgeWeightsParamNames,n_d, n_a, n_z, N_j, 0, pi_z, d_grid, a_grid, z_grid, Phi_aprimeKron, Case2_Type, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, ReturnFnParamNames, PhiaprimeParamNames, FnsToEvaluateParamNames, GeneralEqmEqnParamNames, GEPriceParamNames,[], [], vfoptions);
Params.r=p_eqm_expdisc.r;
Params.tau_u=p_eqm_expdisc.tau_u;
Params.tau_s=p_eqm_expdisc.tau_s;
Params.Tr_beq=p_eqm_expdisc.Tr_beq;
Params.ebendpt1=p_eqm_expdisc.ebendpt1;
Params.ebendpt2=p_eqm_expdisc.ebendpt2;
Params.A=p_eqm_expdisc.A;
Params.emax=p_eqm_expdisc.emax;
Params.beta=p_eqm_expdisc.beta;

[V_expdisc, Policy_expdisc]=ValueFnIter_Case2_FHorz(n_d,n_a,n_z,N_j, d_grid, a_grid, z_grid, pi_z, Phi_aprime, Case2_Type, ReturnFn, Params, DiscountFactorParamNames, ReturnFnParamNames);
StationaryDist_expdisc=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightsParamNames,Policy_expdisc,n_d,n_a,n_z,N_j,pi_z,Params);
AgeConditionalStats_expdisc=LifeCycleProfiles_FHorz_Case2(StationaryDist_expdisc,Policy_expdisc,FnsToEvaluate2,FnsToEvaluateParamNames2,Params,0,n_a,n_z,N_j,0,a_grid,z_grid);

%% Plot the asset life-cycle profiles
figure(1)
plot((1:1:Params.J)+20,AgeConditionalStats_expdisc.Mean)
hold on
plot((1:1:Params.J)+20,AgeConditionalStats_naive.Mean)
plot((1:1:Params.J)+20,AgeConditionalStats_sophisticated.Mean)
hold off
legend('Exponential Discounting','Naive Quasi-Hyperbolic', 'Sophisticated Quasi-Hyperbolic')
xlabel('Age')
ylabel('Assets')
title('Effects of Quasi-Hyperbolic Discounting on Life-Cycle Profile of Assets')
saveas(gcf,'./SavedOutput/Graphs/IIJ2003_LifeCycleProfileOfAssets.png')







