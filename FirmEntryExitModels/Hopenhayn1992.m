% Hopenhayn (1992) - Entry, Exit and Firm Dynamics in Long Run Equilibrium
% While this code implements the framework of Hopenhayn (1992) the model/example itself
% is taken from slides of Chris Edmond: http://www.chrisedmond.net/phd2019/econ90003_lecture19.pdf
% The numerical example begins on slide 27.
% This involves solving a stationary competitive equilibrium.

% Hopyenhayn (1992) does contain an 'Example' on pg 1144, but it is kind of silly.
% There is a distribution of entrants at two productivity levels, 0 and 1, with mass 0.9 of 
% former and mass 0.1 of later (v({0})=0.9 and v({1})=0.1 in notation of the paper). 
% But both of these are below the cutoff level of
% productivity for exit (x*=6.52). So the equilibrium distribution of firms
% is actually undefined (without giving some more info) and the role of entry and exit in shaping it is
% that every period some firms choose to enter, appear at the productivity levels 0 & 1,
% produce (due to timing of entry and exit decisions being beginining and end of period respectively),
% and then exit! [eqm distribution of firms is undefined, it can involve any
% pre-existing firms as long as their productivity is above x*=6.52. For the statistics
% reported in Table 2 the actual distribution of firms is irrelevant. [These statistics would be the 
% same for the continuua of equilibria associated with different stationary distributions, 
% D(p) and w(N) functions; Hopenhayn & Rogerson (1993) provide some discussion of this issue of a continuua of 
% equilibria in a closely related model, although in their case it arises from no entry nor exit in equilibrium,
% rather than entry and immediate exit as here.]

vfoptions.parallel=0
simoptions.parallel=0

%%
n_z=101; % I here call z, what Hopenhayn calls 'varphi'; namely the firm-specific (permanent) productivity level.

Params.beta=0.8;
Params.alpha=2/3;
Params.cf=20;

Params.ce=40;

Params.Dbar=100; % Part of the (inverse) demand function that determines the price.

Params.p=4.3; % Initial guess, this will be determined by general equilibrium condition.
Params.w=1; % Normalization

Params.rho=0.9;
Params.sigma_epsilon=0.2;
Params.logzbar=1.4;

tauchenoptions.parallel=vfoptions.parallel;
Params.q=4;
[z_grid, pi_z]=TauchenMethod((1-Params.rho)*Params.logzbar,Params.sigma_epsilon^2,Params.rho,n_z,Params.q,tauchenoptions); %[states, transmatrix]=TauchenMethod_Param(mew,sigmasq,rho,znum,q,Parallel,Verbose), transmatix is (z,zprime)
% Compute the stationary distribution of this Markov (this could be done more easily/directly/analytically)
pistar_z=ones(size(z_grid))/n_z; % Initial guess
dist=1;
while dist>10^(-9)
    pistar_z_old=pistar_z;
    pistar_z=(pi_z')*pistar_z;
    dist=max(abs(pistar_z-pistar_z_old));
end

%%
n_a=1;
a_grid=1;
n_d=401;
d_grid=((logspace(0,pi,n_d)-1)/(pi-1))'*500; % logarithmically spaced 0 to 500

%% Value Function Problem (with endogenous exit)
% Entry is irrelevant to the value function problem. Exit is relevant, and
% treated differently based on whether it is endogenous or exogenous.

% For exogenous exit, you would simply need to include the 'conditional surivival probability' as another 'DiscountFactorParamNames'
DiscountFactorParamNames={'beta'};
% For endogenous exit, you need to set
vfoptions.endogenousexit=1;
% By default the toolkit assumes that exit decision affects current period
% returns, but in Hopyenhayn (1992) exit decision is made at end-of-period.
% Because of there is no endogenous state in Hopenyhayn (1992) model this
% this is actually irrelevant to the decisions, and all we need is to activate an option that
% 'keeps' the Policy choice of the current period that would have been be made if the
% firm did not (already choose to) exit.
vfoptions.keeppolicyonexit=1;
% Outside of models like Hopenhayn (1992), which have no endogenous state, you are unlikely to use this option.
% We also need to create 'vfoptions.ReturnToExitFn' (and 'vfoptions.ReturnToExitFnParamNames'), as below.

ReturnFn=@(n_val,aprime_val, a_val, s_val, p, alpha, cf) Hopenhayn1992_ReturnFn(n_val,aprime_val, a_val, s_val, p, alpha, cf);
ReturnFnParamNames={'p', 'alpha', 'cf'}; %It is important that these are in same order as they appear in 'Hopenhayn1992_ReturnFn'

% For endogenous exit, also need to define the 'return to exit'.
vfoptions.ReturnToExitFn=@(a_val, s_val) 0;
% Remark: if a more complex 'return to exit' was desired it could be
% created in much the same way as the ReturnFn is created, but depends only
% on (a,z) variables (and any parameters passed using ReturnToExitFnParamNames).
% The following commented out line provides an 'example'
% ReturnToExitFn=@(a_val, s_val) Hopenhayn1992_ReturnToExitFn(a_val, s_val);
vfoptions.ReturnToExitFnParamNames={}; %It is important that these are in same order as they appear in 'Hopenhayn1992_ReturnToExitFn'

% Check that everything is working so far by solving the value function
if vfoptions.parallel==2
    V0=zeros(n_a,n_z,'gpuArray');
else
    V0=zeros(n_a,n_z);
end
vfoptions % print them to screen
[V,Policy,ExitPolicy]=ValueFnIter_Case1(V0, n_d,n_a,n_z,d_grid,a_grid,z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, ReturnFnParamNames, vfoptions);

%% Stationary Distribution of Agents with entry and exit
% Both entry and exit matter for stationary distribution of agents. 
% In principle whether they are endogenous or exogenous is irrelevant.
simoptions.agententryandexit=1;
% In principle whether exit or entry is exogenous is not needed info, but
% for algorithm turns out to be useful info
simoptions.endogenousexit=1;

% Edmond's slides do not spell out what is done in terms of new entrants. I
% will set their distribution equal to the stationary distribution of the
% markov process on productivity.
% simoptions.DistOfNewAgents=pistar_z;
EntryExitParamNames.DistOfNewAgents={'pistar_z'};
Params.pistar_z=pistar_z;
% Note: VFI Toolkit requires the DistOfNewAgents to be a pdf (so unit mass), and then uses
% the 'MassOfNewAgents' to understand how many there will be entering
% relative to existing agents. (MassOfExistingAgents is kept track of.)
Params.Ne=1;
EntryExitParamNames.MassOfNewAgents={'Ne'};
% EntryExitParamNames.MassOfExistingAgents={'N'};
% simoptions.MassOfNewAgents=1;
% Note: this is just an initial guess, as it will anyway need to be
% determined in general equilibrium. (This can be done analytically, but we won't.)

% Exit can be given in one of two forms. As a function, or as a matrix.
% With endogenous entry, it is always best to give as following matrix.
% simoptions.CondlProbOfSurvival=1-ExitPolicy;
EntryExitParamNames.CondlProbOfSurvival={'zeta'};
Params.zeta=1-ExitPolicy;
% For exogenous entry, things like the following might be easier to use,
% simoptions.CondlProbOfSurvivalParamNames={'lambda'};
% simoptions.CondlProbOfSurvival=@(lambda) 1-lambda;
% This conditional probability can be a state-dependent parameter, in which case this input would be a vector or matrix, etc.

% You can optionally give a,
% simoptions.MassOfExistingAgents

simoptions % Print them to screen

% Check that everything is working so far by solving the simulation of
% agent distribution to get the stationary distribution.
% [StationaryDist, MassOfExistingAgents]=StationaryDist_Case1(Policy,n_d,n_a,n_z,pi_z, simoptions,Params,EntryExitParamNames);
StationaryDist=StationaryDist_Case1(Policy,n_d,n_a,n_z,pi_z, simoptions,Params,EntryExitParamNames);
% Params.N=MassOfExistingAgents;

% Note: When using models, such as entry and exit, where the mass of agents is not equal to 1
% the toolkit will automatically keep track of distributions as
% StationaryDist.pdf and StationaryDist.Mass

% Note: it is possible to do just entry or just exit.
% To do just entry, simply set CondlProbOfExit=0;
% To do just exit, simply set DistOfNewAgents=0;
% (It is not envisaged that you are likely to want either of entry or exit without the other)

%% General Equilibrium with endogenous entry and endogenous exit.
% The entry decision imposes a general equilibrium condition. In practice it determines the price level p.
% Notice that this depends on the value function of existing firms and on the 'DistOfNewAgents'.
% If entry is exogenous it would not impose any general equilibrium condition.
% Notice that exit (whether endogenous or exogenous) does not involve any
% general equilibrium conditions.
% Since the stationary equilibrium commands already need to know all the
% vfoptions and simoptions, we don't need to add any further info to highlight that 
% this problem includes entry and exit as these
% already contain everything the stationary equilibrium command needs to know.

%Use the toolkit to find the equilibrium prices
GEPriceParamNames={'p','Ne'};

FnsToEvaluateParamNames(1).Names={'alpha'};
% Note: With entry-exit the mass of the distribution of agents often
% matters. So it becomes an extra input arguement in all functions to be
% evaluated.
FnsToEvaluateFn_1 = @(d_val, aprime_val,a_val,z_val,agentmass,alpha) z_val*(d_val^alpha); % Total output
FnsToEvaluate={FnsToEvaluateFn_1};

% Just to test: (note, is same command as usual, the 'StationaryDist' contains the
% mass, and this is automatically used internally to make all relevant
% changes to the algorithms)
simoptions.keeppolicyonexit=1; % Is not needed for simulation, but is needed for EvalFnOnAgentDist.
AggVars=EvalFnOnAgentDist_AggVars_Case1(StationaryDist, Policy, FnsToEvaluate, Params, FnsToEvaluateParamNames, n_d, n_a, n_z, d_grid, a_grid, z_grid, simoptions.parallel, simoptions, EntryExitParamNames);

% The general equilibrium condition is that the EV^e-ce=0.
% This does not fit standard format for general equilibrium conditions.
heteroagentoptions.specialgeneqmcondn={0,'entry'};
% Certain kinds of general equilibrium conditions that are non-standard can
% be used via heteroagentoptions.specialgeneqmcondn
GeneralEqmEqnParamNames(1).Names={'Dbar'};
GeneralEqmEqn_1 = @(AggVars,p,Dbar) AggVars/Dbar-p(1); %The requirement that the price is determined by the demand eqn
GeneralEqmEqnParamNames(2).Names={'beta','ce'};
GeneralEqmEqn_Entry = @(EValueFn,p,beta,ce) beta*EValueFn-ce; % Free entry conditions (expected returns equal zero in eqm).
GeneralEqmEqns={GeneralEqmEqn_1, GeneralEqmEqn_Entry};
% Note that GeneralEqmEqn_Entry need to be pointed out as special because
% it depends on the distribution of entrants and not the distribution of
% existing agents (all standard general eqm conditions involve the later).

% Worth noting that in the basic case of a stationary general eqm in model
% of Hopenhayn (1992) these two equations can be simplified using some
% basic analytical tricks (see Chris Edmond's slides). We don't do these
% here as want to be able to then easily do transition path, and the
% simplifications would not work for the transition.

heteroagentoptions.verbose=1;
n_p=0;
disp('Calculating price vector corresponding to the stationary eqm')
% tic;
% NOTE: EntryExitParamNames has to be passed as an additional input compared to the standard case.
[p_eqm_initial,p_eqm_index_initial, GeneralEqmCondition_initial]=HeteroAgentStationaryEqm_Case1(V0, n_d, n_a, n_z, n_p, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, ReturnFnParamNames, FnsToEvaluateParamNames, GeneralEqmEqnParamNames, GEPriceParamNames,heteroagentoptions, simoptions, vfoptions, EntryExitParamNames);
% findeqmtime=toc
Params.p=p_eqm_initial.p;
Params.Ne=p_eqm_initial.Ne;

%% In lecture notes Chris Edmonds goes on to plot two graphs.

% First contains the value function and some things like exit cutoff level.
[~,z_starindex]=min(ExitPolicy); % By default Matlab reports index of the first min() when there are more than one element with that value.
z_star=exp(z_grid(z_starindex));
EV=pi_z*V';
figure(1)
plot(exp(z_grid),V, exp(z_grid),EV)
title('Figure from Slide 28 of Chris Edmonds (except these are after exit decision)')
xlim([0,8])
legend('V','E[Vprime|z]')

% Second relates to the distribution of firms, and a weighted distribution by employment.
% To graph this we need to calculate employment, which is by coincidence the only policy.
EmploymentDecisionsIndexes=shiftdim(Policy(1,1,:),1); % Policy also contains the trivial decision about the 'a' endogenous state.
EmploymentDecisions=d_grid(EmploymentDecisionsIndexes)';
DistOfEmployment=(EmploymentDecisions.*StationaryDist.pdf)/sum(EmploymentDecisions.*StationaryDist.pdf);
DistOfEmployment=DistOfEmployment/sum(DistOfEmployment);

figure(2)
subplot(2,1,1); plot(exp(z_grid),StationaryDist.pdf, exp(z_grid),pistar_z, exp(z_grid),DistOfEmployment')
title('Figure from Slide 29 of Chris Edmonds (pdfs; differ due to grid)')
legend('share firms (mu)', 'fbar', 'share of employment')
subplot(2,1,2); plot(exp(z_grid),cumsum(StationaryDist.pdf), exp(z_grid),cumsum(pistar_z), exp(z_grid),cumsum(DistOfEmployment'))
title('Figure from Slide 29 of Chris Edmonds (cdfs)')
% legend('share firms (mu)', 'fbar', 'share of employment')

% I recommend using the cdf. Simply that the pdf is (definitionally) sensitive to the 
% number of grid points used to approximate z. 
% [Alternatively you can plot a nonparametric probability distribution function estimate, 
% instead of just the actual probability density function. If you wish to do this I recommend
% kernel-smoothed non-parametric estimators of the probability distribution funciton. 
% Note distribution vs density.]




