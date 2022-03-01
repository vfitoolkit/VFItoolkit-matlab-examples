% Example based on Guerrieri & Lorenzoni (2017) - Credit Crises, Precautionary Savings, and the Liquidity Trap
% It computes the initial equilibrium, the final equilibrium, the flexible
% price (baseline) transition, and the sticky-wage (New-Keynesian) transition.
% It reproduces most of Figures 1 to 8 of the original paper.
% For details of the model itself, see the original paper: https://ideas.repec.org/a/oup/qjecon/v132y2017i3p1427-1467..html
%
% Note: The values of the parameters reported in Table 1 GL2017 paper are not the actual parameters. 
% The values in this example are instead taken from their codes. See the replication for full explanation.
%
% Relatedly, the GL2017 code imports an 'income_process.mat'. To avoid this
% dependency I here use Tauchen method with hyperparameter q=2.1 which
% gives almost the exact same grid and transition for the exogenous shock process on z.
%
% This example has all 'save' and 'load' commands commented out. You can uncomment them if you want.
%
% This example is hard-coded to depend on having a GPU.

%% To translate Guerrieri & Lorenzoni (2017) into the standard setup of VFI Toolkit I use following:
% d variables: n_it
% aprime variables: b_it+1
% a variables: b_it
% z variables: theta_it, e_it

%% Number of grid points
n_d=41; % Endogenous labor 
n_a=2^10; % Assets
n_theta=13;

%% Parameters
% (mostly from G&L2017, Table 1)
Params.beta=0.9774; %Model period is one-quarter of a year.

Params.gamma=4; % Coefficient of relative risk aversion
Params.eta=1.5; % Curvature of utility from leisure
Params.psi=15.8819; % Coefficient on leisure in utility

Params.pi_eu=0.0573; % Transition to unemployment (the last three comes from their code, not paper)
Params.pi_ue=0.882; % Transition to employment
% Note: you cannot change pi_eu or pi_ue without changing the return
% function because the implied uncondtional probability of being unemployed
% is hard coded there.
Params.rho=0.967; % Persistence of productivity shock (G&L2017 call this rho)
Params.sigmasq_epsilon=0.017; % Variance of the shock to the log-AR(1) process on labour productivity.
Params.tauchenq=2.1; % I have reverse engineered this value from the grid of GL2017 (copy of their grid is included in their codes). They themselves reverse engineered choice of roughly 2.1 so that the variance of the resulting process (variance of 12-state markov logz) is as close as possible to what it should be (variance of true AR(1) logz).
% This choice of tauchenq is crucial to the results/replication of GL2017. It means
% that the min and max productivity shocks in the model, while they have the correct variance, have a range which
% is roughly just +-standard deviation. Because this range of shocks is so (empirically unrealistically) small
% the equilibrium interest rate is higher than it would otherwise be; if you use tauchenq=3 the zero-lower bound 
% on interest rates in the 'new-keynesian sticky-wage' transition lasts for decades.
% (See 'Comment' on blog at vfitoolkit.com: http://www.vfitoolkit.com/updates-blog/2020/transition-paths-example-based-on-guerrieri-lorenzoni-2017/ )

Params.v=0.1670; % Unemployment benefit
Params.B=2.6712; % Bond supply
Params.Bprime=Params.B; % Bond supply is unchanging (for most of the paper); this is needed as it is part of government budget constraint that determines lump-sum tax tau_t. (Obviously B=Bprime must hold in any stationary general eqm.)
% For creating graphs later it will be useful to store both the initial and final values of phi
Params.phi_initial=1.6005;
Params.phi_final=0.8767;
Params.phi=Params.phi_initial; % Borrowing limit

Params.r=0.006; % General equilibrium interest rate (I use this instead of q). This is just an initial guess. [Model is quarterly, so this corresponds to roughly 2.4% (=4*r) at an annual rate. (Ignoring compounding.)]
% Params.tau % Lump-sum taxes, determined in general equilibrium (I implement it directly inside the ReturnFn)

Params.omega=0; % This is not actually needed for anything until we get to the 'Sticky Wage' model (Section 4 of GL2017, pg 1450)
% In the New Keynesian Sticky Wages model this appears as a wedge that
% increases the utility of leisure (intended to capture a fall in labor
% demand as a result of wages not falling). See GL2017 for explanation.

%% Create the grid for exogenous shocks z, and the transition matrix for these (using Tauchen method)

% Create markov process for the exogenous income (based on idea of employment and unemployment states, following Imrohoroglu, 1989).
[theta1_grid,pi_theta1]=discretizeAR1_Tauchen(0,Params.rho,sqrt(Params.sigmasq_epsilon),n_theta-1,Params.tauchenq);
z_grid=[0; exp(theta1_grid)];
% G&L2017, pg 1438 "when first employed, workers draw theta from its unconditional distribution"; so here compute the unconditional distribution
pistar_theta1=ones(n_theta-1,1)/(n_theta-1);
for ii=1:10^4 % There is a more efficient form to do this directly from a formula but I am feeling lazy. %FIX THIS LATER!!!
    pistar_theta1=pi_theta1'*pistar_theta1; 
end

pi_z=[(1-Params.pi_ue), Params.pi_ue*pistar_theta1'; Params.pi_eu*ones(n_theta-1,1),(1-Params.pi_eu)*pi_theta1];
% Rows did not sum to one due to rounding errors at order of 10^(-11), fix this
pi_z=pi_z./sum(pi_z,2);
pistar_z=ones(n_theta,1)/n_theta;
for ii=1:10^4 %  % There is a more efficient way to do this directly from a formula but I am feeling lazy. %FIX THIS LATER!!!
    pistar_z=pi_z'*pistar_z; % Formula could be used to find stationary dist of the employment unemployment process, then just combine with stationary dist of theta1, which is already calculated
end
% "The average level of theta is chosen so that yearly output in the initial steady state is normalized to 1"
z_grid=z_grid/sum(z_grid.*pistar_z);
% Double-check that this is 1
% sum(z_grid.*pistar_z)

% That the "normalized to 1" refers to E[theta] and not E[n*theta] is clear from setting
% v=0.1 to satisfy "For the unemployment benefit, we also follow Shimer
% (2005) and set it to 40% of average labor income." (pg 1438)
% Note, they do not actually ever normalize to 1 in the codes, GL2017 has E[theta]=1.0718


%% Grids
% Set grid for asset holdings
Params.alowerbar=-1.25*Params.phi; % This seems reasonable (No-one can go below -Params.phi in any case). Note that Fischer deflation experiment (second last part of paper) won't work unless this is below -1.1*phi
Params.aupperbar=20; % Not clear exactly what value is appropriate, have gone with this, and checked that increasing it makes no difference.
a_grid=(Params.aupperbar-Params.alowerbar)*(1/(exp(1)-1))*(exp(linspace(0,1,n_a)')-1)+Params.alowerbar;
% GL2017 codes use alowerbar of -2 and aupperbar of 50, but since no-one
% goes anywhere near 50 this seems excessive (even for them 0.9995 of population ends up below 20.8690). 

%Bring model into the notational conventions used by the toolkit
d_grid=linspace(0,1,n_d)'; % Labor supply

n_z=n_theta;

%Create descriptions of SS values as functions of d_grid, a_grid, s_grid &
%pi_s (used to calculate the integral across the SS dist fn of whatever
%functions you define here)
FnsToEvaluate.A = @(d, aprime,a,z) a; % Aggregate assets (which is this periods state)

%Now define the functions for the General Equilibrium conditions
    %Should be written as LHS of general eqm eqn minus RHS, so that 
    %the closer the value given by the function is to zero, the closer 
    %the general eqm condition is to holding.
GeneralEqmEqns.ZeroAggAssets = @(A,B) A-B; %The requirement that the aggregate assets (lending and borrowing; by government and private) equal zero

%% 
DiscountFactorParamNames={'beta'};

ReturnFn=@(d,aprime,a,z, r, gamma, psi, eta, phi, v,B,Bprime,omega) GuerrieriLorenzoni2017_ReturnFn(d,aprime,a,z, r, gamma, psi, eta, phi, v,B,Bprime,omega);

%% Solve the initial stationary equilibrium
% GL2007 refer to this as the initial steady-state equilibrium, which it is not. It is the inital stationary equilibrium. (there are plenty of shocks at the idiosyncratic level, hence not steady-state which means the absence of shocks)
% Comment: The 'steady-state equilibrium' is standard in literature, but I personally consider it misleading and 
% prefer the term 'stationary equilibrium' as being more accurate.
% [Relatedly: there is a difference in a Representative Agent DSGE between the steady-state value and the mean 
% value of the model with shocks! Think carefully about why this is and you will understand my insistence of use 
% of stationary rather than steady-state language.]

%Use the toolkit to find the equilibrium price index
GEPriceParamNames={'r'}; %,'tau'

heteroagentoptions.verbose=1;
disp('Calculating price vector corresponding to the stationary eqm')
% tic;
[p_eqm_initial,p_eqm_index_initial, ~]=HeteroAgentStationaryEqm_Case1(n_d, n_a, n_z, 0, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, [], [], [], GEPriceParamNames,heteroagentoptions);
% findeqmtime=toc
Params.r=p_eqm_initial.r;


%% Now that we know what the equilibrium price is, lets calculate a bunch of other things associated with the equilibrium

disp('Calculating various equilibrium objects')
[~,Policy_initial]=ValueFnIter_Case1(n_d,n_a,n_z,d_grid,a_grid,z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, []);

% PolicyValues=PolicyInd2Val_Case1(Policy,n_d,n_a,n_s,d_grid,a_grid);

StationaryDist_initial=StationaryDist_Case1(Policy_initial,n_d,n_a,n_z,pi_z);

AggVars_initial=EvalFnOnAgentDist_AggVars_Case1(StationaryDist_initial, Policy_initial, FnsToEvaluate,Params, [],n_d, n_a, n_z, d_grid, a_grid,z_grid);
    
% save ./SavedOutput/GuerrieriLorenzoni2017_initial.mat Params p_eqm_initial Policy_initial StationaryDist_initial AggVars_initial n_d n_a n_z

% Some things you might want to take a look at just to see what is going on.
% p_eqm_initial
% AggVars_initial
% [Params.B, Params.Bprime]
% 
% Policy_initial(2,1000:end,11:13)
% plot(shiftdim(Policy_initial(2,:,:),1))
% plot(shiftdim(Policy_initial(2,:,:),1)-(1:1:n_a)'.*ones(n_a,n_z))
% plot(sum(StationaryDist_initial,2))
% plot(cumsum(sum(StationaryDist_initial,2)))

%% Final stationary equilibrium
Params.phi=Params.phi_final;

[p_eqm_final,~,~]=HeteroAgentStationaryEqm_Case1(n_d, n_a, n_z, 0, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, [], [], [], GEPriceParamNames,heteroagentoptions);

Params.r=p_eqm_final.r;
[V_final,Policy_final]=ValueFnIter_Case1(n_d,n_a,n_z,d_grid,a_grid,z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, []);

% PolicyValues=PolicyInd2Val_Case1(Policy,n_d,n_a,n_s,d_grid,a_grid);

StationaryDist_final=StationaryDist_Case1(Policy_final,n_d,n_a,n_z,pi_z);

AggVars_final=EvalFnOnAgentDist_AggVars_Case1(StationaryDist_final, Policy_final, FnsToEvaluate,Params, [],n_d, n_a, n_z, d_grid, a_grid,z_grid);

% save ./SavedOutput/GuerrieriLorenzoni2017_final.mat Params p_eqm_final p_eqm_index_final MarketClearance_final V_final Policy_final StationaryDist_final AggVars_final n_d n_a n_z

%% Compute Annual GDP

% GL2017 describe Figure 1 as "Figure I shows the optimal values of
% consumption and labor supply as a function of the initial level of bond
% holdings". This is incorrect. The x-axis is not the level of bond holdings, it is
% bond-holdings as a fraction of annual output. I follow GL2017 in what
% they plot, adding the footnote to bottom of figure to explain what is
% actually graphed.
% In fact this is true of many of the Figures of GL2017, which label the
% x-axis as either b or B, but are in fact reporting b (or B) divided by
% annual output.
% To do this I need to compute annual GDP: quarterly output is y=theta*n.
% Following few lines do this (together with multiplication by 4 to make it annual)
FnsToEvaluateExtra.Output = @(d, aprime,a,z) d*z; % Output

QuarterlyOutput_initial=EvalFnOnAgentDist_AggVars_Case1(StationaryDist_initial, Policy_initial, FnsToEvaluateExtra,Params, [],n_d, n_a, n_z, d_grid, a_grid,z_grid,2);
QuarterlyOutput_final=EvalFnOnAgentDist_AggVars_Case1(StationaryDist_final, Policy_final, FnsToEvaluateExtra,Params, [],n_d, n_a, n_z, d_grid, a_grid,z_grid,2);

AnnualOutput_initial=4*QuarterlyOutput_initial.Output.Mean;
AnnualOutput_final=4*QuarterlyOutput_final.Output.Mean;

%% Figure 1
figure(1)
l_a=length(n_a);
l_z=length(n_z);
% We want to plot consumption, first set up a function that returns the value of consumption
Fig1FnsToEvaluate.Consumption = @(d, aprime,a,z,r, v, B, Bprime) GuerrieriLorenzoni2017_ConsumptionFn(d, aprime, a, z,r, v, B, Bprime); % Consumption

% Now evaluate it on the grid
ValuesOnGrid=EvalFnOnAgentDist_ValuesOnGrid_Case1(Policy_initial, Fig1FnsToEvaluate, Params, [], n_d, n_a, n_z, d_grid, a_grid, z_grid);
ConsumptionDecision=ValuesOnGrid.Consumption;

subplot(2,1,1); plot(a_grid,ConsumptionDecision(:,2),a_grid,ConsumptionDecision(:,8))
% legend('Mean','10th','25th','Median')
title({'Consumption'})
xlabel('Bond holdings')
% Labour supply is the d variable
subplot(2,1,2); plot(a_grid,d_grid(Policy_initial(1,:,2)),a_grid,d_grid(Policy_initial(1,:,8)))
legend('\theta^2','\theta^8');
title({'Labor Supply'})
xlabel('Bond holdings')
% saveas(gcf,['./SavedOutput/Graphs/GuerrieriLorenzoni_Figure1.pdf'])

% Now a version that reproduces GL2017 exactly
subplot(2,1,1); plot(a_grid/AnnualOutput_initial,ConsumptionDecision(:,2),a_grid/AnnualOutput_initial,ConsumptionDecision(:,8))
% legend('Mean','10th','25th','Median')
title({'Consumption'})
xlabel('Bond holdings as fraction of annual output')
ylabel('(Quarterly)')
% Labour supply is the d variable
subplot(2,1,2); plot(a_grid/AnnualOutput_initial,d_grid(Policy_initial(1,:,2)),a_grid/AnnualOutput_initial,d_grid(Policy_initial(1,:,8)))
legend('\theta^2','\theta^8');
title({'Labor Supply'})
xlabel('Bond holdings as fraction of annual output')
ylabel('(Quarterly)')
% saveas(gcf,['./SavedOutput/Graphs/GuerrieriLorenzoni_Figure1_GL2017.pdf'])

%% Figure 4
figure(4)
% Unclear from GL2017 paper what is plotted in the top panel.
% From their codes it is clear that it is 'mean based on the unconditional
% distribution of the exogenous shock' (note that this is not the same
% as the 'mean based on the distribution of agents at that asset level').
subplot(2,1,1); plot(a_grid,sum(a_grid(shiftdim(Policy_initial(2,:,:),1)).*pistar_z',2)-a_grid,a_grid,sum(a_grid(shiftdim(Policy_final(2,:,:),1)).*pistar_z',2)-a_grid)
subplot(2,1,2); plot(a_grid,sum(StationaryDist_initial,2), a_grid,sum(StationaryDist_final,2))
% saveas(gcf,['./SavedOutput/Graphs/GuerrieriLorenzoni_Figure4.pdf'])

% Now a version that reproduces GL2017 exactly
% Note: they are reporting x-axes that are 'as a fraction of initial annual
% output' for the final equilibrium results, not the final annual output.
% (Although in practice initial and final annual output are almost equal)
% For the top panel in fact everything (both axes) is as a fraction of
% initial annual output.
subplot(2,1,1); plot(a_grid/AnnualOutput_initial,(sum(a_grid(shiftdim(Policy_initial(2,:,:),1)).*pistar_z',2)-a_grid)/AnnualOutput_initial,a_grid/AnnualOutput_initial,(sum(a_grid(shiftdim(Policy_final(2,:,:),1)).*pistar_z',2)-a_grid)/AnnualOutput_initial)
subplot(2,1,2); plot(a_grid/AnnualOutput_initial,sum(StationaryDist_initial,2), a_grid/AnnualOutput_initial,sum(StationaryDist_final,2))
% saveas(gcf,['./SavedOutput/Graphs/GuerrieriLorenzoni_Figure4_GL2017.pdf'])


%%
% Free up space on gpu
clear ConsumptionDecision
% clear Policy_initial
clear PolicyValues PolicyValuesPermute
clear StationaryDist_final

% load ./SavedOutput/GuerrieriLorenzoni2017_final.mat V_final
% load ./SavedOutput/GuerrieriLorenzoni2017_initial.mat StationaryDist_initial

%% Compute the transition path
% For this we need the following extra objects: PricePathOld, PriceParamNames, ParamPath, ParamPathNames, T, V_final, StationaryDist_init
% (already calculated V_final & StationaryDist_init above)

% Number of time periods to allow for the transition (if you set T too low
% it will cause problems, too high just means run-time will be longer).
T=25

% We want to look at a one off unanticipated path of phi. ParamPath & PathParamNames are thus given by
ParamPath.phi=Params.phi_final*ones(T,1); % ParamPath is matrix of size T-by-'number of parameters that change over path'
temp=linspace(Params.phi_initial,Params.phi_final,7); ParamPath.phi(1:6)=temp(2:7); % At t=0, is inital stationary distribution, then falls over the following 6 periods to equal 0.525, remains there
% (the way ParamPath is set is designed to allow for a series of changes in the parameters)

% We need to give an initial guess for the price path on interest rates
PricePath0.r=[linspace(-0.01, p_eqm_final.r, floor(T/3))'; p_eqm_final.r*ones(T-floor(T/3),1)]; % For each price/parameter that will be determined in general eqm over transition path PricePath0 is matrix of size T-by-1

% Rewrite the aggregate variable to be next period bonds rather than
% current bonds as this is the actual timing of the decision which the
% interest rate (r) effects
FnsToEvaluate.Aprime = @(d, aprime,a,z) aprime; % Aggregate assets decisions
% Rewrite the General Eqm conditions as rules for updating the price
transpathoptions.GEnewprice=1; % If you do not do this the codes can still solve, but take much longer as they must figure out an updating rule for themselves.
GeneralEqmEqn.ZeroNetAssets = @(r,Aprime,Bprime) r-0.1*(Aprime-Bprime); % New interest rate is previous minus 0.1 times excess of bonds (I just guessed 0.1 as being pretty conservative, remember that the transition path will anyway do 0.1 new + 0.9 old price when updating at each iteration)

% [transpathoptions.GEnewprice=1 means that the GeneralEqmEqns should be
% expressed as how to generate a new guess for the price based on the
% current guess; transpathoptions.GEnewprice=0 means the GeneralEqmEqns
% should be expressed as for the standard general eqm conditions, namely
% equations that take the value of 0 in general eqm.]

% Now just run the TransitionPath_Case1 command (all of the other inputs
% are things we had already had to define to be able to solve for the
% initial and final equilibria)
transpathoptions.weightscheme=1
transpathoptions.verbose=1
PricePath=TransitionPath_Case1(PricePath0, ParamPath, T, V_final, StationaryDist_initial, n_d, n_a, n_z, pi_z, d_grid,a_grid,z_grid, ReturnFn,  FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, ReturnFnParamNames, FnsToEvaluateParamNames, GeneralEqmEqnParamNames, transpathoptions);

% save ./SavedOutput/GuerrieriLorenzoni2017_transpath1.mat PricePath n_d n_a n_z

% For later we will keep another copy
PricePath_Flex=PricePath;

%% Figure 3
% load ./SavedOutput/GuerrieriLorenzoni2017_initial.mat
% load ./SavedOutput/GuerrieriLorenzoni2017_final.mat
% load ./SavedOutput/GuerrieriLorenzoni2017_transpath1.mat 

Fig3FnsToEvaluate.output = @(d, aprime,a,z) d*z; % y_it=n_it*theta_it Note that since gov budget is balanced every period it neither adds nor subtracts (unemployment benefits + interest payments on B=lump-sum tax revenue)
Fig3FnsToEvaluate.debt = @(d, aprime,a,z) -a*(a<0); % debt is (minus of) negative assets

AggVars_initial=EvalFnOnAgentDist_AggVars_Case1(StationaryDist_initial, Policy_initial, Fig3FnsToEvaluate,Params, [],n_d, n_a, n_z, d_grid, a_grid,z_grid);
AggVarsPath=EvalFnOnTransPath_AggVars_Case1(Fig3FnsToEvaluate, [],PricePath, ParamPath, Params, T, V_final, StationaryDist_initial, n_d, n_a, n_z, pi_z, d_grid, a_grid,z_grid, DiscountFactorParamNames, ReturnFn, ReturnFnParamNames, transpathoptions);

Output_pch=([AggVars_initial.output.Mean; AggVarsPath(:,1)]-AggVars_initial.output.Mean)/AggVars_initial.output.Mean;

figure(3)
% Borrowing limit
subplot(2,2,1); plot(0:1:T,[Params.phi_initial; ParamPath])
title('borrowing constraint')
% household debt-to-GDP ratio
subplot(2,2,2); plot(0:1:T,[AggVars_initial.debt.Mean./AggVars_initial.output.Mean(1); AggVarsPath(:,2)./AggVarsPath(:,1)])
title('household debt-to-annual-GDP ratio')
% interest rate
subplot(2,2,3); plot(0:1:T,100*(((1+[p_eqm_initial; PricePath]).^4)-1)) % converts to annual rate by compounding (not just multiplying by 4)
title('annual interest rate')
% output
subplot(2,2,4); plot(0:1:T,100*Output_pch) % 100* to present one percentage point as 1
title('output')
% ylabel('percent deviation from inital output in stationary eqm')
% saveas(gcf,['./SavedOutput/Graphs/GuerrieriLorenzoni_Figure3.pdf'])

% Now a version that reproduces GL2017 exactly.
% Borrowing limit
subplot(2,2,1); plot(0:1:T,[Params.phi_initial; ParamPath]./AnnualOutput_initial)
title('borrowing constraint as fraction-of-initial-annual-output')
% household debt-to-GDP ratio
subplot(2,2,2); plot(0:1:T,[AggVars_initial.debt.Mean; AggVarsPath(1:end,2)]./AggVars_initial.output.Mean)
title('household debt-to-initial-annual-GDP ratio')
% interest rate
subplot(2,2,3); plot(0:1:T,100*4*[p_eqm_initial; PricePath]) % converts to annual rate by compounding (not just multiplying by 4)
title('annual interest rate')
% output
subplot(2,2,4); plot(0:1:T,100*Output_pch) % 100* to present one percentage point as 1
% ylabel('percent deviation from inital output in stationary eqm')
% saveas(gcf,['./SavedOutput/Graphs/GuerrieriLorenzoni_Figure3_GL2017.pdf'])


%% Figure 5
Fig5FnsToEvaluate.Consumption = @(d, aprime,a,z,r, v, B, Bprime) GuerrieriLorenzoni2017_ConsumptionFn(d, aprime, a, z,r, v, B, Bprime); % Consumption

AggVars_initial=EvalFnOnAgentDist_AggVars_Case1(StationaryDist_initial, Policy_initial, Fig5FnsToEvaluate,Params, [],n_d, n_a, n_z, d_grid, a_grid,z_grid);
AggVarsPath_GE=EvalFnOnTransPath_AggVars_Case1(Fig5FnsToEvaluate, [],PricePath, ParamPath, Params, T, V_final, StationaryDist_initial, n_d, n_a, n_z, pi_z, d_grid, a_grid,z_grid, DiscountFactorParamNames, ReturnFn, ReturnFnParamNames, transpathoptions);
% Only debt limit reduction path
UnchangedPricePath=p_eqm_initial*ones(T,1);
AggVarsPath_partial_onlydebtlimit=EvalFnOnTransPath_AggVars_Case1(Fig5FnsToEvaluate, [],UnchangedPricePath, ParamPath, Params, T, V_final, StationaryDist_initial, n_d, n_a, n_z, pi_z, d_grid, a_grid,z_grid, DiscountFactorParamNames, ReturnFn, ReturnFnParamNames, transpathoptions);
% Only interest rate path
UnchangedParamPath=Params.phi_initial*ones(T,1);
AggVarsPath_partial_onlyinterestrate=EvalFnOnTransPath_AggVars_Case1(Fig5FnsToEvaluate, [],PricePath,PricePathNames, UnchangedParamPath, ParamPathNames, Params, T, V_final, StationaryDist_initial, n_d, n_a, n_z, pi_z, d_grid, a_grid,z_grid, DiscountFactorParamNames, ReturnFn, ReturnFnParamNames, transpathoptions);

AggVarsPath_GE_pch=100*(AggVarsPath_GE-AggVars_initial.Consumption.Mean)./AggVars_initial.Consumption.Mean;
AggVarsPath_partial_onlydebtlimit_pch=100*(AggVarsPath_partial_onlydebtlimit-AggVars_initial.Consumption.Mean)./AggVars_initial.Consumption.Mean;
AggVarsPath_partial_onlyinterestrate_pch=100*(AggVarsPath_partial_onlyinterestrate-AggVars_initial.Consumption.Mean)./AggVars_initial.Consumption.Mean;

figure(5)
plot(0:1:T,[0;AggVarsPath_GE_pch], 0:1:T,[0;AggVarsPath_partial_onlydebtlimit_pch], 0:1:T,[0;AggVarsPath_partial_onlyinterestrate_pch])
title('Consumption Response Deviation')
legend('General eqm response','Partial eqm response to debt limit reduction','Partial eqm response to interest rate changes')
% saveas(gcf,['./SavedOutput/Graphs/GuerrieriLorenzoni_Figure5.pdf'])


%% Figure 6
% Unlike the other figures relating to the transition, which just require
% the distribution of agents at each time period, this figure requires
% following individuals along the transition path based on where they
% started. Hence will use SimPanelValues_TransPath_Case1() rather than
% EvalFnOnTransPath_AggVars_Case1(); actually the later could/should also be used
% here, but want to show the different options available as part of VFI Toolkit.

% For each of the four transition paths simulate 100 paths drawing from the relevant initial percentile, then take mean.
% (Alternative would be a much bigger panel based on drawing from the actual inital stationary distribution, then 
% just restrict panel based on percentiles and take means, but this would involve much more computation time)
simoptions.numbersims=100;

% We are going to want the values for consumption
Fig6FnsToEvaluate.Consumption = @(d, aprime,a,z,r, v, B, Bprime) GuerrieriLorenzoni2017_ConsumptionFn(d, aprime, a, z,r, v, B, Bprime); % Consumption

% First, figure out the asset values that correspond to the percentiles
assetdist=cumsum(sum(StationaryDist_initial,2));
[~, prctileindexes]=min(abs(assetdist-(1:1:100)/100)); % prctilevalues_doublecheck should be approx equal to prctilevalues
% prctilevalues_doublecheck=assetdist(prctileindexes); % should give [0.01,0.02, etc up to 1]

% 1st percentile (note, the 1st percentile are going to be those closest to the borrowing constraint)
% Set up the appropriate inital distribution and simulate a panel data set (over the transition) from this.
InitialDist_1stpercentile=zeros(n_a,n_z,'gpuArray');
InitialDist_1stpercentile(prctileindexes(1),:)=StationaryDist_initial(prctileindexes(1),:)./sum(StationaryDist_initial(prctileindexes(1),:)); % Normalized version of agents holding the 1st-percentile amount of assets, I make sure they have the appropriate distribution over the exogenous shock dimension.
% Everything else is just completely standard
AggVars_initial=EvalFnOnAgentDist_AggVars_Case1(InitialDist_1stpercentile, Policy_initial, Fig6FnsToEvaluate,Params, Fig6FnsToEvaluateParamNames,n_d, n_a, n_z, d_grid, a_grid,z_grid);
SimPanelValues=SimPanelValues_TransPath_Case1(PricePath, ParamPath, T, V_final, InitialDist_1stpercentile, n_d, n_a, n_z, pi_z, d_grid,a_grid,z_grid, ReturnFn, Fig6FnsToEvaluate, Params, DiscountFactorParamNames, ReturnFnParamNames, Fig6FnsToEvaluateParamNames, transpathoptions,simoptions);
% The line in the figure is just the mean for each time period of these (I am guessing), but expressed as 
% percent deviation from steady state. [Not obvious if I should take mean and then percent deviation, or 
% take percent deviation and then mean; have gone with the former.]
SimPanelValues=shiftdim(SimPanelValues,1);
% Fig6_1stPercentileTrace=(mean(SimPanelValues,2)-mean(SimPanelValues(1,:)))/mean(SimPanelValues(1,:));
Fig6_1stPercentileTrace=(mean(SimPanelValues,2)-AggVars_initial.Consumption.Mean)/AggVars_initial.Consumption.Mean;

% Now just repeat for 10th, 20th and 50th percentiles
InitialDist_10thpercentile=zeros(n_a,n_z,'gpuArray');
InitialDist_10thpercentile(prctileindexes(10),:)=StationaryDist_initial(prctileindexes(10),:)./sum(StationaryDist_initial(prctileindexes(10),:)); % Normalized version of agents holding the 1st-percentile amount of assets, I make sure they have the appropriate distribution over the exogenous shock dimension.
AggVars_initial=EvalFnOnAgentDist_AggVars_Case1(InitialDist_10thpercentile, Policy_initial, Fig6FnsToEvaluate,Params, Fig6FnsToEvaluateParamNames,n_d, n_a, n_z, d_grid, a_grid,z_grid);
SimPanelValues=SimPanelValues_TransPath_Case1(PricePath, ParamPath, T, V_final, InitialDist_10thpercentile, n_d, n_a, n_z, pi_z, d_grid,a_grid,z_grid, ReturnFn, Fig6FnsToEvaluate, Params, DiscountFactorParamNames, ReturnFnParamNames, Fig6FnsToEvaluateParamNames, transpathoptions,simoptions);
SimPanelValues=shiftdim(SimPanelValues,1);
Fig6_10thPercentileTrace=(mean(SimPanelValues,2)-AggVars_initial.Consumption.Mean)/AggVars_initial.Consumption.Mean;
InitialDist_20thpercentile=zeros(n_a,n_z,'gpuArray');
InitialDist_20thpercentile(prctileindexes(20),:)=StationaryDist_initial(prctileindexes(20),:)./sum(StationaryDist_initial(prctileindexes(20),:)); % Normalized version of agents holding the 1st-percentile amount of assets, I make sure they have the appropriate distribution over the exogenous shock dimension.
AggVars_initial=EvalFnOnAgentDist_AggVars_Case1(InitialDist_20thpercentile, Policy_initial, Fig6FnsToEvaluate,Params, Fig6FnsToEvaluateParamNames,n_d, n_a, n_z, d_grid, a_grid,z_grid);
SimPanelValues=SimPanelValues_TransPath_Case1(PricePath, ParamPath, T, V_final, InitialDist_20thpercentile, n_d, n_a, n_z, pi_z, d_grid,a_grid,z_grid, ReturnFn, Fig6FnsToEvaluate, Params, DiscountFactorParamNames, ReturnFnParamNames, Fig6FnsToEvaluateParamNames, transpathoptions,simoptions);
SimPanelValues=shiftdim(SimPanelValues,1);
Fig6_20thPercentileTrace=(mean(SimPanelValues,2)-AggVars_initial.Consumption.Mean)/AggVars_initial.Consumption.Mean;
InitialDist_50thpercentile=zeros(n_a,n_z,'gpuArray');
InitialDist_50thpercentile(prctileindexes(50),:)=StationaryDist_initial(prctileindexes(50),:)./sum(StationaryDist_initial(prctileindexes(50),:)); % Normalized version of agents holding the 1st-percentile amount of assets, I make sure they have the appropriate distribution over the exogenous shock dimension.
AggVars_initial=EvalFnOnAgentDist_AggVars_Case1(InitialDist_50thpercentile, Policy_initial, Fig6FnsToEvaluate,Params, Fig6FnsToEvaluateParamNames,n_d, n_a, n_z, d_grid, a_grid,z_grid);
SimPanelValues=SimPanelValues_TransPath_Case1(PricePath, ParamPath, T, V_final, InitialDist_50thpercentile, n_d, n_a, n_z, pi_z, d_grid,a_grid,z_grid, ReturnFn, Fig6FnsToEvaluate, Params, DiscountFactorParamNames, ReturnFnParamNames, Fig6FnsToEvaluateParamNames, transpathoptions,simoptions);
SimPanelValues=shiftdim(SimPanelValues,1);
Fig6_50thPercentileTrace=(mean(SimPanelValues,2)-AggVars_initial.Consumption.Mean)/AggVars_initial.Consumption.Mean;

figure(6)
plot(0:1:T-1, [0; Fig6_1stPercentileTrace], 0:1:T-1, [0;Fig6_10thPercentileTrace], 0:1:T-1, [0;Fig6_20thPercentileTrace], 0:1:T-1, [0;Fig6_50thPercentileTrace])
title('Consumption Response by Percentile in Initial Wealth Distribution')
legend('1st percentile','10th percentile','20th percentile','50th percentile')
% saveas(gcf,['./SavedOutput/Graphs/GuerrieriLorenzoni_Figure6.pdf'])


%% Figure 7
Fig7FnsToEvaluate.Employment = @(d, aprime,a,z) d; %n_it in notation of GL2017

Employment_initial=EvalFnOnAgentDist_AggVars_Case1(StationaryDist_initial, Policy_initial, Fig7FnsToEvaluate,Params, [],n_d, n_a, n_z, d_grid, a_grid,z_grid,2);
EmploymentPath_GE=EvalFnOnTransPath_AggVars_Case1(Fig7FnsToEvaluate, [],PricePath, ParamPath, Params, T, V_final, StationaryDist_initial, n_d, n_a, n_z, pi_z, d_grid, a_grid,z_grid, DiscountFactorParamNames, ReturnFn, ReturnFnParamNames,transpathoptions);

EmploymentPath_GE_pch=100*(EmploymentPath_GE-Employment_initial.Employment.Mean)./Employment_initial.Employment.Mean;

figure(7)
plot(0:1:T,[0; EmploymentPath_GE_pch])
title('Employment Response')
% saveas(gcf,['./SavedOutput/Graphs/GuerrieriLorenzoni_Figure7.pdf'])


%% Sticky Wages
% What GL2017 refer to as sticky wages actually enters the model as an increase in
% the utility of leisure. (They justify this as decreasing hours worked
% without decreasing the actual real-age-per-unit-of-time-worked which is
% intended to capture a decrease in labour demand by (unmodelled) firms.)

% Only two changes: first is that the return function needs to be modified
% for this by including omega. (This was actually already implemented in baseline for this code, so no further change.)
% Second is general equilibrium conditions, these now need to enforce r>=0
% and determine omega>0 as that which is needed to acheive that r=0 as a general
% eqm outcome, whenever omega=0 would otherwise mean r<0.
% Then compute transition path (omega=0 in initial and final states, so these remain unchanged).
% The key is all in writing the general equilibrium conditions in right way
% so that omega=0, except when this would lead to r<0, and in these cases
% need to pick omega so that r=0.

% We need to give an initial guess for the price path on interest rates and
% omega. Lets just start with the flexible prices general eqm path for
% interest rates, and with omega=0 for all periods as our initial guess.
PricePath0.r=PricePath.r;
PricePath0.omega=zeros(T,1);

% Rewrite the General Eqm conditions as rules for updating the price
transpathoptions.GEnewprice=1; % If you do not do this the codes can still solve, but take much longer as they must figure out an updating rule for themselves.
GeneralEqmEqnParamNames(1).Names={'Bprime'};
% Only change to this is to enforce that will only try to decrease interest rate when r>=0 is satisfied, never when it is less than zero. When r<0
% instead as 0.01 (note that transition path codes will use 0.9*old+0.1*new anyway, so essentially adding 0.001, a tenth of one percentage point, 
% to interest rate)
GeneralEqmEqn_1 = @(AggVars,p,Bprime) (p(1)>=0)*(p(1)-0.1*(AggVars(1)-Bprime))+(p(1)<0)*(p(2)>0)*(p(1)+2*abs(p(1))+0.0001); % New interest rate is previous minus 0.1 times excess of bonds (I just guessed 0.1 as being pretty conservative, remember that the transition path will anyway do 0.1 new + 0.9 old price when updating at each iteration)
GeneralEqmEqnParamNames(2).Names={};
GeneralEqmEqn_2 = @(AggVars,p,Bprime) (p(1)<0)*(p(2)+0.003)+(p(1)>=0)*(p(2)>0.003)*(p(2)-0.003); % If r>=0 then send omega towards zero (as in this case it returns max{0,omega-0.005}). r<0 then increase omega (in this case it returns omega+0.005). (Note: the choice of +0.03 in principle will be good or bad depending on the update rule for the transition path; given that weight on old price is 0.9 this will shift omega up by (1-0.9)*0.02 whenever the interest rate is negative)
GeneralEqmEqns={GeneralEqmEqn_1,GeneralEqmEqn_2};
% Remark: This approach to how to update omega is different from that used
% by GL2017. They instead force the ZLB on the interest rate r, and then use
% omega to ensure that goods markets clear (Y=C; remember this is an
% endowment economy). omega is updated based on how far the model is from
% goods market clearance.
% Remark: The above updating rules took a few attempts to come up with (the r>=0
% constraint is a 'knife-edge' and so it was not trivial to find something
% that settles down rather than 'jump' back and forth between r<0 and r>0).

% [transpathoptions.GEnewprice=1 means that the GeneralEqmEqns should be
% expressed as how to generate a new guess for the price based on the
% current guess; transpathoptions.GEnewprice=0 means the GeneralEqmEqns
% should be expressed as for the standard general eqm conditions, namely
% equations that take the value of 0 in general eqm.]

% Now just run the TransitionPath_Case1 command (all of the other inputs
% are things we had already had to define to be able to solve for the
% initial and final equilibria)
transpathoptions.weightscheme=1; % This is anyway the default
transpathoptions.verbose=1;
transpathoptions.tolerance=10^(-4); % will run until r and omega settle to four digits
transpathoptions % show what the options have been set to
PricePath_NK=TransitionPath_Case1(PricePath0, ParamPath, T, V_final, StationaryDist_initial, n_d, n_a, n_z, pi_z, d_grid,a_grid,z_grid, ReturnFn,  FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, ReturnFnParamNames, FnsToEvaluateParamNames, GeneralEqmEqnParamNames, transpathoptions);

%% Figure 8 (I do an additional Figure 17 that shows the 'wedges' and employment)
Fig8FnsToEvaluate.output = @(d, aprime,a,z) d*z; % y_it=n_it*theta_it Note that since gov budget is balanced every period it neither adds nor subtracts (unemployment benefits + interest payments on B=lump-sum tax revenue)
Fig8FnsToEvaluate.employment = @(d, aprime,a,z) d; %n_it in notation of GL2017

AggVars_initial=EvalFnOnAgentDist_AggVars_Case1(StationaryDist_initial, Policy_initial, Fig8FnsToEvaluate,Params, [],n_d, n_a, n_z, d_grid, a_grid,z_grid);
% Difference between following two lines is PricePath vs PricePath_NK
AggVarsPath_Flex=EvalFnOnTransPath_AggVars_Case1(Fig8FnsToEvaluate, [],PricePath, ParamPath, Params, T, V_final, StationaryDist_initial, n_d, n_a, n_z, pi_z, d_grid, a_grid,z_grid, DiscountFactorParamNames, ReturnFn, ReturnFnParamNames, transpathoptions);
AggVarsPath_NK=EvalFnOnTransPath_AggVars_Case1(Fig8FnsToEvaluate, [],PricePath_NK, ParamPath, Params, T, V_final, StationaryDist_initial, n_d, n_a, n_z, pi_z, d_grid, a_grid,z_grid, DiscountFactorParamNames, ReturnFn, ReturnFnParamNames, transpathoptions);

OutputPath_pch_Flex=(AggVarsPath_Flex(:,1)-AggVars_initial.output.Mean)./AggVars_initial.output.Mean;
OutputPath_pch_NK=(AggVarsPath_NK(:,1)-AggVars_initial.output.Mean)./AggVars_initial.output.Mean;
% My extras
EmploymentPath_pch_Flex=(AggVarsPath_Flex(:,2)-AggVars_initial.employment.Mean)./AggVars_initial.employment.Mean; % This has already been calculated.
EmploymentPath_pch_NK=(AggVarsPath_NK(:,2)-AggVars_initial.employment.Mean)./AggVars_initial.employment.Mean;

figure(8)
% interest rate
subplot(2,1,1); plot(0:1:T,4*100*[p_eqm_initial; PricePath_Flex], 0:1:T,4*100*[p_eqm_initial; PricePath_NK(:,1)])
title('annual interest rate')
legend('flex price', 'NK fix price (ZLB)')
% output
subplot(2,1,2); plot(0:1:T,100*[0; OutputPath_pch_Flex],0:1:T,100*[0; OutputPath_pch_NK])
title('output')
% ylabel('percent deviation from inital output in stationary eqm')
% saveas(gcf,['./SavedOutput/Graphs/GuerrieriLorenzoni_Figure8.pdf'])

saveas(gcf,['./SavedOutput/Graphs/GuerrieriLorenzoni_Figure8_Tauchenq3.pdf'])

% Not actually part of GL2017, but let's take a look at the path for omega (the wedge on labour supply that 'enforces' the zero lower bound on interest rates)
figure(17)
% wedges
subplot(2,1,1); plot(0:1:T,[0; PricePath_NK(:,2)])
title('omega (NK wedge on wages)')
% employment
subplot(2,1,2); plot(0:1:T,100*[0; EmploymentPath_pch_Flex],0:1:T,100*[0; EmploymentPath_pch_NK])
title('employment')
legend('flex price', 'NK fix price (ZLB)')
% ylabel('percent deviation from inital output in stationary eqm')
% saveas(gcf,['./SavedOutput/Graphs/GuerrieriLorenzoni_ExtraFigure1.pdf'])





