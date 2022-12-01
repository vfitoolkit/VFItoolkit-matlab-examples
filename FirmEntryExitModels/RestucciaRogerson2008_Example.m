% Replication of Restuccia & Rogerson (2008) - Policy Distortions and Aggregate Productivity with Heterogeneous Establishments
%
% Note: rho is defined differently here to in the paper. Here it does not
% include the (exogenous) conditional probability of firm exit. In RR2008
% it does. This is done to emphasize how to include exit/death in models
% more generally when using VFI Toolkit.
%
% For explanation of model, the calibration choices, and what the policy
% experiments are aiming to understand, see paper.


%% Set some basic variables

n_s=100; % Firm-specific Productivity level
n_tau=3; % Negative (subsidy), zero, and positive (tax).

%Parameters
Params.beta=0.96; % Discount rate
% Firms
Params.alpha=0.283;
Params.gamma=0.567;
Params.delta=0.08; % Depreciation rate of physical capital
Params.cf=0; % Fixed cost of production
% Firm entry and exit
Params.ce=1; % Fixed cost of entry (this is a normalization)
Params.lambda=0.1; % Probability of firm exit

Params.Ne=1; % RR2008 call this E, the mass of new potential new entrant distribution.

% The actual 'distortionary policy rate'
Params.taurate=0; % This is the rate for the tax.
Params.subsidyrate=0; % This is the rate for the subsidy.

Params.w=1; % This is just an initial guess since w (wage) is determined in general eqm.

%% Grid on the productivity level, s, as well as the distribution over these h(s).
% The grid is here created by the next few lines of text.
% The original numbers which can be found in the codes provided by Restuccia & Rogerson (2008) for their model. 
% The original source is Rossi-Hansberg and Wright (AER, 2007) All industries for 2000,
% which in turn is created from data from the US Census of Businesses. The
% following commented out line does the import that creates what is stored here as text.
% load ./OriginalCodes/files/establishment_dist.txt
establishment_dist=[4, 0.482359949; 9, 0.21665686; 14, 0.092591623;...
    19, 0.049473226; 24, 0.031397291; 29, 0.021542855; 34, 0.015801151; 39, 0.011847888; 44, 0.009307797;...
    49, 0.007457313; 59, 0.011325843; 69, 0.007940958; 79, 0.005988763; 89, 0.004677303; 99, 0.003753013;...
    124, 0.00683308; 149, 0.004509741; 174, 0.003126242; 199, 0.002215598; 224, 0.001685936; 249, 0.001293212;...
    299, 0.001867145; 349, 0.001285596; 399, 0.000873831; 449, 0.000649621; 499, 0.000511097; 599, 0.000705634;...
    699, 0.000463018; 799, 0.000331475; 899, 0.0002469; 999, 0.000184065; 1249, 0.000311482; 1499, 0.000191522;...
    1749, 0.000130432; 1999, 9.10802E-05; 2249, 7.14044E-05; 2499, 5.42673E-05; 2999, 6.96589E-05; ...
    3499, 4.42707E-05; 3999, 3.09419E-05; 4499, 1.96759E-05; 4999, 1.60263E-05; 9999, 5.06178E-05; ...
    10000, 1.45982E-05];
% establishment_dist(:,1); %(upper range of number of employees)
% establishment_dist(:,2); %(fraction of establishment in each group)

% Note: Table 1 of RR2008 reports an 's range', this is actually implicit
% from the numbers in establishment_dist and the parameters gamma and alpha.
s_grid=exp(linspace(log(1),log(establishment_dist(end,1)^(1-Params.gamma-Params.alpha)),n_s));
% The stationary distribution, called h(s) by RR2008, will here be called pistar_s
cumsum_pistar_s=interp1(establishment_dist(:,1),cumsum(establishment_dist(:,2)),s_grid.^(1/(1-Params.gamma-Params.alpha)));
% The first few points of s_grid have to be extrapolated (as outside the
% range of the actual data). Following line implements this, and similar
% for the max value. I just divide equally across these first few (which is same as RR2008 do).
temp=~isnan(cumsum_pistar_s);
tempind=find(temp,1,'first');
cumsum_pistar_s(1:tempind)=cumsum((1/tempind)*cumsum_pistar_s(tempind)*ones(1,tempind));
cumsum_pistar_s(end)=1;
pistar_s=(cumsum_pistar_s-[0,cumsum_pistar_s(1:end-1)])';
% Note: This stationary distribution is created by interpolating establishment_dist(:,2)
% (interpolation of the cumulative distribution tends to perform better in
% practice than interpolation of the pdf, as cdf tends to be smoother)
% Here this is done using 1D interpolation (the matlab default spline interpolation), RR2008 do this manually by evenly
% dividing the data across the s_grid points that are 'in-between' establishment_dist points.
% (You can see the difference if you use original RR2008 codes to create 'ns' and 'hs', and then run following line
% plot(establishment_dist(:,1),cumsum(establishment_dist(:,2)),s_grid.^(1/(1-Params.gamma-Params.alpha)), cumsum(pistar_s), z, cumsum(hs))
% The difference is not large, but the use of interpolation clearly provide better approximation.

% Note: if you want the actual establishment sizes implicit in s_grid to
% compare to data just run: s_grid.^(1/(1-Params.gamma-Params.alpha))

% Comment: The calibration of the model is awkward, in the sense that it is
% the distribution of potential entrants that is calibrated to match the US
% Census Bureau data on existing firms, while conceptually it is the
% stationary distribution of (exisiting) firms that should match the US
% Census Bureau data on existing firms. Although this later would require a
% much more complicated calibration procedure (moment matching based on simulation).

%% Grids
% Grids for d variables: there are none (because all the decisions being
% made by the firms have closed form solutions, so there are none that the
% toolkit actually needs to solve for)
n_d=0;
d_grid=[];
% Grids for a variables: there are none
n_a=1; % Still need to declare an a variable, the toolkit is not really designed for situation of none. Since it only has a single value and won't really be used anywhere it won't actually do anything anyway.
a_grid=1;
% Grids for z variables: there are two
% s_grid
tau_grid=[-1; 0; 1]; % Negative will be subsidy, zero, positive will be tax
n_z=[n_s,length(tau_grid)];
z_grid=[s_grid'; tau_grid];
pi_z=eye(prod(n_z),prod(n_z)); % Ones on the diagonal, zeros elsewhere. These are essentially permanent states, and there is zero probability of transition.

%% Figure 1

figure(1)
semilogx(s_grid.^(1/(1-Params.gamma-Params.alpha)),cumsum(pistar_s),'-',establishment_dist(:,1),cumsum(establishment_dist(:,2)),'o','LineWidth',2)
title('Distribution of Plants by Employment')
set(gca,'XTick',[1 10 100 1000 10000])
set(gca,'XTickLabel',{'1','10','100','1,000','10,000'})
xlabel('Number of Employees (log scale)')
ylabel('Cummulative Distribution of Establishments')
legend('Model','Data','Location','NorthWest')
axis([0 14000 0 1.1])
% saveas(gcf,['./RestucciaRogerson2008_Figure1.pdf'])
% Remark: the 'sagging curves' of model between data points is a result of
% the graph being done with log scale, while the interpolation is done
% based on (linear scale) number of employees directly.

%% Odds and ends

% Because asset markets are complete the households problem leads to the
% standard result that the households can just be treated as a
% Representative Household, and equilibrium in the asset markets requires
% that the interest rate, i (referred to in RR2008 paper as R, but here
% follow their codes and refer to as i)
Params.i=1/Params.beta-1; % This is standard general eqm result in complete market models, comes from consumption euler eqn together with requirements of stationary eqm.
% The net return to capital in equilibrium will thus be
Params.r=Params.i+Params.delta; % That the gross return is just 1/beta-1 and equals i (that the gross return to capital equals the interest rate is a requirement of capital market clearance in model)

% This (i=1/beta-1) means HHs supply amount of K that clears market, and
% labour supply is perfectly inelastic (leisure is not valued by HHs) 
% means that any value of w clears labour market.

%% Solve the Value Function
% Entry is irrelevant to the value function problem. Exit is relevant, and treated differently based on whether it is endogenous or exogenous.

% For exogenous exit, you simply need to include the 'conditional surivival probability' as another 'DiscountFactorParamNames'
Params.oneminuslambda=1-Params.lambda; % This is now the conditional probability of survival.
Params.rho=1/(1+Params.i); % The factor at which firms discount the future depends on both the risk-free interest rate and the risk/probability of exit
DiscountFactorParamNames={'rho','oneminuslambda'};

ReturnFn=@(aprime, a, z1,z2,w,r,alpha,gamma,taurate,subsidyrate,cf) RestucciaRogerson2008_ReturnFn(aprime, a, z1,z2,w,r,alpha,gamma,taurate,subsidyrate,cf);

% Check that everything is working so far by solving the value function
[V,Policy]=ValueFnIter_Case1(n_d,n_a,n_z,[],a_grid,z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, []);

% If you wanted to look at the value fn
% figure(2)
% surf(shiftdim(V,1))
% title('Value fn')

%% Simulate Agent Distribution, with Entry and Exit
% I use different notation to RR2008. The pdf of potential
% entrants I call upsilon (they used g). I call the 'conditional entry
% decision' ebar (they used xbar). I call the 'mass' by which this is
% multiplied Ne (they used E). Note that Ne is not the actual mass of
% entrants because of the role of ebar.
% So the actual distribution of entrants is then Ne*upsilon.*ebar (note that
% this is not a pdf). The actual mass of entrants can be calculated by
% repeatedly taking the sum of Ne*upsilon.*ebar.

% Need to 'activate' these options
simoptions.agententryandexit=1
% Note: it is possible to do just entry or just exit.
% To do just entry, simply set CondlProbOfExit=0;
% To do just exit, simply set DistOfNewAgents=0;
% (It is not envisaged that you are likely to want either of entry or exit without the other)
% Note that from the perspective of the simulation of agent distribution
% whether these are endogenous or exogenous entry/exit decisions is irrelevant.

pistar_tau=[0;1;0];

% Because they are not a default part of agent simulation, you need to pass the entry/exit aspects as part of simoptions.
EntryExitParamNames.DistOfNewAgents={'upsilon'};
%Params.upsilon=kron(pistar_tau,pistar_s); % Note: these should be in 'reverse order'
Params.upsilon=pistar_s.*(pistar_tau');
EntryExitParamNames.CondlEntryDecisions={'ebar'};
Params.ebar=ones([n_a,n_z]); % Takes value of one for enter, zero for not-enter. This is just an initial guess as the actual decisions are determined as part of general equilibrium.
% Note: VFI Toolkit requires the DistOfNewAgents to be a pdf (so unit mass), and then uses
% the 'MassOfNewAgents' to understand how many there will be entering relative to existing agents.
EntryExitParamNames.MassOfNewAgents={'Ne'}; % This is implied by the need for there to be a stationary distribution, so the mass of firm entry must equal the mass of firm exit, which is lambda*1

EntryExitParamNames.CondlProbOfSurvival={'oneminuslambda'};
% This conditional probability can be a state-dependent parameter, in which case this input would be a vector or matrix, etc.

% Check that everything is working so far by solving the simulation of
% agent distribution to get the stationary distribution.
StationaryDist=StationaryDist_Case1(Policy,n_d,n_a,n_z,pi_z, simoptions, Params, EntryExitParamNames);

% If you wanted to look at the pdf:
% surf(shiftdim(StationaryDist.pdf,1))

%% General equilibrium conditions (and aggregates needed to evaluation them)

% Endogenous entry adds a general equilibrium condition, specifically a 'free entry' condition.
% The only additional information needed for the general equilibrium in terms of entry
% and exit is to know which if any of the general equilibrium conditions are 'free entry'
% conditions and so need to be evaluated on potential entrants
% distribution.

%Use the toolkit to find the equilibrium prices
GEPriceParamNames={'ebar','ce'};

% Caclulating the general equilibrium does not require any aggregate variables
FnsToEvaluate=struct();

heteroagentoptions.specialgeneqmcondn={'condlentry','entry'};
% A 'condlentry' general equilibrium condition will take values of greater
% than zero for firms that decide to enter, less than zero for first that
% decide not to enter (or more accurately, after entry decision they draw
% their state, and then decide to cancel/abort their entry).

GEPriceParamNames={'w'}; 
% Note that this parameter does not directly appear in any of the general eqm conditions, only indirectly by it's effect on the value fn.
% Note that 'ebar', the conditional entry decision, is also determined as part of general eqm, and so in some sense is a general eqm parameter. 
% But since this is unavoidably the case for conditional entry there is no need to declare it.

GeneralEqmEqns.CondlEntry = @(ValueFn,beta) beta*ValueFn-0; % % Conditional entry condition
% CondlEntry: first input must be ValueFn, then any parameters
GeneralEqmEqns.Entry = @(EValueFn,beta,ce) beta*EValueFn-ce; % Free entry conditions (expected returns equal zero in eqm); note that the first 'General eqm price' is ce, the fixed-cost of entry.
% CondlEntry: EValueFn is a 'reserved' name, you can input it and any parameters

% In principle, 'Ne', the mass of (potential) new entrants is also a
% parameter to be determined in general equilibrium, and hence would also
% appear in GEPriceParamNames. The relevant condition would be labour
% market clearance 1-Nbar=0. The labour supply in the model is equal to 1
% (households are endowed with a unit endowment of labour, and this is
% supplied perfectly inelastically because households do not value
% leisure). Nbar is the labour demand and comes from integral of nbar over
% the distribution of firms. To implement this we would need to add the following:
% GEPriceParamNames={'w','Ne'};
% FnsToEvaluate.nbar = @(aprime,a,z1,z2,mass,alpha,gamma,r,w,taurate) (((1-taurate*z2)*z1*gamma)/w)^(1/(1-gamma)) *((alpha/r)^((1-gamma)/(1-gamma-alpha)) *(gamma/w)^(gamma/(1-gamma-alpha)) *(z1*(1-taurate*z2))^(1/(1-alpha-gamma)))^(alpha/(1-gamma)); % which evaluates to Nbar in the aggregate
% GeneralEqmEqn.LabourMarket = @(nbar) 1-nbar;
% heteroagentoptions.specialgeneqmcondn={'condlentry','entry',0};
% (Note, together with the previous two this gives us three general eqm condtions)
% Because RR2008 model is linear in Ne (not the case for most models of entry), we 
% can instead simply impose this afterwards, by setting Ne=1/Nbar. This is done
% here just after computing the general equilibrium.
% Obviously endogenizing labour supply would change this, and this 'short-cut' version would no longer work.

heteroagentoptions.verbose=1;
n_p=0;
disp('Calculating price vector corresponding to the stationary eqm')
% tic;
% NOTE: EntryExitParamNames has to be passed as an additional input compared to the standard case.
[p_eqm,p_eqm_index, GeneralEqmCondition]=HeteroAgentStationaryEqm_Case1(0, n_a, n_z, n_p, pi_z, [], a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, [], [], [], GEPriceParamNames,heteroagentoptions, simoptions, [], EntryExitParamNames);
% findeqmtime=toc
Params.w=p_eqm.w;
Params.ebar=p_eqm.ebar;

% RR2008 gets baseline solution of w=1.8955, xbar is all ones (based on running the code they provide on Review of Economic Dynamics website).
% I get w=1.9074, ebar is all ones.

% Calculate some things in the general eqm
[V,Policy]=ValueFnIter_Case1(n_d,n_a,n_z,[],a_grid,z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, []);
StationaryDist=StationaryDist_Case1(Policy,n_d,n_a,n_z,pi_z, simoptions, Params, EntryExitParamNames);

% Impose the labour market clearance, which involves calculating Ne. See comments above (about 10-20 lines above).
FnsToEvaluate.nbar = @(aprime,a,z1,z2,agentmass,alpha,gamma,r,w,taurate,subsidyrate) ((1-((z2>=0)*taurate+(z2<0)*subsidyrate)*z2)*z1)^(1/(1-alpha-gamma)) *(alpha/r)^(alpha/(1-gamma-alpha)) *(gamma/w)^((1-alpha)/(1-gamma-alpha)); % which evaluates to Nbar in the aggregate
AggValues=EvalFnOnAgentDist_AggVars_Case1(StationaryDist, Policy, FnsToEvaluate, Params, [], n_d, n_a, n_z, d_grid, a_grid, z_grid, [], simoptions, EntryExitParamNames);
InitialNe=Params.Ne;
Params.Ne=1/AggValues.nbar.Aggregate; % AggValues is presently equal to Nbar. This line is imposing/satisfying the labour market clearance condition.
StationaryDist.mass=StationaryDist.mass*(Params.Ne/InitialNe); % Take advantage of linearity of the stationary distribution in new entrants distribution.

%% Table 1
% Not exactly replicating anything as this is just the parameter values...

%Table 1
FID = fopen('./RestucciaRogerson2008_Table1.tex', 'w');
fprintf(FID, 'Benchmark calibration to US data \\\\ \n');
fprintf(FID, '\\begin{tabular*}{1.00\\textwidth}{@{\\extracolsep{\\fill}}lll} \n \\hline \n');
fprintf(FID, 'Parameters & Value & Target \\\\ \\hline \n');
fprintf(FID, '$\\alpha$ & %8.3f & Capital income share \\\\ \n', Params.alpha);
fprintf(FID, '$\\gamma$ & %8.3f & Labour income share \\\\ \n', Params.gamma);
fprintf(FID, '$\\beta$ & %8.2f & Real rate of return \\\\ \n', Params.beta);
fprintf(FID, '$\\delta$ & %8.2f & Investment to output rate \\\\ \n', Params.delta);
fprintf(FID, '$c_e$ & %8.1f & Normalization \\\\ \n', Params.ce);
fprintf(FID, '$c_f$ & %8.1f & Benchmark rate \\\\ \n', Params.cf);
fprintf(FID, '$\\lambda$ & %8.1f & Annual exit rate \\\\ \n', Params.lambda);
fprintf(FID, '$s$ range & [%d, %8.2f] & Relative establishment sizes \\\\ \n', s_grid(1), s_grid(end));
fprintf(FID, '$h(s)$ & see Fig. 1 & Size distribution of establishments \\\\ \n');
fprintf(FID, '\\hline \n \\end{tabular*} \n');
% fprintf(FID, '\\begin{minipage}[t]{1.00\\textwidth}{\\baselineskip=.5\\baselineskip \\vspace{.3cm} \\footnotesize{ \n');
% fprintf(FID, 'Note: Do not attempt to replicate Co-worker mean as I do not know definition. Avg size of entering firms is defined in terms of nprime, while avg size of exiting firms is defined in terms of n. \n');
% fprintf(FID, '}} \\end{minipage}');
fclose(FID);

%% Table 2

% One could code this directly easily enough since there are closed form
% solutions for kbar and nbar in terms of (s,tau). But will use the VFI
% Toolkit commands to show how to apply them.

FnsToEvaluate.kbar = @(aprime,a,z1,z2,agentmass,alpha,gamma,r,w,taurate,subsidyrate) (alpha/r)^((1-gamma)/(1-gamma-alpha)) *(gamma/w)^(gamma/(1-gamma-alpha)) *(z1*(1-((z2>=0)*taurate+(z2<0)*subsidyrate)*z2))^(1/(1-alpha-gamma));
FnsToEvaluate.nbar = @(aprime,a,z1,z2,agentmass,alpha,gamma,r,w,taurate,subsidyrate) ((1-((z2>=0)*taurate+(z2<0)*subsidyrate)*z2)*z1)^(1/(1-alpha-gamma)) *(alpha/r)^(alpha/(1-gamma-alpha)) *(gamma/w)^((1-alpha)/(1-gamma-alpha)); % which evaluates to Nbar in the aggregate
FnsToEvaluate.output = @(aprime,a,z1,z2,agentmass, alpha,gamma,r,w,taurate,subsidyrate) ((1-((z2>=0)*taurate+(z2<0)*subsidyrate)*z2))^((alpha+gamma)/(1-gamma-alpha))*z1^(1/(1-gamma-alpha)) *(alpha/r)^(alpha/(1-gamma-alpha)) *(gamma/w)^(gamma/(1-gamma-alpha));
% To use agentmass as input you must use that exact name, and it must be first input after z

ValuesOnGrid=EvalFnOnAgentDist_ValuesOnGrid_Case1_Mass(StationaryDist.mass, Policy, FnsToEvaluate, Params, [],EntryExitParamNames, n_d, n_a, n_z, [], a_grid, z_grid, [], simoptions);

ProbDensityFns=EvalFnOnAgentDist_pdf_Case1(StationaryDist, Policy, FnsToEvaluate, Params, [], n_d, n_a, n_z, [], a_grid, z_grid, [], simoptions, EntryExitParamNames);

% s_grid.^(1/(1-Params.gamma-Params.alpha))
nbarValues=ValuesOnGrid.nbar(:,:,:);
nbarValues=ValuesOnGrid.nbar(:,:,:);
normalize_employment=nbarValues(1,1,2); % Normalize so that smallest occouring value of nbar in the baseline is equal to 1.
nbarValues=nbarValues./normalize_employment;

Partion1Indicator=logical(nbarValues<5);
Partion2Indicator=logical((nbarValues>=5).*(nbarValues<50));
Partion3Indicator=logical(nbarValues>=50);

% Check that the following is equal to prod(n_z), so 300
sum(sum(Partion1Indicator+Partion2Indicator+Partion3Indicator))

ShareOfEstablishments(1)=sum(sum(StationaryDist.pdf(Partion1Indicator)));
ShareOfEstablishments(2)=sum(sum(StationaryDist.pdf(Partion2Indicator)));
ShareOfEstablishments(3)=sum(sum(StationaryDist.pdf(Partion3Indicator)));

Output_pdf=ProbDensityFns.output;
ShareOfOutput(1)=sum(sum(sum(Output_pdf(Partion1Indicator))));
ShareOfOutput(2)=sum(sum(sum(Output_pdf(Partion2Indicator))));
ShareOfOutput(3)=sum(sum(sum(Output_pdf(Partion3Indicator))));

Labour_pdf=ProbDensityFns.nbar;
ShareOfLabour(1)=sum(sum(sum(Labour_pdf(Partion1Indicator))));
ShareOfLabour(2)=sum(sum(sum(Labour_pdf(Partion2Indicator))));
ShareOfLabour(3)=sum(sum(sum(Labour_pdf(Partion3Indicator))));

Capital_pdf=ProbDensityFns.kbar;
ShareOfCapital(1)=sum(sum(sum(Capital_pdf(Partion1Indicator))));
ShareOfCapital(2)=sum(sum(sum(Capital_pdf(Partion2Indicator))));
ShareOfCapital(3)=sum(sum(sum(Capital_pdf(Partion3Indicator))));

AverageEmployment(1)=sum(sum(nbarValues(Partion1Indicator).*StationaryDist.pdf(Partion1Indicator)))/sum(sum(StationaryDist.pdf(Partion1Indicator)));
AverageEmployment(2)=sum(sum(nbarValues(Partion2Indicator).*StationaryDist.pdf(Partion2Indicator)))/sum(sum(StationaryDist.pdf(Partion2Indicator)));
AverageEmployment(3)=sum(sum(nbarValues(Partion3Indicator).*StationaryDist.pdf(Partion3Indicator)))/sum(sum(StationaryDist.pdf(Partion3Indicator)));

%Table 2
FID = fopen('./RestucciaRogerson2008_Table2.tex', 'w');
fprintf(FID, 'Distribution statistics of benchmark economy \\\\ \n');
fprintf(FID, '\\begin{tabular*}{1.00\\textwidth}{@{\\extracolsep{\\fill}}llcr} \n \\hline \\hline \n');
fprintf(FID, ' & \\multicolumn{3}{l}{Establishment size (number of employees)} \\\\ \\cline{2-4} \n');
fprintf(FID, ' & $<$5 & 5 to 49 & $\\geq$50 \\\\ \\hline \n');
fprintf(FID, 'Share of establishments & %8.2f & %8.2f & %8.2f \\\\ \n', ShareOfEstablishments);
fprintf(FID, 'Share of output         & %8.2f & %8.2f & %8.2f \\\\ \n', ShareOfOutput);
fprintf(FID, 'Share of labour         & %8.2f & %8.2f & %8.2f \\\\ \n', ShareOfLabour);
fprintf(FID, 'Share of capital        & %8.2f & %8.2f & %8.2f \\\\ \n', ShareOfCapital);
fprintf(FID, 'Share of employment     & %8.2f & %8.2f & %8.2f \\\\ \n', AverageEmployment);
fprintf(FID, '\\hline \n \\end{tabular*} \n');
% fprintf(FID, '\\begin{minipage}[t]{1.00\\textwidth}{\\baselineskip=.5\\baselineskip \\vspace{.3cm} \\footnotesize{ \n');
% fprintf(FID, 'Note: Do not attempt to replicate Co-worker mean as I do not know definition. Avg size of entering firms is defined in terms of nprime, while avg size of exiting firms is defined in terms of n. \n');
% fprintf(FID, '}} \\end{minipage}');
fclose(FID);

%% Calculate a bunch of things related to those reported in Table 3 just to show how.
FnsToEvaluate.subsidy = @(aprime,a,z1,z2,agentmass, alpha,gamma,r,w,taurate,subsidyrate) (z2<0)*subsidyrate* ((1-((z2>=0)*taurate+(z2<0)*subsidyrate)*z2))^((alpha+gamma)/(1-gamma-alpha))*z1^(1/(1-gamma-alpha)) *(alpha/r)^(alpha/(1-gamma-alpha)) *(gamma/w)^(gamma/(1-gamma-alpha)); % (z2<0)*subsidyrate)* output
% Following is just 'indicator for subsidised' times output, needed to calculate Ys
FnsToEvaluate.outputofsubsidised = @(aprime,a,z1,z2,agentmass, alpha,gamma,r,w,taurate,subsidyrate) (z2<0)*((1-((z2>=0)*taurate+(z2<0)*subsidyrate)*z2))^((alpha+gamma)/(1-gamma-alpha))*z1^(1/(1-gamma-alpha)) *(alpha/r)^(alpha/(1-gamma-alpha)) *(gamma/w)^(gamma/(1-gamma-alpha));

AggValues=EvalFnOnAgentDist_AggVars_Case1(StationaryDist, Policy, FnsToEvaluate, Params, [], n_d, n_a, n_z, d_grid, a_grid, z_grid, [], simoptions, EntryExitParamNames);

Output.Y=AggValues.output.Aggregate;
Output.N=AggValues.nbar.Aggregate;
Output.K=AggValues.kbar.Aggregate;
Output.KdivY=Output.K/Output.Y;
% Params.w
Output.mass=StationaryDist.mass; % E in notation of RR2008
Output.TFP=(Output.Y/Output.N)./((Output.K/Output.N)^Params.alpha); % RR2008 call this 'A'.
Output.Ys_divY=AggValues.outputofsubsidised.Aggregate/AggValues.output.Aggregate; % The variable Ys/Y represents the output share of establishments that are receiving a subsidy
Output.SdivY=AggValues.subsidy.Aggregate/AggValues.output.Aggregate; % The variable S /Y is the total subsidies paid out to establishments receiving subsidies as a fraction of output
Output.tau_s=Params.subsidyrate; % The variable tau_s is the size of the subsidy required to generate a steady-state capital stock equal to that in the distortion-free economy




