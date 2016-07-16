%%% Example using a variant of the Basic RBC model (following Aruoba, Fernandez-Villaverde, & Rubio-Ramirez, 2006)

% Aruoba, Fernandez-Villaverde, & Rubio-Ramirez (2006) discuss solving this model with a wide variety of different numerical solution methods.

disp('Running BasicRealBusinessCycleModel_NumericalErrors.m')

% WARNING: Marcet-DenHaan statistics are not correctly implemented
disp('WARNING: Marcet-DenHaan statistics are not correctly implemented')

% Run code that solves the Basic RBC model.
BasicRealBusinessCycleModel

%% If you are just interested in how to solve the value function problem for the Basic Real Business Cycle Model you can stop here.
% The rest of the code reproduces some relevant parts of the Figures and Tables of Aruoba, Fernanadez-Villaverde, & Rubio-Ramirez (2006).

% The following four are used later as defining the areas relevant for graphs (and for Euler eqn errors)
[~,point7K_ss]=max(a_grid>0.7*K_ss);
[~,onepoint3K_ss]=max(a_grid>1.3*K_ss);
[~,zlow]=max(z_grid>-0.065);
[~,zhigh]=max(z_grid>0.065);

%% The Policy Functions
AFVRR2006_figs=figure;

% Figure 1 (Fig 3 if UseAlternativeParams=1)
subplot(3,2,1); plot(a_grid,d_grid(Policy(1,:,ceil(n_z/2))))
xlim([a_grid(point7K_ss),a_grid(onepoint3K_ss)]);

% Figure 2 (Fig 4 if UseAlternativeParams=1)
subplot(3,2,2); plot(a_grid,a_grid(Policy(2,:,ceil(n_z/2)))-(1-delta)*a_grid)
xlim([a_grid(point7K_ss),a_grid(onepoint3K_ss)]);


%% Simulated Densities of Output, Capital, and Consumption
% 1000 simulations of 500 points (I use a 100 point burn in, then do not say if they burn-in or just start in steady-state?)
NSims=1000;
simoptions.burnin=100;
simoptions.simperiods=500;
simoptions.parallel=2;
HistBins=500;

StationaryDist=zeros(n_a,n_z,NSims,'gpuArray');
for ii=1:NSims
    StationaryDist(:,:,ii)=StationaryDist_Case1_Simulation(Policy,n_d,n_a,n_z,pi_z,simoptions);
end

%Judging from the y-axes of the Figures it appears the 'densities' are calculated as histograms formed by summing across all of the simulations.
SteadyStateHist=sum(StationaryDist,3)*simoptions.simperiods; %Multiply by SimPeriods as otherwise they are the actual densities being summed
CapitalSteadyStateHist=zeros(n_a,1,'gpuArray');
OutputSteadyStateHist=zeros(HistBins,1,'gpuArray');
ConsumptionSteadyStateHist=zeros(HistBins,1,'gpuArray');
OutputSteadyStateValues=zeros(n_a,n_z,'gpuArray'); %Create matrix of the value of output for each point of domain of dist
ConsumptionSteadyStateValues=zeros(n_a,n_z,'gpuArray');
for a_c=1:n_a
    for z_c=1:n_z
        CapitalSteadyStateHist(Policy(2,a_c,z_c))=CapitalSteadyStateHist(Policy(2,a_c,z_c))+SteadyStateHist(a_c,z_c);
        OutputSteadyStateValues(a_c,z_c)=exp(z_grid(z_c))*(a_grid(a_c)^alpha)*(d_grid(Policy(1,a_c,z_c))^(1-alpha));
        ConsumptionSteadyStateValues(a_c,z_c)=exp(z_grid(z_c))*(a_grid(a_c)^alpha)*(d_grid(Policy(1,a_c,z_c))^(1-alpha))+(1-delta)*a_grid(a_c)-a_grid(Policy(2,a_c,z_c));
    end
end
minOutput=min(min(OutputSteadyStateValues));
maxOutput=max(max(OutputSteadyStateValues));
minConsumption=min(min(ConsumptionSteadyStateValues));
maxConsumption=max(max(ConsumptionSteadyStateValues));
OutputGrid=linspace(minOutput,maxOutput+(1/HistBins)*(maxOutput-minOutput),HistBins);
ConsumptionGrid=linspace(minConsumption,maxConsumption+(1/HistBins)*(maxConsumption-minConsumption),HistBins);
for a_c=1:n_a
    for z_c=1:n_z
        FoundBin=0;
        counter=1;
        while FoundBin==0
            if counter==HistBins
                OutputSteadyStateHist(counter)=OutputSteadyStateHist(counter)+SteadyStateHist(a_c,z_c);
                FoundBin=1;
            elseif OutputGrid(counter)<=OutputSteadyStateValues(a_c,z_c) && OutputSteadyStateValues(a_c,z_c)<OutputGrid(counter+1)
                OutputSteadyStateHist(counter)=OutputSteadyStateHist(counter)+SteadyStateHist(a_c,z_c);
                FoundBin=1;
            else
                counter=counter+1;
            end
        end
        FoundBin=0;
        counter=1;
        while FoundBin==0
            if counter==HistBins
                ConsumptionSteadyStateHist(counter)=ConsumptionSteadyStateHist(counter)+SteadyStateHist(a_c,z_c);
                FoundBin=1;
            elseif ConsumptionGrid(counter)<=ConsumptionSteadyStateValues(a_c,z_c) && ConsumptionSteadyStateValues(a_c,z_c)<ConsumptionGrid(counter+1)
                ConsumptionSteadyStateHist(counter)=ConsumptionSteadyStateHist(counter)+SteadyStateHist(a_c,z_c);
                FoundBin=1;
            else
                counter=counter+1;
            end
        end
    end
end

% Fig 5, when UseAlternativeParams=1
subplot(3,2,3); plot(OutputGrid, OutputSteadyStateHist)

% Fig 6, when UseAlternativeParams=1
subplot(3,2,4); plot(a_grid,CapitalSteadyStateHist)

% Fig 7, when UseAlternativeParams=1
subplot(3,2,5); plot(ConsumptionGrid, ConsumptionSteadyStateHist)


%% Calculate the DenHaan-Marcet Chi-squared statistic
% y_t is the vector (c_t;z_t)
% Write model as f(y_t)=E_t[phi(y_{t+1},y_{t+2},...)]
% So for this model
% fy_t is the vector (MU_c; z_t)
% Et_phi is a vector ( ; rho*z_{t-1}+e)

%NOTE: I use different simulations for the chi-squared (here the time
%dimension is necessary, previously it was not).
% First we generate a time series of indexes for the a & z variables (of
% size (periods,num_a_vars+num_z_vars))
simoptions.seedpoint=[ceil(n_a/2),ceil(n_z/2)];
simoptions.polindorval=PolIndOrVal;
TimeSeriesIndexesSims=zeros(2,500,NSims);
for ii=1:NSims
    TimeSeriesIndexesSims(:,:,ii)=SimTimeSeriesIndexes_Case1(Policy,n_d,n_a,n_z,pi_z,simoptions);
end

% %Define the functions which we wish to create time series for (from the TimeSeriesIndexes)
% TimeSeriesFn_1 = @(d_val,d_ind,aprime_val,aprime_ind,a_val,a_ind,z_val,z_ind) a_val; %Capital Stock
% TimeSeriesFn_2 = @(d_val,d_ind,aprime_val,aprime_ind,a_val,a_ind,z_val,z_ind) z_val; %Investment
% TimeSeriesFn={TimeSeriesFn_1, TimeSeriesFn_2};
% TimeSeries=TimeSeries_Case1(TimeSeriesIndexes,Policy, TimeSeriesFn, n_d, n_a, n_z, d_grid, a_grid, z_grid);

m=2; q=1;
DenHaanMarcetStat=zeros(NSims,1);
B_T=zeros(m*q,NSims);
A_T=zeros(m*q,m*q,NSims);
for ii=1:NSims
    TimeSeriesIndexes=TimeSeriesIndexesSims(:,:,ii);


    
    fy_t=zeros(m,SimPeriods-2);
    E_t_phi=zeros(m,SimPeriods-2);
    u_tplus1=zeros(m,SimPeriods-2);
    
    hx_t=ones(q,SimPeriods-2); % q-dimensional

    B_t=zeros(m*q,SimPeriods-2);
    A_t=zeros(m*q,m*q,SimPeriods-2);
    
    for t=1:SimPeriods-2
        a_t_c=TimeSeriesIndexes(1,t);
        a_t=a_grid(a_t_c); %=TimeSeries(t,1)
        z_t_c=TimeSeriesIndexes(2,t);
        z_t=z_grid(z_t_c);
        a_tplus1_c=TimeSeriesIndexes(1,t+1); %=Policy(2,a_t_c,z_t_c);
        a_tplus1=a_grid(a_tplus1_c);
        l_t=d_grid(Policy(1,a_t_c,z_t_c));
        c_t=exp(z_t)*(a_t^alpha)*(l_t^(1-alpha))+(1-delta)*a_t-a_tplus1;
        
        z_tplus1_c=TimeSeriesIndexes(2,t+1);
        z_tplus1=z_grid(z_tplus1_c);
        a_tplus2_c=TimeSeriesIndexes(1,t+2); %=Policy(2,a_tplus1_c,zprime_c);
        a_tplus2=a_grid(a_tplus2_c);
        l_tplus1=d_grid(Policy(1,a_tplus1_c,z_tplus1_c));
        c_tplus1=exp(z_tplus1)*(a_tplus1^alpha)*(l_tplus1^(1-alpha))+(1-delta)*a_tplus1-a_tplus2;
        
        MPK_t=alpha*exp(z_t)*(a_t^(alpha-1))*(l_t^(1-alpha));
        
        fy_t(:,t)=[theta*(c_t^(theta*(1-tau)-1))*((1-l_t)^((1-theta)*(1-tau))); rho*z_t];
        E_t_phi(:,t)=[beta*(c_tplus1^(theta*(1-tau)-1))*((1-l_tplus1)^((1-theta)*(1-tau)))*(1+MPK_t-delta); z_tplus1];
        u_tplus1(:,t)=fy_t(:,t)-E_t_phi(:,t); % m-dimensional
        
        %hx_t(:,t)
        
        B_t(:,t)=kron(u_tplus1(:,t),hx_t(:,t));
        A_t(:,:,t)=kron(u_tplus1(:,t),hx_t(:,t))*kron(u_tplus1(:,t),hx_t(:,t))';

    end
    
    
    %B
    B_T(:,ii)=(1/(SimPeriods-2))*sum(B_t,2);
    %A
    A_T(:,:,ii)=(1/(SimPeriods-2))*sum(A_t,3);
    
    DenHaanMarcetStat(ii)=(SimPeriods-2)*B_T(:,ii)'*(A_T(:,:,ii)^(-1))*B_T(:,ii);
    
end

%Having calculated the DenHaan-Marcet stats, we now want to look at the 
Chi2cdfvalues=chi2cdf(DenHaanMarcetStat/(SimPeriods-2),1);
MoreThan5Percent=sum(Chi2cdfvalues>0.95);
LessThan5Percent=sum(Chi2cdfvalues<0.05);

% Table 3 (Table 4 if UseAlternativeParams=1)
fprintf('DenHann-Marcet chi-squared test \n')
fprintf('            Less than 5%%    More than 5%% \n')
fprintf('Value fn    %8.2f           %8.2f   \n', [LessThan5Percent, MoreThan5Percent])


%% Calculate the Euler eqn errors

EulerEqnErrors=zeros(n_a,n_z,'gpuArray');
u_c_tplus1=zeros(1,n_z,'gpuArray');
R=zeros(1,n_z,'gpuArray');
for a_c=1:n_a
    a_t=a_grid(a_c);
    for z_c=1:n_z
        a_tplus1_c=Policy(2,a_c,z_c);
        a_tplus1=a_grid(a_tplus1_c);
        l_t=d_grid(Policy(1,a_c,z_c));
        c_t=exp(z_grid(z_c))*(a_t^alpha)*(l_t^(1-alpha))+(1-delta)*a_t-a_tplus1;
        
        for zprime_c=1:n_z
            a_tplus2_c=Policy(2,a_tplus1_c,zprime_c);
            a_tplus2=a_grid(a_tplus2_c);
            l_tplus1=d_grid(Policy(1,a_tplus1_c,zprime_c));
            c_tplus1=exp(z_grid(zprime_c))*(a_tplus1^alpha)*(l_tplus1^(1-alpha))+(1-delta)*a_tplus1-a_tplus2;
            u_c_tplus1(zprime_c)=theta*(c_tplus1^(theta*(1-tau)-1))*((1-l_tplus1)^((1-theta)*(1-tau)));
            R(zprime_c)=1+alpha*exp(z_grid(zprime_c))*(a_tplus1^(alpha-1))*(l_tplus1^(1-alpha))-delta;
        end
        E_t=(u_c_tplus1.*R)*pi_z(z_c,:)';
        EulerEqnErrors(a_c,z_c)=1-((beta*E_t/(theta*(1-l_t)^((1-theta)*(1-tau))))^(1/(theta*(1-tau)-1)))/c_t;
    end
end


% Figure 8/9 (Fig 10/11 if UseAlternativeParams=1)
subplot(3,2,6); plot(a_grid,log10(abs(EulerEqnErrors(:,ceil(n_z/2)))))
ylim([-10,-1]); xlim([a_grid(point7K_ss),a_grid(onepoint3K_ss)]);

% Table 5 (Table 6 if UseAlternativeParams=1)
AbsMaxEE=max(max(log10(abs(EulerEqnErrors(point7K_ss:onepoint3K_ss,zlow:zhigh)))));
IntegEE=log10(abs(sum(sum(EulerEqnErrors.*(SteadyStateHist/(NSims*SimPeriods))))));
fprintf('Table 5: Euler errors (Abs(log10)) \n')
fprintf('            Absolute max Euler Error     Integral of the Euler Errors \n')
fprintf('Value fn:   %8.2f                        %8.6f \n', [AbsMaxEE, IntegEE])



%%
print(AFVRR2006_figs,'-djpeg', './SavedOutput/AFVRR2006_figs.jpg')




