%%% Example using a variant of the Basic RBC model (following Aruoba, Fernandez-Villaverde, & Rubio-Ramirez, 2006)

disp('Running BasicRealBusinessCycleModel_BusinessCycleStatistics.m')

% Run code that solves the Basic RBC model.
BasicRealBusinessCycleModel

%% We will generate some simulated data and then use this to calculate the standard business cycle statistics
% First we generate a time series of indexes for the d,a & z variables (of
% size (num_d_vars+num_a_vars+num_a_vars+num_z_vars,periods)) (the first
% num_a_vars is for aprime).
simoptions.parallel=2;
simoptions.seedpoint=[ceil(n_a/2),ceil(n_z/2)];
simoptions.burnin=1000;
simoptions.simperiods=1000;
TimeSeriesIndexes=SimTimeSeriesIndexes_Case1(Policy,n_d,n_a,n_z,pi_z,simoptions);

%Define the functions which we wish to create time series for (from the TimeSeriesIndexes)
TimeSeriesFn_1 = @(d_val,d_ind,aprime_val,aprime_ind,a_val,a_ind,z_val,z_ind) a_val; %Capital Stock
TimeSeriesFn_2 = @(d_val,d_ind,aprime_val,aprime_ind,a_val,a_ind,z_val,z_ind) aprime_val-(1-delta)*a_val; %Investment
TimeSeriesFn={TimeSeriesFn_1, TimeSeriesFn_2};

TimeSeries=TimeSeries_Case1(TimeSeriesIndexes,Policy, TimeSeriesFn, n_d, n_a, n_z, d_grid, a_grid, z_grid,simoptions);

StdBusCycleStats=[mean(TimeSeries(1,:)),mean(TimeSeries(2,:)); var(TimeSeries(1,:)),var(TimeSeries(2,:))]



