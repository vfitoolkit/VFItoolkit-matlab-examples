% Example using a variant of the Basic RBC model (following Aruoba, Fernandez-Villaverde, & Rubio-Ramirez, 2006)

disp('Running BasicRealBusinessCycleModel_BusinessCycleStatistics.m')

% Run code that solves the Basic RBC model.
BasicRealBusinessCycleModel

%% We will generate some simulated data and then use this to calculate the standard business cycle statistics
% Set some options, the following are actually just the defaults anyway
simoptions.burnin=1000;
simoptions.simperiods=10000;

% Define the functions which we wish to create time series for (from the TimeSeriesIndexes)
FnsToEvaluate.K = @(d,aprime,a,z) a; %Capital Stock
FnsToEvaluate.I = @(d,aprime,a,z,delta) aprime-(1-delta)*a; %Investment
% Note that the inputs are the states (in order) followed by any other parameter values

TimeSeries=TimeSeries_Case1(Policy, FnsToEvaluate, Params, n_d, n_a, n_z, d_grid, a_grid, z_grid,pi_z,simoptions);

StdBusCycleStats=[mean(TimeSeries.K),mean(TimeSeries.I); var(TimeSeries.K),var(TimeSeries.I)]

%% We could look at the simulated time series directly
figure(2)
plot(TimeSeries.K)
title('Simulated time series for aggregate physical capital (K)')
