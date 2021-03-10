%% ANALYSISmainfile
% this is the main file from which all data analysis is called


%% ======================================================================================================================================================== %%
%% SETTING PARAMETERS
%--------------------%
%% I   Define location of EPSP xls files and the storage location
savepath = '/Users/.../'; % wherever you want to store your results (your .eps plots)
filepath = '/Users/.../Raw Data/'; %%% the path where your EPSP recoding files are stored
filename = 'Experiment1_EPSP_measurements.xls';     %the EPSP recording (.xls or .xlsx) file that you want to analyze
loadpathEPSPdata = [filepath,filename];


%% II   SAVING and figure appearance
% save data and figures?
save = 0;                % if 0 -> don't save
filled = 0;              % not filled is compatible with Illustrator

% make new folder if saving is enabled
if save == 0 
    mkdir(savepath);
end


%% III  Set parameters for stability analyses
% -> these are the 'optimal stability criteria'
blocksize = 25           % size of the blocks that are averaged
mean_stable = 3;         % how many SEM is the mean allowed to fluctuate (used by hardingham, 2010)
SD_stable = 0.30;        % how many % is the SD allowed to fluctuate     (used by hardingham, 2006)

min_initial_stable = 4;  % how many blocks need to be stable including the 1st sweeps? -> 100 consecutive sweeps at beginning of recording!


%% IV   Parameters for MonteCarlo Simulations
% when you run the MonteCarlo Simulations, which range of parameters do you wish so simulate, how many runs and what range of results should be plotted?
n_runs = 1000;                           % idealy 10k

% range of quantal parameters to be simulated by MonteCarlo
N_range = [1 : 20];                    
P_range = [0.1 : 0.1 : 0.9];           
Q_range = [0.1 : 0.1 : 1.5];            
QuParam_ranges = {N_range, P_range, Q_range};

% range of quantal parameters to be plotted in the figure
N_plotting_range = [1 : 16];            % which range of N would you like to plot?
P_plotting_range = [0.1 : 0.1 : 0.9];   % which range of P would you like to plot?
Q_plotting_range = [0.1 : 0.1 : 1];     % which range of Q would you like to plot?
QuParam_Plotting_ranges = {N_plotting_range, P_plotting_range, Q_plotting_range};


%% ======================================================================================================================================================== %%
%% RUNNING THE ANALYSES
%-----------------------%

%% I     Stability Protocol
    % if at least 100 consecutive sweeps are stable with respect to the 1st bin, you take those, 
    % if not, the longest consecutive stable epoch is chosen
    
[SSN_stable_EPSPS, SSN_EPSP_STATS] = Stability_analysis(loadpathEPSPdata, savepath, filename, blocksize,mean_stable,SD_stable,min_initial_stable, filled, save);
        % SSN_stable_EPSPS = logical array of length all sweeps of experiment, showing the selected stable sweeps
        % SSN_EPSP_STATS   = block center, mean, SD, SEM of the EPSP that fall into that block

%% II    Plot histograms of stable sweeps
    binwidth = [0.02, 0.05, 0.075, 0.1, 0.15, 0.2];         % which binwidth would you like to plot?
    movmean_winsize = [7 5 5 3 3 3];                        % n bins for moving average corresponding to binwidth of histograms
Plot_stable_histogram(loadpathEPSPdata, savepath, filename, binwidth, movmean_winsize, SSN_stable_EPSPS, 1); 


%% III   Extracting the statistics from the selected stable sweeps
[EPSP_distribution_stats, EPSP_ProbDist, noise_distribution_stats, noise_ProbDist] = EPSP_statistics(loadpathEPSPdata, SSN_stable_EPSPS);
        % EPSP_distribution_stats = mean, SD, skewness, n sweeps of the stable EPSP distribution
        % EPSP_ProbDist           = stable EPSPs turned into ProbDist object, used later to extract mean and SD in the NPQ_approximation
        % same for noise
        
%% IV    Run quantal analysis using SMAQ
[quantal_parameters, ConfidenceIntervalsCentralMoments] = SMAQ(savepath, filename, loadpathEPSPdata, EPSP_distribution_stats, SSN_stable_EPSPS, EPSP_ProbDist, save,filled);
        %quantal_parameters  = unique solution for N, P, Q found by quantal approx method
        %ConfidenceIntervalsCentralMoments = analytical estimate of lower and upper bound for mean, SD, skewness giving 95% confidence (bootstrapping)

%% V     Save stability data and quantal approximation results
Stability_Criteria = table(blocksize, mean_stable, SD_stable);
if save == 1 
    Save_file(loadpathEPSPdata, savepath, filename, SSN_stable_EPSPS, EPSP_distribution_stats, noise_distribution_stats, ConfidenceIntervalsCentralMoments, quantal_parameters, Stability_Criteria);
end
   
%% VI    Assess the reliability of SMAQ with Bayesian-inspired Monte Carlo simulations, as reported in the paper
% Simulating the entire range of quantal parameters to find the range of N,
% P, Q that could have also produced our SMAQ results:
[statistics_All_Model_MC_simulations,...
    All_Model_MC_simulations,...
    histograms_database_N, histograms_database_P, histograms_database_Q,...
    histograms_database_models,...
    QuParam_ranges] = ...
        MonteCarloSim_SMAQ_Reliability_main_file(savepath, filename, QuParam_ranges, noise_distribution_stats, n_runs, save);
     
MonteCarloSim_SMAQ_Reliability_binning_analysis(statistics_All_Model_MC_simulations,...
    All_Model_MC_simulations,...
    histograms_database_N, histograms_database_P, histograms_database_Q,...
    histograms_database_models,...
    savepath, filename, QuParam_ranges, QuParam_Plotting_ranges,quantal_parameters,n_runs, save);
