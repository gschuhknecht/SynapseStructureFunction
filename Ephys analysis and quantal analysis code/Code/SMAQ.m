function[quantal_parameters, ConfidenceIntervalsCentralMoments] = SMAQ(savepath, filename, loadpathEPSPdata, EPSP_distribution_stats, stable_period, EPSP_ProbDist, save,filled)
%% Implementation of Statistical Moments Analysis of Quanta (SMAQ)
% Finds the analytical SMAQ solution and writes it in quantal_parameters 
% Generates the plot with the graphical solution of SMAQ for the correct range of P, data saved in parameter_space
% The confidence intervals associated with the empiral measuremtents of mean and SD are computed from bootstrapping and by fitting a Gaussian
% (don't use those), CI for skewness is computed from the SES equation; all these are exported in Confidence_Intervals
% Plots the analytical error bounds associated with the equations for N, Ngamma, and Q (using bootstrapping CIs for mean and SD, since they come from non-normal distribution)

%% load data and get stable EPSPs
% we need the data for finding empirical confidence intervals of mean and SD using bootstrapping
data = xlsread(loadpathEPSPdata);
EPSP = data(:,2);
stable_EPSP = EPSP(stable_period(:,1) == 1);

%% retrieve statistics of EPSP distribution
%  retrieved from the EPSP_distribution_stats table, which contains empirical measures of mean, SD, skewness, number of sweeps,  median
EPSP_distribution_stats = table2array(EPSP_distribution_stats);
mu = EPSP_distribution_stats(1,1);
sigma = EPSP_distribution_stats(1,2);
gamma = EPSP_distribution_stats(1,3);
n_sweeps = EPSP_distribution_stats(1,4);

%% determine the confidence intervals
% This equation calculates the CI for skewness, based on n_sweeps
% Documentation of GraphPadPrism: https://www.graphpad.com/support/faqid/1577/
SES = sqrt((6*n_sweeps*(n_sweeps-1)) / ((n_sweeps-2)*(n_sweeps+1)*(n_sweeps+3)));   % Standard Error of Skewness
CI_gam = SES * 1.96;                                                                % 'Margin of Error of Skewness'
        % => (The CI of skewness = skewness +- the margin of error (-> the margin goes both up and down, see below)

% the EPSP data is not normally-distributed, therefore the CIs for mean and SD are defined empirically by bootstrapping:
CI_mean_BS = bootci(10000, @mean, stable_EPSP); % default alpha is 0.05 -> 95% CI
CI_std_BS = bootci(10000, @std, stable_EPSP); % default alpha is 0.05 -> 95% CI

CI_mu_low = CI_mean_BS(1);
CI_mu_high = CI_mean_BS(2);
CI_sigma_low = CI_std_BS(1);
CI_sigma_high = CI_std_BS(2);

CI_limits_BS = [CI_mean_BS, CI_std_BS];

    %----------------------------------------------------------------------------------------------------------------------------%
    % this code extracts the CIs for mean and SD from the probdist object, but assuming a normal distribution
    % only done to save these CIs into the table as a comparison with the CIs derived from BootStrapping, these CIs are not used for calculations
          CI_limits_NormalDist = paramci(EPSP_ProbDist,'Alpha',0.05);    % for 95% CI, alpha is 0.05, can be adjusted
        % CI_mu_low = CI_limits(1,1);         % CI of mu
        % CI_mu_high = CI_limits(2,1);
        % 
        % CI_sigma_low = CI_limits(1,2);      % CI of mean
        % CI_sigma_high = CI_limits(2,2);
    %----------------------------------------------------------------------------------------------------------------------------%

%% write CI into table for exporting to stability function saving in Excel
ConfidenceIntervalsCentralMoments = array2table([CI_limits_BS, CI_limits_NormalDist, [(gamma + CI_gam); (gamma - CI_gam)]],'VariableNames',{'CI_mu_bootstrap','CI_sigma_bootstrap','CI_mu_fitGaussian','CI_sigma_fitGaussian','CI_gamma'})

%% Derive the full analytical solution
N = - ((mu^2) / ((gamma * mu - sigma)*sigma));
P = (gamma*mu - sigma) / (gamma*mu - 2*sigma);
Q = ( (-gamma * mu * sigma + (2 * sigma^2)) / mu);

quantal_parameters = table(N,P,Q)