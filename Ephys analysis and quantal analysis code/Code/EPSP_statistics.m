function [EPSP_stats, EPSP_ProbDist, noise_stats, noise_ProbDist] = EPSP_statistics(loadpathEPSPdata, stable_period)
%% Extracts the mean, SD, skewness, and number of stable sweeps from the EPSPs and the NOISE and brings them to the main script
%  mean, SD, skewness, n stable sweeps, median are written in table (EPSP_stats, noise stats)
%  additionally, a fitdist object is created from the EPSPs and noise distributions (EPSP_ProbDist, noise_ProbDist)

%----------------------------------------------------------------------------------------------------------------------------------
% IMPORTANT: when the sample (EPSPs, noise) is drawn from a distribution (as in our case) matlab implements a 'bias-correction',
% that needs to be activated by setting flag=0 in skewness function:
% y = skewness(X,flag)
% X = sample
% flag = 0 -> skewness is corrected for finite sample size (bias depends on the sample size)
% THIS CORRECTION IS STATE OF THE ART, AND IS IMPLEMENTED IN THE SKEWNESS ALGORITHMS OF OTHER STATISTICS PROGRAMS AS WELL!
% it is activated in all skewness functions around this pipeline
%----------------------------------------------------------------------------------------------------------------------------------


%% load data
data = xlsread(loadpathEPSPdata);
EPSP = data(:,2);
noise = data(:,3);
sweeps = (1 : length(EPSP))';
    
%% write EPSP statistics for the stable period of rec. into new arrays
stable_sweeps = sweeps(stable_period(:,1) == 1);
stable_EPSP = EPSP(stable_period(:,1) == 1);
stable_noise = noise(stable_period(:,1) == 1);

mu  = mean(stable_EPSP);
sigma = std(stable_EPSP);
gamma = skewness(stable_EPSP,0);    % IMPORTANT: 0 corrects for finite sampling of distribution
n_sweeps = length(stable_sweeps);
med = median(stable_EPSP);

EPSP_stats = table(mu,sigma,gamma,n_sweeps, med);

%% same for NOISE
mu_noise  = mean(stable_noise);
sigma_noise = std(stable_noise);
gamma_noise = skewness(stable_noise,0);         % IMPORTANT: 0 corrects for finite sampling of distribution
med_noise = median(stable_noise);

noise_stats = table(mu_noise,sigma_noise,gamma_noise,n_sweeps, med_noise);

%% fit a normal distribution object to the stable EPSP histogram
EPSP_ProbDist = fitdist(stable_EPSP,'normal');
noise_ProbDist = fitdist(stable_noise,'normal');
end