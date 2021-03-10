function[statistics_All_Model_MC_simulations,...
    All_Model_MC_simulations,...
    histograms_database_N, histograms_database_P, histograms_database_Q,...
    histograms_database_models,...
    QuParam_ranges] = MonteCarloSim_SMAQ_Reliability_main_file(savepath, filename, QuParam_ranges, noise_distribution_stats, n_runs, saving)
%% MonteCarloSim_QuantApproxReliability
% Monte Carlo Simulation that generates histograms from all permutations of a
% range of N, P, Q. Histograms are then analysed with SMAQ to recover the
% quantal parameters that gave rise to the histograms.
% Then the results of SMAQ are compared with the true quantal parameters behind the histograms to
% investigate the reliability of SMAQ.
%
% Statistical question: 'If SMAQ predcits - e.g. a certain N - , which
% binomial models ever produced this results, and with which probability?
% So given a certain N as result of SMAQ, which models could have produced
% this result with 95% probability?
%
% The simulation uses the same noise level and number of sweeps of the experimental
% histograms

% the final data matrices are stored in the folder that is open

%% Range of N, P, Q to be simulated, imported from main script
N_range = QuParam_ranges{1};
P_range = QuParam_ranges{2};
Q_range = QuParam_ranges{3};

%% what is the experimentally recorded noise level?
noise_distribution_stats = table2array(noise_distribution_stats);
noise = noise_distribution_stats(2);
n_sweeps = noise_distribution_stats(4);

Gmean = 0; % noise offset, should be zero

%% saving -> make directory
if saving == 1 
    mkdir(savepath,'MonteCarloSim_SMAQ_Reliability_mfiles');
    savename = [savepath,'MonteCarloSim_SMAQ_Reliability_mfiles/' filename '_'];
end

%% generate all the combinations of parameters:
test_models = combvec(N_range,P_range,Q_range); % lists the permutations of N,P,Q that will be tested
fprintf(' \n \n number of models to be tested: %d\n \n \n', size(test_models,2)); % display number of models to be tested

%% set up binning for the solutions provided by SMAQ:
% the results of SMAQ are binned according to:
edges_N = [0.5 : 1 : 50.5];     % all bins are same
edges_P = [0.05 : 0.1 : 0.95];  % all bins are same
edges_Q = [0.05 : 0.1 : 3.05];  % all bins are same
    
%% set up matrices that store results
histograms_database_models = NaN(size(test_models,2),3);            % stores the underlying binomial models

histograms_database_N = NaN(size(test_models,2),size(edges_N,2)-1); % stores binned solutions for N of Quantal Approximation Method
histograms_database_P = NaN(size(test_models,2),size(edges_P,2)-1); % stores binned solutions for P of Quantal Approximation Method
histograms_database_Q = NaN(size(test_models,2),size(edges_Q,2)-1); % stores binned solutions for Q of Quantal Approximation Method

All_Model_MC_simulations = zeros(n_runs,6,size(test_models,2)) ;    %initialize the final results matrix

%% run the simulation
% run through all the MODELS
for n = 1 : size(test_models,2)
    fprintf('\n run number %d: \n', n);    % outputs the number of the current model
    tic;
        
    % the quantal parameters that are modeled in this run are pulled from the test_models matrix
    N = test_models(1,n);
    P = test_models(2,n);
    Q = test_models(3,n);
    
    % call the function that implements the Monte-Carlo Simulation 
    MC_results = MonteCarlo_BinomialModel_SMAQ(N, P, Q, noise, Gmean, n_sweeps, n_runs);
        % simulation gives results for [mu sigma gamma Nsim Psim Qsim] and
        % stores them in an array of size (n_runs, 6) 
    
    % write the results of the simulation into the results matrix
    All_Model_MC_simulations(:,:,n) = MC_results;
        % the results matrices of the models are stacked 'behind' each other in a 3D matrix
        % this matrix will serve as the final look-up table

    histograms_database_models(n,:) = [N P Q];
        % write into database which model we just simulated, will also be a
        % look-up table later, to recover what the  underlying binomial model were
    
    % BINNING THE SOLUTIONS OF SMAQ
    histograms_database_N(n,:) = histcounts(All_Model_MC_simulations(:,4,n),edges_N);  % binned solutions for N of SMAQ
    histograms_database_P(n,:) = histcounts(All_Model_MC_simulations(:,5,n),edges_P);  % binned solutions for P of SMAQ
    histograms_database_Q(n,:) = histcounts(All_Model_MC_simulations(:,6,n),edges_Q);  % binned solutions for Q of SMAQ
        % The SMAQ solutions are binned now, so solutions outside the bins specified above will be neglected...
        
    toc;    
end


%% save statistics of SMAQ predictions for N,P,Q over all simulations
statistics_All_Model_MC_simulations = zeros(3,5,size(test_models,2)) ;
% enter simulated quantal parameters of model in 1st column
statistics_All_Model_MC_simulations(1:3, 1,:) = test_models(:,:);
% enter means of SMAQ solutions for N, P, Q in 2nd column
statistics_All_Model_MC_simulations(1:3, 2,:) = mean(All_Model_MC_simulations(:,4 : 6 ,:));
% enter medians of SMAQ solutions for N, P, Q in 3rd column
statistics_All_Model_MC_simulations(1:3, 3,:) = median(All_Model_MC_simulations(:,4 : 6 ,:));
% enter low (2.5) percentile of SMAQ solutions for N, P, Q in 4th column
statistics_All_Model_MC_simulations(1:3, 4,:) = prctile(All_Model_MC_simulations(:,4 : 6 ,:),2.5);
% enter high (97.5) percentile of SMAQ solutions for N, P, Q in 5th column
statistics_All_Model_MC_simulations(1:3, 5,:) = prctile(All_Model_MC_simulations(:,4 : 6 ,:),97.5);

%% save the simulation results
% write the parameters of the simulation into an excel sheet:

% noise, sweeps and runs of the simulation
Exp_parameters = array2table([n_sweeps, n_runs, noise], 'VariableNames',{'sweeps','runs','noiseSD'});

% ranges of quantal parameters modeled: 
Sim_parameter_ranges = [array2table(['N';'P';'Q'],'VariableNames',{'Q_var'}),...
    array2table([N_range(size(N_range,1)), N_range(size(N_range,2));...
    P_range(size(P_range,1)), P_range(size(P_range,2));...
    Q_range(size(Q_range,1)), Q_range(size(Q_range,2))],...
    'VariableNames',{'lower_bound','upper_bound'})];

if saving == 1 
    %save simulation results in .mat files
    save([savename ,'statistics_All_Model_MC_simulations.mat'],'statistics_All_Model_MC_simulations'); 
%   save([savename ,'All_Model_MC_simulations.mat'],'All_Model_MC_simulations');
    save([savename ,'histograms_database_N.mat'],'histograms_database_N');
    save([savename ,'histograms_database_P.mat'],'histograms_database_P');
    save([savename ,'histograms_database_Q.mat'],'histograms_database_Q');
    save([savename ,'histograms_database_models.mat'],'histograms_database_models');
    
    save([savename ,'QuantalParameterRanges.mat'],'QuParam_ranges');
    
    % save excel sheet
    writetable(Exp_parameters, [savename ,'MC_Simulation_parameters.xlsx'],'Sheet',1,'Range','A1:C2'); % noise, sweeps and runs of the simulation
    writetable(Sim_parameter_ranges, [savename ,'MC_Simulation_parameters.xlsx'],'Sheet',2,'Range','A1:C4'); % ranges of quantal parameters modeled 
end

end