function[MC_results] = MonteCarlo_BinomialModel_SMAQ(N, P, Q, noise, Gmean, n_sweeps, n_runs)
%% Monte-Carlo Simulation of binomial model and SMAQ
% this function simulates n_runs histograms of n_sweeps each, using N,P,Q, noise, then SMAQ is used on each of the histograms

        % model generates n_runs histograms (dim1) of n_sweeps (dim2)
        % adding Gaussian noise to each observation
        model = Q * binornd(N,P,[n_runs,n_sweeps]) + normrnd(Gmean,noise,[n_runs,n_sweeps]);

        %% get statistics of the generated distribution
        % the stats of each run of the simulation are computed
        mu  = mean(model,2);
        sigma = std(model,0,2);     
        gamma = skewness(model,0,2);      % BIAS-CORRECTION ACTIVATED, IMPORTANT!     skewness over 2nd dim -> sweeps of histograms      

        %% use SMAQ to find N,P,Q for that distribution
        % on each run we perform SMAQ to recover N,P,Q 
        % .-wise multiplication and subtraction
        Nsim = - ((mu.^2) ./ ((gamma .* mu - sigma).*sigma));
        Psim = (gamma.*mu - sigma) ./ (gamma.*mu - 2.*sigma);
        Qsim = ( (-gamma .* mu .* sigma + (2 .* sigma.^2)) ./ mu);
            
        %% saving all results in a matrix to be exported back to the script:   
        MC_results(:,1:6) = [mu, sigma, gamma, Nsim, Psim, Qsim];
end