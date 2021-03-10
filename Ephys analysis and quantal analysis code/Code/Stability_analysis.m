function [stable_EPSPs,EPSP_STATS] = Stability_analysis(loadpathEPSPdata, savepath, filename, blocksize, mean_stable, SD_stable, min_initial_stable, filled, save)
%% This protocol looks 1st for 100 stable sweeps RELATIVE to the 1st bin, if this is not the case, then (2nd) the longest consecutive epoch is selected.
% developed for the SSN meeting 2018, gave the best results in terms of histograms and stable sweeps per histogram

%% load data
data = xlsread(loadpathEPSPdata);
EPSP = data(:,2);
noise = data(:,3);
sweeps = (1 : length(EPSP))';

%% calculate EPSP_STATS and SD for blocks of data
n_blocks = floor(length(EPSP)/blocksize);   % blocksize determines how many blocks you have, floor rounds down

%EPSP array is split into a matrix with column length of the block length 
Y = reshape(EPSP(1:(n_blocks*blocksize)),[blocksize,n_blocks]);

% statistics on the blocks
M = mean(Y,1);      % compute the mean column-wise
SD = std(Y,1);      % compute the SD column-wise
SEM =  std(Y,1)/sqrt(blocksize);
X = (blocksize/2 : blocksize : (n_blocks*blocksize) - blocksize/2);     % center of each bin/block

% write mean, SD, and CV^-2 into new matrix
EPSP_STATS = [];
EPSP_STATS(:,1) = X'; % block center
EPSP_STATS(:,2) = M'; % EPSP_STATS
EPSP_STATS(:,3) = SD'; % SD
EPSP_STATS(:,4) = SEM'; % SEM


%% Test for stability relative to a reference block, which moves across the matrix (every block is reference block once)
stability_matrix = zeros(n_blocks,n_blocks);

for reference_block = 1 : n_blocks      % every block becomes ref block once
    for n = 1 : n_blocks                % comparing each block to ref block
        if abs(EPSP_STATS(n,2) - EPSP_STATS(reference_block,2)) <= mean_stable * EPSP_STATS(reference_block,4)...
                && abs(EPSP_STATS(n,3) - EPSP_STATS(reference_block,3)) <= SD_stable * EPSP_STATS(reference_block,3)
                        
                stability_matrix(n,reference_block) = 1;    % if the block is stable relative to the ref block, you enter a 1
        else
        end
    end
end

%% Mark all maximum stable blocks next to the reference block
% finds the blocks neighboring the reference block that are stable, only consecutive sequences of stable sweeps should be used in the further analysis
consec_stability_matrix = zeros(n_blocks,n_blocks);

for reference_block = 1 : n_blocks      % going through all blocks
    
    % starting from the reference block (which is stable with regard to itself), move forward and backward, and enter 1s in the new
    % matrix when the blocks were stable, this is stopped, as soon as you encounter the 1st unstable block
    % -> gives us the CONSECUTIVE stable blocks around the reference block
    for n = reference_block : 1 : n_blocks  % going forward from reference block to end
        if stability_matrix(n,reference_block) == 1
            consec_stability_matrix(n,reference_block) = 1;
        else
            break
        end

    end
    
    for n = reference_block : -1 : 1        % going backward from reference block to beginning (1)
        if stability_matrix(n,reference_block) == 1
            consec_stability_matrix(n,reference_block) = 1;
        else
            break
        end

     end
    
end

S = sum(consec_stability_matrix,1); % sum up the stable blocks per column
F = find(S==max(S));                % locate all the blocks with the max stable sweeps


%% ------------------------------------------------------------------------
%% SELECTION OF STABLE EPOCH 
stable_EPSPs = zeros(length(EPSP),1);   % list of sweeps, with 1 indicating they are stable, this matrix is used by the other scripts later on

%% [1] if the 1st bin has enough consecutive sweeps (usually 100), it is selected as the stable epoch
if S(1) >= min_initial_stable
    EPSP_STATS(:,5) = consec_stability_matrix(:,1);  %F(1) is 1st entry in F, that gives the 1st max stable epoch
        % the 5th column of EPSP_STATS becomes the indicator which blocks are selected as 'stable'
    
    for f = 1 : n_blocks    % run through all blocks and expand the entry (0 or 1) that it has to the stable_EPSP array
        stable_EPSPs((f-1) * blocksize + 1 : f * blocksize, 1) = EPSP_STATS(f,5);
            % the 1s and 0s get expanded into the stable_EPSP array, so for each sweep you get an entry of 0 or 1...
    end
    
%% [2] If not, just take the longest consecutive stable epoch  
else 
    EPSP_STATS(:,5) = consec_stability_matrix(:,F(1));  %F(1) is 1st entry in F; F gives all columns with the max stable epochs (if there are more than one)
    for f = 1 : n_blocks    % run through all blocks and expand the entry (0 or 1) that it has to the stable_EPSP array
        stable_EPSPs((f-1) * blocksize + 1 : f * blocksize, 1) = EPSP_STATS(f,5);
    end
    
end

%% ------------------------------------------------------------------------
% Plotting

%% preparation for plotting
ind_stable = stable_EPSPs == 1;     % converts array to ind-array
n_sweeps = (1:length(EPSP))';       % becomes X in plotting

%% PLOTTING
figure('Name',[filename '_SSN_Stability'],'color',[1.00, 1.00, 1.00]);
set(0,'DefaultAxesFontSize',16);
set(gca, 'LineWidth', 3,'box','off','TickDir','out','FontWeight', 'bold');
hold on;
xlabel('count');
ylabel('amplitude [mV]');

%% Plotting UNSTABLE noise and EPSPs
% unstable noise
i = scatter(n_sweeps(~ind_stable,1), noise(~ind_stable,1),100);
i.MarkerEdgeColor = [0.8 0.8 0.8];
if filled == 1
    i.MarkerFaceColor = [0.95 0.95 0.95];
    else
end
i.LineWidth = 2;

% unstable EPSPs
i = scatter(n_sweeps(~ind_stable,1), EPSP(~ind_stable,1),100); 
i.MarkerEdgeColor = [0.65 0.65 0.65];
if filled == 1
    i.MarkerFaceColor = [0.8 0.8 0.8];
    else
end
i.LineWidth = 2;

%% Plotting STABLE noise and EPSPs
% stable noise
i = scatter(n_sweeps(ind_stable,1), noise(ind_stable,1),100); 
i.MarkerEdgeColor = [0.5 0.5 0.5];
if filled == 1
        i.MarkerFaceColor = [0.85 0.85 0.85];
    else
end
i.LineWidth = 2;

% stable EPSPs
i = scatter(n_sweeps(ind_stable,1), EPSP(ind_stable,1),100); 
i.MarkerEdgeColor = [0.1 0.1 0.1];
if filled == 1
    i.MarkerFaceColor = [0.4 0.4 0.4];
    else
end
i.LineWidth = 2;

%% ERRORBARS for EPSPs give the STANDARD DEVIATION
% 68.27% of values lie within +1 AND -1 SD from mean (so within a range of 2 SD)
% -> errorbar plots one SD in each direction (up and down)

% stable EPSPs
ind_stable_blocks = EPSP_STATS(:,5) == 1;
e = errorbar(EPSP_STATS(ind_stable_blocks,1), EPSP_STATS(ind_stable_blocks,2), EPSP_STATS(ind_stable_blocks,3),'o');
    %        position of bin center           EPSP mean                        standard deviation
e.LineWidth = 4; 
e.CapSize = 20;
e.MarkerSize = 12;
e.Color = [0.2 0.7 0.9];
if filled == 1
    e.MarkerFaceColor = [0.4 0.95 1];
    else
end

e = errorbar(EPSP_STATS(~ind_stable_blocks,1), EPSP_STATS(~ind_stable_blocks,2), EPSP_STATS(~ind_stable_blocks,3),'o');
e.LineWidth = 4; 
e.CapSize = 20;
e.MarkerSize = 12;
e.Color = [0.6 0.8 0.9];
if filled == 1
    e.MarkerFaceColor = [0.6 0.9 1];
    else
end

% errorbar for reference block that gave the most stable epochs
% if the 1st 100 sweeps were selected, the 1st block was used as reference block for the histogram and the reference (plotted in RED)
% block that gave the max number of stable sweeps was likely a different one (the GREEN errorbars will not be the 1st block, but somewhere else...)

if S(1) >= min_initial_stable   % if the 1st bin is ref block, plot it in different color, but also plot the ref block that gave the most stable sweeps
    % bin that gave the most stable sweeps in green as usual
    e = errorbar(EPSP_STATS(F(1),1), EPSP_STATS(F(1),2), EPSP_STATS(F(1),3),'o');
    e.LineWidth = 4; 
    e.CapSize = 20;
    e.MarkerSize = 12;
    e.Color = [0.1 0.8 0.1];
    
    % 1st bin in red (gives rise to 1st 100 stable sweeps), plotted on top
    f = errorbar(EPSP_STATS(1,1), EPSP_STATS(1,2), EPSP_STATS(1,3),'o');
    f.LineWidth = 4; 
    f.CapSize = 20;
    f.MarkerSize = 12;
    f.Color = [0.9 0.2 0];    
    
    if filled == 1
        e.MarkerFaceColor = [0.5 1 0.5];
        f.MarkerFaceColor = [1 0.5 0.5];
        else
    end    
    
else  % if the 1st bin is not ref block, only plot the ref bin in green
    e = errorbar(EPSP_STATS(F(1),1), EPSP_STATS(F(1),2), EPSP_STATS(F(1),3),'o');
    e.LineWidth = 4; 
    e.CapSize = 20;
    e.MarkerSize = 12;
    e.Color = [0.1 0.8 0.1];
    if filled == 1
        e.MarkerFaceColor = [0.5 1 0.5];
        else
    end
end

%% savefigure 
if save == 1
    mkdir([savepath 'Stability_Analysis']);
    print([savepath 'Stability_Analysis/' filename '_Stability_Graph'] , '-painters','-depsc');
end

end

