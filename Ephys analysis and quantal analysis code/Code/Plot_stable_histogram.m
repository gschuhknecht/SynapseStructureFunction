function [] = Plot_stable_histogram(loadpathEPSPdata, savepath, filename, binwidth,movmean_winsize, stable_epoch, save)
%% plot histogram for the stable data
% this script plots the selected stable sweeps and noise as histograms with a range of bin sizes, which need to be specified
% also a moving average is plotted, the points of which need to be specified too
% IMPORTANT for moving average: the matlab function movmean automatically truncated the window size at the endpoints when there are not enough elements
% to fill the window. When the window is truncated, the average is taken over only the elements that fill the window.
% -> at the end points, there is truncation, and the mean does not fall down to zero

%% setting parameters
data = xlsread(loadpathEPSPdata);
EPSP = data(:,2);
noise = data(:,3);

stable_EPSPs = EPSP(stable_epoch == 1);
stable_noise = noise(stable_epoch == 1);

figure('Name',[filename],'color',[1.00, 1.00, 1.00],'Position',[0 0 1000 1000]);
counter = 1;

    for a = binwidth % running through the different binwidths
        
        %% calculate moving average
        [H,edges] = histcounts(stable_EPSPs,'Binwidth',a);
        S = movmean(H,movmean_winsize(counter));
        centers = edges(1:end-1); % S is plotted using stairs, delete the last edge of the histogram to get synchronous plotting
        
        % important: duplicate the last value of the moving mean, so you get a histogram that is the same length as the EPSP histogram
        % (for cosmetic reasons, stairs is plotted differently to histogram)
        S  = [S S(end)];        
        centers  = [centers centers(end) + a];

        %% Plotting
        subplot(2,ceil(length(binwidth)/2),counter)
        set(0,'DefaultAxesFontSize',16);
        set(gca, 'LineWidth', 3,'box','off','TickDir','out','FontWeight', 'bold');
        hold on;
        xlabel('amplitude [mV]');
        ylabel('count');
        title({['bin size ' num2str(a) 'mV'] ; [num2str(movmean_winsize(counter)) '. mov average']});
        
        % noise histogram
        n = histogram(stable_noise,'Binwidth',a,'Linewidth',2);
        n.FaceColor = [0.85 0.85 0.85];
        n.EdgeColor = [0.5 0.5 0.5];

        % EPSP histogram
        e = histogram(stable_EPSPs,'Binwidth',a,'Linewidth',2);
        e.FaceColor = [0.2 0.2 0.2];
        e.EdgeColor = [0.1 0.1 0.1];

        b = stairs(centers, S,'LineWidth',3);
        b.Color = [1 0.5 0.1];

        counter = counter + 1;
        
    end 
        if save == 1
            mkdir([savepath 'Stable_Histograms']);
            print([savepath 'Stable_Histograms/' filename '_Stable_Histograms'] , '-painters','-depsc');
        end
end

