%% Load data

RPL_1 = load('Ripples_probe1.mat');
RPL_2 = load('Ripples_probe2.mat');
RPL_3 = load('Ripples_probe3.mat');

%% Plot cross correlations

min_corr = 0.3;
maxlag = 0.05;
binsize = 0.005;

figure('units','normalized','outerposition',[.3 .2 .3 .3]);
tiledlayout(1,3,'TileSpacing','compact');

[cross_corr12,bin_centers] = cross_correlogram_AF_v0(RPL_1.RPL.t(RPL_1.RPL.corr>min_corr),RPL_2.RPL.t(RPL_2.RPL.corr>min_corr),maxlag,binsize);
[cross_corr13,bin_centers] = cross_correlogram_AF_v0(RPL_1.RPL.t(RPL_1.RPL.corr>min_corr),RPL_3.RPL.t(RPL_3.RPL.corr>min_corr),maxlag,binsize);
[cross_corr23,bin_centers] = cross_correlogram_AF_v0(RPL_2.RPL.t(RPL_2.RPL.corr>min_corr),RPL_3.RPL.t(RPL_3.RPL.corr>min_corr),maxlag,binsize);

nexttile;   plot(bin_centers,cross_corr12,'k'); hold on;    plot([0 0],ylim,'r--'); title('Probes 1-2');    xlabel('Time lag (s)');  ylabel('Fraction');
nexttile;   plot(bin_centers,cross_corr13,'k'); hold on;    plot([0 0],ylim,'r--'); title('Probes 1-3');    xlabel('Time lag (s)');  ylabel('Fraction');
nexttile;   plot(bin_centers,cross_corr23,'k'); hold on;    plot([0 0],ylim,'r--'); title('Probes 2-3');    xlabel('Time lag (s)');  ylabel('Fraction');
