%% SCRIPT FOR THE GENERAL EVALUATION OF A RECORDING FROM A PROBE (A.F. May 2024)
%  Inludes the evaluation of LFP, IMU and Single Units
%  Run within the Sorted_Units_AF Folder

%=== Load data
probe_number = 2;
load(['LFP_probe',num2str(probe_number),'.mat']);                   % LFP data, saved by loadTrodesLFP_AF_v1
load(['SU_kilosort4_outdir_probe',num2str(probe_number),'.mat']);   % SINGLE units, saved by loadSpikeData_AF_v1
load('IMU_data.mat');                                               % IMU data, saved by loadTrodesAnalog_AF_v1

%=== Parameters
Fs_LFP = 500;           % Sampling frequency for working on the LFP  
Fs_imu = NP_imu.Fs;     % Sampling frequency of the IMU data

%=== Transform the LFP into double and interpolate at 1kHz even sampling
t_LFP = linspace(red_out.t_ds(1),red_out.t_ds(end),(red_out.t_ds(end)-red_out.t_ds(1))*Fs_LFP)';
LFP = interp1(red_out.t_ds,double(red_out.lfp).*red_out.voltage_scaling,t_LFP);
N_channels = size(LFP,2);
LFP_mn = normalize(mean(LFP,2));

%=== Calculate and normalize the PSDs of the LFP
[psd,f_psd] = pwelch(LFP,[],[],[],Fs_LFP);  psd = psd./sum(psd,1);

%=== Interpolate LFP at the accelerometer time samples (useful later)
LFP_interp = interp1(t_LFP,LFP,NP_imu.t,'linear','extrap');  
[PS_LFP,freq_SG] = cwt(mean(LFP_interp,2),Fs_imu);

%=== Populate the s cell with the single units and calculate firing rate
NP_unit = table();
NP_unit = [NP_unit; out.good_units];
NP_unit.fr = cellfun(@(x) 1/mean(diff(x/1e6)),NP_unit.spikeTimes_usec);
NP_unit = NP_unit(NP_unit.fr<1,:);      % Exclude neurons with high-firing rate
n_cells = size(NP_unit,1);              % Number of units
s = cell(n_cells,1);                    % Single units' spikes
Rate = zeros(length(t_LFP),n_cells);    % Smoothed Firing Rate
for nc = 1:n_cells
    s{nc,1} = NP_unit.spikeTimes_usec{nc,1}/1e6;                    % Each row contains the time (in s) of the spikes from a single unit
    too_close = find(diff(s{nc,1})<0.001);                          % Clean up duplicated spikes
    if ~isempty(too_close),disp([num2str(numel(too_close)), ' cleaned spikes']);s{nc,1}(too_close+1) = [];end
    Rate(:,nc) = kernel_rate_AF_v1(s{nc,1},0.1,t_LFP);              % Define Rate
end

%=== Convert NP unit to structure
NP_unit = table2struct(NP_unit);
[~,correlation_sorted_ids] = sort([NP_unit.fr]);    % Get sorting based on firing frequency

%=== Calculate Population Spike Density
all_s = sort(vertcat(s{:,1}));
t_rate = [t_LFP(1):0.1:t_LFP(end)];
all_rate = kernel_rate_AF_v1(all_s,1,t_rate)/n_cells;

%% CROSS CORRELATIONS AND RANKING

% %=== Calculate cross correlations between cells
% cell_pairs = nchoosek(1:n_cells,2);
% num_pairs = size(cell_pairs,1);
% max_corr = zeros(num_pairs,1);
% parfor i=1:num_pairs
%     temp_corr = xcorr(Rate(:,cell_pairs(i,1)),Rate(:,cell_pairs(i,2)),round(Fs_LFP*.5));
%     [~,tmp_idx] = max(abs(temp_corr));
%     max_corr(i) = temp_corr(tmp_idx);
% end
% 
% %=== Find the best pair and then use it as seed for ranking correlations
% [~,top_cell_pair] = max(abs(max_corr));                                                                 % Get the pair with max correlation
% seed_cell = cell_pairs(top_cell_pair,1);                                                                % Get the first cell of the pair
% cell_id_for_corr_sorting = setdiff(cell_pairs(any(cell_pairs == seed_cell,2),:),seed_cell,'stable');    % Ids of all cells paired with the seed cell
% corr_vl_for_corr_sorting = max_corr(any(cell_pairs == seed_cell,2));                                    % Correlation values for the paired cells
% [~,tmp_idx] = sort(corr_vl_for_corr_sorting);                                                           % Now sort correlation values
% correlation_sorted_ids = [seed_cell; cell_id_for_corr_sorting(tmp_idx)];                                % Get the corresponding cell sequence
% 
%=== Check autocorrelation function
Cross_c = [];
max_psd = [];
for nc=1:n_cells
    [cross_c,ctrs_theta] = cross_correlogram_AF_v0(s{nc,1},s{nc,1},100,0.5);
    Cross_c = [Cross_c;cross_c'];
    [psd,f_psd] = pwelch(cross_c,[],[],[],1/0.5);  max_psd = [max_psd; max(psd)];
end
[~,tmp_idx] = sort(max_psd,'descend');
%tmp_idx = [1:n_cells]
figure('units','normalized','outerposition',[.4 .2 0.3 .5]);
tiledlayout(1,2,'TileSpacing','compact');
nexttile;   imagesc([ctrs_theta(1),ctrs_theta(end)],[],log(Cross_c(tmp_idx,:)));    colormap('viridis');
nexttile;   plot(ctrs_theta,smoothdata(nanmean(Cross_c,1),'movmean',5)); set(gca,'YScale', 'log');

%=== Try to run PCA
[coeff,score,latent] = pca(Rate);
theta = cart2pol(coeff(:,1),coeff(:,2));
[~,correlation_sorted_ids] = sort(theta);

%=== Find periods of high VS LOW activity
rate_th = otsuthresh(histcounts(all_rate));
[actv_vector,actv_num,actv_str,actv_stp] = RateSegm_AF_v0(all_rate',rate_th,100,0.1);

%=== Calculate for each cell median time spike to takeoff
t2act = zeros(n_cells,1);
for nc =1:n_cells
    tgt_seq = knnsearch(s{nc},t_rate(actv_str)');
    t2act(nc) = median(t_rate(actv_str)'-s{nc}(tgt_seq));
end
[~,correlation_sorted_ids] = sort(t2act);


%% Process IMU data to extract flights, look at average LFP

%=== Extract flight periods from accelerometer
a_abs_NP = vecnorm(NP_imu.acc,2,2);                         % Absolute acceleration
a_flt_NP = bandpass(a_abs_NP,[7 9],Fs_imu);                 % Filtered at the wing-beat frequency
[up,lo] = envelope(a_flt_NP,round(0.06*Fs_imu),'peak');     % Upper and lower envelopes 
env = normalize(up - lo,'range');                           % Amplitude of the envelope
env_th = otsuthresh(histcounts(env));                       % Threshold (based on Otsu method). Can be set at 0.35
wBeats = movsum(env>env_th,2*Fs_imu)>Fs_imu/5;              % Euristic criterion for flight detection

%=== Plot wing-beats and average LFP
figure('units','normalized','outerposition',[0 .2 1 .5]);
tiledlayout(2,1,'TileSpacing','compact');
ax(1) = nexttile;       area(NP_imu.t,wBeats*1,'FaceAlpha',0.3,'LineStyle','none');  hold on;
plot(NP_imu.t,normalize(a_flt_NP,'range',[-1 1]),'r');  plot(NP_imu.t,normalize(env,'range'),'k');  hold off;   ylabel('Accelerometer');        xticks([]);
ax(2) = nexttile;    area(NP_imu.t,wBeats*1,'FaceAlpha',0.3,'LineStyle','none');  hold on;
plot(t_LFP,normalize(mean(LFP,2),'range'),'b');                                                                 ylabel('Average LFP (norm)');   xticks([]);
linkaxes(ax,'x');   xlim('tight');  xlabel('Time (s)');

%=== Plot raster with sorted units, LFP and Spike Density
figure('units','normalized','outerposition',[0 .3 1 .6]);
tiledlayout(11,1,'TileSpacing','none');
cx(1) = nexttile(1,[4 1]);
for nc= 1:n_cells
    plot(s{correlation_sorted_ids(nc)}, nc*ones(size(s{correlation_sorted_ids(nc)})), 'k|','MarkerSize', round(n_cells*0.01));   hold on;          % Raster for each cell
end
area(NP_imu.t,wBeats*n_cells,0,'FaceColor',[0 0 1],'FaceAlpha',0.5,'LineStyle','none');    % Plot flights
ylim([0 n_cells]);  ylabel('Unit #');   xticks([]); set(gca,'TickLength',[0 0]);
cx(2) = nexttile(5,[1 1]);  plot(NP_imu.t,a_flt_NP);    ylabel('Accelerometer');    set(gca,'TickLength',[0 0]);
cx(3) = nexttile(6,[2 1]);  plot(t_LFP,LFP_mn); xticks([]); ylabel('LFP (norm)');    set(gca,'TickLength',[0 0]);
cx(4) = nexttile(8,[2 1]);
imagesc([NP_imu.t(1),NP_imu.t(end)],[freq_SG(1),freq_SG(end)],imgaussfilt(abs(PS_LFP),[2 10])); shading interp;  colormap(hot);
set(gca, 'YScale', 'log','YDir', 'normal','TickLength',[0 0]);   ylim([1 50]);  yticks([1 5 10 20 50]);    ylabel('Freq (Hz)');
cx(5) = nexttile(10,[2 1]);
plot(t_rate,all_rate,'k');  hold on; area(t_rate,actv_vector,'FaceAlpha',0.3,'LineStyle','none');   ylabel('Spike Density');   set(gca,'TickLength',[0 0]);
linkaxes(cx,'x');   xlim('tight');  xlabel('Time (s)');

%=== Calculate the cross-correlation between average LFP and Accelerometer
figure('units','normalized','outerposition',[.1 .2 .2 .3]);
[acc_r,acc_lags] = xcorr(mean(LFP_interp,2),a_abs_NP,round(3*Fs_imu));
plot(acc_lags/Fs_imu,acc_r);    hold on;    plot([0 0],ylim);   xlabel('Lag(s)');   ylabel('Cross-correlation');    xlim('tight');
[~,max_corr_loc] = max(acc_r);  title(['Delay: ',num2str(acc_lags(max_corr_loc)/Fs_LFP,1),'ms']);   

%% Inspect LFP on probe channels

%=== Warning sign if first two channels are not collinear
flag_diff = diff(red_out.channelPositions(:,2));
if flag_diff(1), warning('First two channels are not collinear!'); end

%=== Map channel locations into NP matrix (4 sites per row)
col_idx = ((red_out.channelPositions(:,1)/8)+5)/2;
row_idx = (red_out.channelPositions(:,2)-red_out.channelPositions(1,2))/20+1;

%=== Check difference between collinear channels
collinear_idx = find(~diff(red_out.channelPositions(:,2)));
med_abs_dev = zeros(size(collinear_idx));
for i=1:numel(collinear_idx)
    med_abs_dev(i) = mad(LFP(:,collinear_idx(i))-LFP(:,collinear_idx(i)+1));
end

%=== Initialize relevant matrices
NP_rms = NaN(row_idx(end),4);   NP_rms(sub2ind(size(NP_rms),row_idx,col_idx)) = rms(LFP,1);    
NP_avg = NaN(row_idx(end),4);   NP_avg(sub2ind(size(NP_avg),row_idx,col_idx)) = mean(LFP,1);
NP_rng = NaN(row_idx(end),4);   NP_rng(sub2ind(size(NP_rng),row_idx,col_idx)) = range(LFP,1);    
NP_cor = NaN(row_idx(end),4);   NP_cor(sub2ind(size(NP_cor),row_idx,col_idx)) = corr(a_abs_NP,LFP_interp);
NP_tta = NaN(row_idx(end),4);   NP_tta(sub2ind(size(NP_tta),row_idx,col_idx)) = sum(psd(f_psd>4 & f_psd<8,:));   
NP_rpl = NaN(row_idx(end),4);   NP_rpl(sub2ind(size(NP_rpl),row_idx,col_idx)) = sum(psd(f_psd>100 & f_psd<200,:));
NP_dif = NaN(row_idx(end),1);   NP_dif(unique(row_idx)) = med_abs_dev;

%=== Look at basic statistics across channels
figure('units','normalized','outerposition',[.7 .1 .3 .6]);
tiledlayout(1,7,'TileSpacing','compact');
nexttile;   imagesc(NP_avg,'AlphaData',~isnan(NP_avg));   set(gca,'YTickLabel',[]);   set(gca,'XTickLabel',[]);   title('LFP Avg.');   
nexttile;   imagesc(NP_rng,'AlphaData',~isnan(NP_rng));   set(gca,'YTickLabel',[]);   set(gca,'XTickLabel',[]);   title('LFP Range');
nexttile;   imagesc(NP_rms,'AlphaData',~isnan(NP_rms));   set(gca,'YTickLabel',[]);   set(gca,'XTickLabel',[]);   title('LFP RMS');
nexttile;   imagesc(NP_dif,'AlphaData',~isnan(NP_dif));   set(gca,'YTickLabel',[]);   set(gca,'XTickLabel',[]);   title('Coll. Diff'); 
nexttile;   imagesc(NP_cor,'AlphaData',~isnan(NP_cor));   set(gca,'YTickLabel',[]);   set(gca,'XTickLabel',[]);   title('Corr with Acc.');
nexttile;   imagesc(NP_tta,'AlphaData',~isnan(NP_tta));   set(gca,'YTickLabel',[]);   set(gca,'XTickLabel',[]);   title('Rel Theta Pwr');
nexttile;   imagesc(NP_rpl,'AlphaData',~isnan(NP_rpl));   set(gca,'YTickLabel',[]);   set(gca,'XTickLabel',[]);   title('Rel Ripple Pwr');
cb = colorbar;   set(cb, 'YTick', []); 

%=== Look at the LFP on all the channels (downsampled at 10Hz)
t_10Hz = [t_LFP(1):0.1:t_LFP(end)];
LFP_2look = interp1(t_LFP,LFP,t_10Hz,'linear','extrap');
scale_f = 0.5*median(range(LFP,1));
figure('units','normalized','outerposition',[0 0 .7 1]);
plot(t_10Hz,LFP_2look+flip([1:N_channels]*scale_f),'k');    hold on;
area(t_10Hz,logical(interp1(NP_imu.t,double(wBeats),t_10Hz))*N_channels*scale_f,'FaceAlpha',0.5,'LineStyle','none','FaceColor',[0 0 1]);
xlim('tight');  ylim('tight');  xlabel('Time (s)'); ylabel('LFP');
Y_labels = num2str(flip([row_idx,col_idx,red_out.channelMap'],1));                    
yticks([1:N_channels]*scale_f);   yticklabels(Y_labels);


