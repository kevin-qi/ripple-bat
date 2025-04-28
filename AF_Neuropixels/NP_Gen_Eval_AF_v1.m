function NP_Gen_Eval_AF_v1(probe_number)
%% SCRIPT FOR THE GENERAL EVALUATION OF A RECORDING FROM A PROBE (A.F. May 2024)
%  Inludes the evaluation of LFP, IMU and Single Units
%  Run within the Sorted_Units_AF Folder, or any folder that contains LFP_probe, SU_kilosort4_outdir_probe and IMU_data mat files

%=== Load data
load(['LFP_probe',num2str(probe_number),'.mat']);                   % LFP data, saved by loadTrodesLFP_AF_v1
load(['SU_kilosort4_outdir_probe',num2str(probe_number),'.mat']);   % SINGLE units, saved by loadSpikeData_AF_v1
load('IMU_data.mat');                                               % IMU data, saved by loadTrodesAnalog_AF_v1
if isfile([fileparts(fileparts(cd)),'\','Unique_Identifier.mat'])
    load([fileparts(fileparts(cd)),'\','Unique_Identifier.mat']);   % Load session identifier
elseif ~isempty(dir(fullfile(cd,'Extracted_Behavior_*')))
    behavioral_file = dir(fullfile(cd,'Extracted_Behavior_*'));     % Load session identifier (alternative way)
    load(behavioral_file.name);       
    unique_ID = options.unique_ID;  
else
    unique_ID = {'Dataset_x','zzzzz','yymmdd'};                     % Load session identifier (generic)
end

%=== Parameters
Fs_LFP = 500;               % Sampling frequency for working on the LFP  
Fs_imu = NP_imu.Fs;         % Sampling frequency of the IMU data
useMUA_for_spikeDensity =1; % If also using MUA for spike density
options.savedata = 1;       % If saving the results of analysis
options.savefigures = 1;    % Save Figures  
fig_count = 1;              % Save Figures

%=== Create analysis folder for storing the results
if options.savedata
    analysis_directory=fullfile(pwd,['Evaluation_Probe_',num2str(probe_number)]);
    if ~exist(analysis_directory,'dir')
        mkdir(analysis_directory);
    end
end

%% GENERAL ASSIGNMENTS AND PROCESSING

%=== Transform the LFP into double and interpolate at 1kHz or 500 Hz even sampling
t_LFP = linspace(red_out.t_ds(1),red_out.t_ds(end),(red_out.t_ds(end)-red_out.t_ds(1))*Fs_LFP)';
LFP = interp1(red_out.t_ds,double(red_out.lfp).*red_out.voltage_scaling,t_LFP);

%=== Get the number of channels and the mean LFP
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
%NP_unit = NP_unit(NP_unit.fr<1,:);      % Exclude neurons with high-firing rate
n_cells = size(NP_unit,1);              % Number of units
s = cell(n_cells,1);                    % Single units' spikes
Rate = zeros(length(t_LFP),n_cells);    % Smoothed Firing Rate
for nc = 1:n_cells
    s{nc,1} = NP_unit.spikeTimes_usec{nc,1}/1e6;                    % Each row contains the time (in s) of the spikes from a single unit
    too_close = find(diff(s{nc,1})<0.001);                          % Clean up duplicated spikes
    if ~isempty(too_close)
        %disp([num2str(numel(too_close)), ' cleaned spikes']);
        s{nc,1}(too_close+1) = [];
    end
    Rate(:,nc) = kernel_rate_AF_v1(s{nc,1},0.1,t_LFP);              % Define Rate
end

NP_unit = table2struct(NP_unit);                    % Convert NP unit to structure
[~,correlation_sorted_ids] = sort([NP_unit.fr]);    % Get sorting based on firing frequency

%=== Populate the mua cell with the single units and calculate firing rate
MUA_unit = table();
MUA_unit = [MUA_unit; out.mua_units];
MUA_unit.fr = cellfun(@(x) 1/mean(diff(x/1e6)),MUA_unit.spikeTimes_usec);
%NP_unit = NP_unit(NP_unit.fr<1,:);      % Exclude mua with high-firing rate
n_mua = size(MUA_unit,1);              % Number of mua units
s_mua = cell(n_mua,1);                    % Single units' spikes
Rate_mua = zeros(length(t_LFP),n_mua);    % Smoothed Firing Rate
for nc = 1:n_mua
    s_mua{nc,1} = MUA_unit.spikeTimes_usec{nc,1}/1e6;                    % Each row contains the time (in s) of the spikes from a single unit
    too_close = find(diff(s_mua{nc,1})<0.001);                          % Clean up duplicated spikes
    if ~isempty(too_close)
        %disp([num2str(numel(too_close)), ' cleaned spikes']);
        s_mua{nc,1}(too_close+1) = [];
    end
    Rate_mua(:,nc) = kernel_rate_AF_v1(s_mua{nc,1},0.1,t_LFP);              % Define Rate
end
MUA_unit = table2struct(MUA_unit);                    % Convert MUA unit to structure

%=== Calculate Population Spike Density
all_s = sort(vertcat(s{:,1}));
if useMUA_for_spikeDensity
all_s = sort([vertcat(s{:,1});vertcat(s_mua{:,1})]);
end
t_rate = [t_LFP(1):0.1:t_LFP(end)];
all_rate = kernel_rate_AF_v1(all_s,1,t_rate)/n_cells;

%=== Extract flight periods from accelerometer
a_abs_NP = vecnorm(NP_imu.acc,2,2);                         % Absolute acceleration
a_flt_NP = bandpass(a_abs_NP,[7 9],Fs_imu);                 % Filtered at the wing-beat frequency
[up,lo] = envelope(a_flt_NP,round(0.06*Fs_imu),'peak');     % Upper and lower envelopes 
env = normalize(up - lo,'range');                           % Amplitude of the envelope
env_th = otsuthresh(histcounts(env));                       % Threshold (based on Otsu method). Can be set at 0.35
wBeats = movsum(env>env_th,2*Fs_imu)>Fs_imu/5;              % Euristic criterion for flight detection

%=== Calculate the cross-correlation between average LFP and Accelerometer
[acc_r,acc_lags] = xcorr(mean(LFP_interp,2),a_abs_NP,round(3*Fs_imu));

%=== Warning sign if first two channels are not collinear
flag_diff = diff(red_out.channelPositions(:,2));
if flag_diff(1), warning('First two channels are not collinear!'); end

%=== Map channel locations into NP matrix (4 sites per row)
col_idx = ((red_out.channelPositions(:,1)/8)+5)/2;
row_idx = (red_out.channelPositions(:,2)-red_out.channelPositions(1,2))/20+1;

%=== Check difference between collinear channels
collinear_idx = find(~diff(red_out.channelPositions(:,2)));
if ~isempty(collinear_idx)
    med_abs_dev = zeros(size(collinear_idx));
    for i=1:numel(collinear_idx)
        med_abs_dev(i) = mad(LFP(:,collinear_idx(i))-LFP(:,collinear_idx(i)+1));
    end
else
    med_abs_dev = zeros(size(row_idx));
end

%=== Get average location of the single units and MUA and count them
unit_loc_x = red_out.channelPositions(knnsearch(red_out.channelPositions(:,1),cellfun(@(x) mean(x(:,1)),out.good_units.spikePos_um)),1);
unit_loc_y = red_out.channelPositions(knnsearch(red_out.channelPositions(:,2),cellfun(@(x) mean(x(:,2)),out.good_units.spikePos_um)),2);
mua_loc_x = red_out.channelPositions(knnsearch(red_out.channelPositions(:,1),cellfun(@(x) mean(x(:,1)),out.mua_units.spikePos_um)),1);
mua_loc_y = red_out.channelPositions(knnsearch(red_out.channelPositions(:,2),cellfun(@(x) mean(x(:,2)),out.mua_units.spikePos_um)),2);
unit_loc_row = (unit_loc_y-red_out.channelPositions(1,2))/20+1;
mua_loc_row = (mua_loc_y-red_out.channelPositions(1,2))/20+1;
[n_spikes,row_edges] = histcounts(unit_loc_row,(row_idx(1)-0.5):(row_idx(end)+0.5));
ch_spikes = n_spikes(ismember(row_edges(1:end-1)+0.5,unique(row_idx)'));
n_mua = histcounts(mua_loc_row,(row_idx(1)-0.5):(row_idx(end)+0.5));
ch_mua = n_mua(ismember(row_edges(1:end-1)+0.5,unique(row_idx)'));

%=== Create table for rapid assessment of the LFP config
Probe = table();
Probe.ch = red_out.channelMap';
Probe.xPos = red_out.channelPositions(:,1);
Probe.yPos = red_out.channelPositions(:,2);
Probe.col = ((Probe.xPos/8)+5)/2;
Probe.row = (Probe.yPos-red_out.yPos(1))/20+1;
Probe.depths = red_out.channelPositions(:,2);                   % This is a remnant of an old version of the function
Probe.theta = sum(psd(f_psd>4 & f_psd<8,:))'; 
Probe.ripple = sum(psd(f_psd>100 & f_psd<200,:))'; 
profile = varfun(@mean,Probe,'GroupingVariables','depths');
profile.spikes = ch_spikes';

%=== Initialize relevant matrices
NP_rms = NaN(row_idx(end),4);   NP_rms(sub2ind(size(NP_rms),row_idx,col_idx)) = rms(LFP,1);    
NP_avg = NaN(row_idx(end),4);   NP_avg(sub2ind(size(NP_avg),row_idx,col_idx)) = mean(LFP,1);
NP_rng = NaN(row_idx(end),4);   NP_rng(sub2ind(size(NP_rng),row_idx,col_idx)) = range(LFP,1);    
NP_cor = NaN(row_idx(end),4);   NP_cor(sub2ind(size(NP_cor),row_idx,col_idx)) = corr(a_abs_NP,LFP_interp);
NP_tta = NaN(row_idx(end),4);   NP_tta(sub2ind(size(NP_tta),row_idx,col_idx)) = Probe.theta;   
NP_rpl = NaN(row_idx(end),4);   NP_rpl(sub2ind(size(NP_rpl),row_idx,col_idx)) = Probe.ripple;
NP_dif = NaN(row_idx(end),1);   NP_dif(unique(row_idx)) = med_abs_dev;
NP_spk = NaN(row_idx(end),1);   NP_spk(unique(row_idx)) = ch_spikes;
NP_mua = NaN(row_idx(end),1);   NP_mua(unique(row_idx)) = ch_mua;

%=== Probe layout for visualization only
NP_matrix = zeros(480,4);
NP_matrix(sub2ind(size(NP_matrix),Probe.row,Probe.col)) = 1;
Probe_tip = NaN(15,4);  Probe_tip(1:5,4) = 0;   Probe_tip(6:10,3:4) = 0;    Probe_tip(11:15,2:4) = 0;   
NP_matrix = [Probe_tip;NP_matrix];

%% Process IMU data to extract flights, look at average LFP

%=== Visualize probe (with tip for clarity)
figure('units','normalized','outerposition',[.3 .3 .1 .6]);
tiledlayout(1,3,'TileSpacing','compact');
nexttile;   axis('off');    % Empty tile
nexttile;
imagesc(NP_matrix,'AlphaData',~isnan(NP_matrix));
set(gca,'YTickLabel',[],'XTickLabel',[],'YDir','normal'); axis('off');
nexttile;   axis('off');    % Empty tile
sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3},', Probe ',num2str(probe_number)],'Interpreter','none');
fig_count = saveFig(analysis_directory,cell2mat(unique_ID),fig_count,options.savefigures);

%=== Plot wing-beats and average LFP
figure('units','normalized','outerposition',[0 .2 1 .5]);
tiledlayout(2,1,'TileSpacing','compact');
ax(1) = nexttile;       area(NP_imu.t,wBeats*1,'FaceAlpha',0.3,'LineStyle','none');  hold on;
plot(NP_imu.t,normalize(a_flt_NP,'range',[-1 1]),'r');  plot(NP_imu.t,normalize(env,'range'),'k');  hold off;   ylabel('Accelerometer');    legend('','Signal','Envelope');        xticks([]);
ax(2) = nexttile;    area(NP_imu.t,wBeats*1,'FaceAlpha',0.3,'LineStyle','none');  hold on;
plot(t_LFP,normalize(mean(LFP,2),'range'),'b');                                                                 ylabel('Average LFP (norm)');   xticks([]);
linkaxes(ax,'x');   xlim('tight');  xlabel(['Time (s), LFP sampling: ',num2str(red_out.sampling_freq,'%.2f'),' Hz']);
sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3},', Probe ',num2str(probe_number)],'Interpreter','none');
fig_count = saveFig(analysis_directory,cell2mat(unique_ID),fig_count,options.savefigures);

%=== Plot raster with sorted units, LFP and Spike Density
figure('units','normalized','outerposition',[0 .3 1 .6]);
tiledlayout(11,1,'TileSpacing','none');
cx(1) = nexttile(1,[4 1]);
for nc= 1:n_cells
    plot(s{correlation_sorted_ids(nc)}, nc*ones(size(s{correlation_sorted_ids(nc)})), 'k|','MarkerSize', max(round(n_cells*0.01),1));   hold on;          % Raster for each cell
end
area(NP_imu.t,wBeats*n_cells,0,'FaceColor',[0 0 1],'FaceAlpha',0.5,'LineStyle','none');    % Plot flights
ylim([0 n_cells]);  ylabel('Unit #');   xticks([]); set(gca,'TickLength',[0 0]);
cx(2) = nexttile(5,[1 1]);  plot(NP_imu.t,a_flt_NP);    ylabel('Accelerometer');    set(gca,'TickLength',[0 0]);
cx(3) = nexttile(6,[2 1]);  plot(t_LFP,LFP_mn); xticks([]); ylabel('LFP (norm)');    set(gca,'TickLength',[0 0]);
cx(4) = nexttile(8,[2 1]);
imagesc([NP_imu.t(1),NP_imu.t(end)],[freq_SG(1),freq_SG(end)],imgaussfilt(abs(PS_LFP),[2 10])); shading interp;  colormap(hot);
set(gca, 'YScale', 'log','YDir', 'normal','TickLength',[0 0]);   ylim([1 50]);  yticks([1 5 10 20 50]);    ylabel('Freq (Hz)');
cx(5) = nexttile(10,[2 1]);
plot(t_rate,all_rate,'k');  ylabel('Spike Density');   set(gca,'TickLength',[0 0]);
linkaxes(cx,'x');   xlim('tight');  xlabel('Time (s)');
sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3},', Probe ',num2str(probe_number)],'Interpreter','none');
fig_count = saveFig(analysis_directory,cell2mat(unique_ID),fig_count,options.savefigures);

%=== Show the cross-correlation between average LFP and Accelerometer
figure('units','normalized','outerposition',[.1 .2 .2 .4]);
plot(acc_lags/Fs_imu,acc_r);    hold on;    plot([0 0],ylim);   xlabel('Lag(s)');   ylabel('Cross-correlation');    xlim('tight');
[~,max_corr_loc] = max(acc_r);  %title(['Delay: ',num2str(acc_lags(max_corr_loc)/Fs_LFP,1),'ms']); 
sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3},', Probe ',num2str(probe_number)],'Interpreter','none');
fig_count = saveFig(analysis_directory,cell2mat(unique_ID),fig_count,options.savefigures);

%% Inspect LFP on probe channels

%=== Look at basic statistics across channels
figure('units','normalized','outerposition',[.3 .1 .3 .6]);
tiledlayout(1,9,'TileSpacing','compact');
nexttile;   imagesc(NP_avg,'AlphaData',~isnan(NP_avg));   set(gca,'YTickLabel',[],'XTickLabel',[],'YDir','normal'); axis('off');   title('LFP Avg.');       
nexttile;   imagesc(NP_rng,'AlphaData',~isnan(NP_rng));   set(gca,'YTickLabel',[],'XTickLabel',[],'YDir','normal'); axis('off');   title('LFP Range');   
nexttile;   imagesc(NP_rms,'AlphaData',~isnan(NP_rms));   set(gca,'YTickLabel',[],'XTickLabel',[],'YDir','normal'); axis('off');   title('LFP RMS');     
nexttile;   imagesc(NP_dif,'AlphaData',~isnan(NP_dif));   set(gca,'YTickLabel',[],'XTickLabel',[],'YDir','normal'); axis('off');   title('Coll. Diff');  
nexttile;   imagesc(NP_cor,'AlphaData',~isnan(NP_cor));   set(gca,'YTickLabel',[],'XTickLabel',[],'YDir','normal'); axis('off');   title('Corr with Acc.');  
nexttile;   imagesc(NP_tta,'AlphaData',~isnan(NP_tta));   set(gca,'YTickLabel',[],'XTickLabel',[],'YDir','normal'); axis('off');   title('Rel Theta Pwr');   
nexttile;   imagesc(NP_rpl,'AlphaData',~isnan(NP_rpl));   set(gca,'YTickLabel',[],'XTickLabel',[],'YDir','normal'); axis('off');   title('Rel Ripple Pwr');  
nexttile;   imagesc(NP_spk,'AlphaData',~isnan(NP_spk));   set(gca,'YTickLabel',[],'XTickLabel',[],'YDir','normal'); axis('off');   title('Units'); 
nexttile;   imagesc(NP_mua,'AlphaData',~isnan(NP_mua));   set(gca,'YTickLabel',[],'XTickLabel',[],'YDir','normal'); axis('off');   title('MUA');            
cb = colorbar;   set(cb, 'YTick', []);  colormap('viridis'); 
sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3},', Probe ',num2str(probe_number)],'Interpreter','none');
fig_count = saveFig(analysis_directory,cell2mat(unique_ID),fig_count,options.savefigures);

%=== Look at basic statistics across channels
figure('units','normalized','outerposition',[.3 .1 .1 .5]);
scatter(normalize(profile.mean_theta,'range',[-1 1]), profile.depths,50,'r','square','filled');    hold on;
scatter(normalize(profile.mean_ripple,'range',[-1 1]),profile.depths,50,'b','^','filled');
scatter(normalize(profile.spikes,'range',[-1 1]),     profile.depths,50,'g','o','filled');
ylabel('Depth from tip (um)');    legend('Theta','Ripple','Spikes','Location', 'eastoutside');  xlim(1.5*[-1 1]); xticks([]);   ylim('tight');
sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3},', Probe ',num2str(probe_number)],'Interpreter','none');
fig_count = saveFig(analysis_directory,cell2mat(unique_ID),fig_count,options.savefigures);

%=== Look at the LFP on all the channels (downsampled at 10Hz)
t_10Hz = [t_LFP(1):0.1:t_LFP(end)];
LFP_2look = interp1(t_LFP,LFP,t_10Hz,'linear','extrap');
scale_f = 0.5*median(range(LFP,1));
figure('units','normalized','outerposition',[0 0 1 1]);
tiledlayout(10,1,'TileSpacing','compact');
ax(1) = nexttile(1,[9 1]);
plot(t_10Hz,LFP_2look+flip([1:N_channels]*scale_f),'k');    hold on;
area(t_10Hz,logical(interp1(NP_imu.t,double(wBeats),t_10Hz))*N_channels*scale_f,'FaceAlpha',0.1,'LineStyle','none','FaceColor',[0 0 1]);
ylim('tight');  ylabel('LFP');  xticks([]);
Y_labels = num2str(flip([row_idx,col_idx,red_out.channelMap'],1));                    
yticks([1:N_channels]*scale_f);   yticklabels(Y_labels);    set(gca, 'YDir','reverse');
ax(2) = nexttile;
plot(t_rate,all_rate,'k');
linkaxes(ax,'x');   xlim('tight');  xlabel('Time (s)'); box('off');   
sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3},', Probe ',num2str(probe_number)],'Interpreter','none');
fig_count = saveFig(analysis_directory,cell2mat(unique_ID),fig_count,options.savefigures);
