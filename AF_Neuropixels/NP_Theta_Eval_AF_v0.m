function NP_Theta_Eval_AF_v0(probe_number)
%% SCRIPT FOR THE GENERAL EVALUATION OF A RECORDING FROM A PROBE (A.F. May 2024)
%  Inludes the evaluation of LFP, IMU and Single Units
%  Run within the Sorted_Units_AF Folder, or any folder that contains LFP_probe, SU_kilosort4_outdir_probe and IMU_data mat files

%=== Load data
load(['LFP_probe',num2str(probe_number),'.mat']);                                   % LFP data, saved by loadTrodesLFP_AF_v1
load(['SU_kilosort4_outdir_probe',num2str(probe_number),'.mat']);                   % SINGLE units, saved by loadSpikeData_AF_v1
load('IMU_data.mat');                                                               % IMU data, saved by loadTrodesAnalog_AF_v1
bhvFile = dir('Extracted_Behavior_*');  load([bhvFile.folder '/' bhvFile.name]);    % Behavioral file
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

%=== Create analysis folder for storing the results
options.savedata = 0;                                   % If creating folder for saving the data or not
options.savefigures = 1;                                % Save Figures
fig_count = 1;                                          % Id of the first figure
if ~exist('folder_name') && options.savedata
    analysis_directory = [replace(pwd,'NP_datasets','NP_datasets_Analysis_Figures'),'\NP&BHv_Analysis_',datestr(now, 'yymmdd_HHMM')];
    if ~exist(analysis_directory,'dir');mkdir(analysis_directory);end
elseif exist('folder_name') && options.savedata
    analysis_directory = [replace(pwd,'NP_datasets',folder_name),'\NP&BHv_Analysis_',datestr(now, 'yymmdd_HHMM')];
    if ~exist(analysis_directory,'dir');mkdir(analysis_directory);end
elseif ~options.savedata
    analysis_directory = 'none';
    options.savefigures = 0;
end

%% GENERAL ASSIGNMENTS AND PROCESSING

%=== Transform the LFP into double and interpolate at even sampling
t_LFP = linspace(red_out.t_ds(1),red_out.t_ds(end),(red_out.t_ds(end)-red_out.t_ds(1))*Fs_LFP)';
LFP = interp1(red_out.t_ds,double(red_out.lfp).*red_out.voltage_scaling,t_LFP);

%=== Get the number of channels and the mean LFP (normalized)
N_channels = size(LFP,2);
LFP_mn = normalize(mean(LFP,2));

%=== Calculate and normalize the PSDs of the LFP
[psd,f_psd] = pwelch(LFP,[],[],[],Fs_LFP);  psd = psd./sum(psd,1);

%=== Interpolate LFP at the accelerometer time samples (useful later)
LFP_interp = interp1(t_LFP,LFP,NP_imu.t,'linear','extrap');  
[PS_LFP,freq_SG] = cwt(mean(LFP_interp,2),Fs_imu);

%=== Extract flight periods from accelerometer
a_abs_NP = vecnorm(NP_imu.acc,2,2);                         % Absolute acceleration
a_flt_NP = bandpass(a_abs_NP,[7 9],Fs_imu);                 % Filtered at the wing-beat frequency
[up,lo] = envelope(a_flt_NP,round(0.06*Fs_imu),'peak');     % Upper and lower envelopes 
env = normalize(up - lo,'range');                           % Amplitude of the envelope
env_th = otsuthresh(histcounts(env));                       % Threshold (based on Otsu method). Can be set at 0.35
wBeats = movsum(env>env_th,2*Fs_imu)>Fs_imu/5;              % Euristic criterion for flight detection
a_mdg = movmean(abs(a_abs_NP-1),Fs_imu*0.5);                % Mean deviation from g

%=== Transform start/stop samples of flight to the IMU sampling
f_smp = f_smp(t(f_smp(:,1))>0,:);
f_num = size(f_smp,1);
for i=1:f_num
   f_smp(i,1) = knnsearch_fast_AF_v0(NP_imu.t',t(f_smp(i,1)),0.05);
   f_smp(i,2) = knnsearch_fast_AF_v0(NP_imu.t',t(f_smp(i,2)),0.05);
end

%=== LFP and PSD during flight and during rest
LFP_r = LFP_interp(~wBeats,:);      [PSD_r,f_PSD_r] = pwelch(LFP_r,[],[],[],Fs_imu);  PSD_r = PSD_r./sum(PSD_r,1);  f_PSD_smpl_r = mean(diff(f_PSD_r));
LFP_f = LFP_interp( wBeats,:);      [PSD_f,f_PSD_f] = pwelch(LFP_f,[],[],[],Fs_imu);  PSD_f = PSD_f./sum(PSD_f,1);  f_PSD_smpl_f = mean(diff(f_PSD_f));
LFP_i = LFP_interp(a_mdg<0.04,:);   [PSD_i,f_PSD_i] = pwelch(LFP_i,[],[],[],Fs_imu);  PSD_i = PSD_i./sum(PSD_i,1);  f_PSD_smpl_i = mean(diff(f_PSD_i));

%=== PSD of the accelerometer signal during flight
[PSD_a,f_PSD_a] = pwelch(a_abs_NP(wBeats),[],[],[],Fs_imu);  PSD_a = PSD_a./sum(PSD_a,1);  f_PSD_smpl_a = mean(diff(f_PSD_a));

%=== Average, smooth and normalize PSDs
scaled_PSD_r = normalize(smoothdata(mean(PSD_r,2),'movmedian',1/f_PSD_smpl_r),'range');
scaled_PSD_f = normalize(smoothdata(mean(PSD_f,2),'movmedian',1/f_PSD_smpl_f),'range');
scaled_PSD_i = normalize(smoothdata(mean(PSD_i,2),'movmedian',1/f_PSD_smpl_i),'range');
scaled_PSD_a = normalize(smoothdata(mean(PSD_a,2),'movmedian',1/f_PSD_smpl_a),'range');

%=== Populate the s cell with the single units and calculate firing rate
NP_unit = table();
NP_unit = [NP_unit; out.good_units];
NP_unit.fr = cellfun(@(x) 1/mean(diff(x/1e6)),NP_unit.spikeTimes_usec);
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

%=== Calculate Population Spike Density
all_s = sort(vertcat(s{:,1}));
t_rate = [t_LFP(1):0.1:t_LFP(end)];
all_rate = kernel_rate_AF_v1(all_s,1,t_rate)/n_cells;

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
Probe.theta = sum(PSD_r(f_PSD_r>4 & f_PSD_r<11,:))'; 
Probe.theta_i = sum(PSD_i(f_PSD_i>4 & f_PSD_i<11,:))';
Probe.theta_f = sum(PSD_f(f_PSD_f>4 & f_PSD_f<11,:))'; 
Probe.ripple = sum(PSD_r(f_PSD_r>100 & f_PSD_r<200,:))'; 
profile = varfun(@mean,Probe,'GroupingVariables','depths');
profile.spikes = ch_spikes';

%=== Find channel with highest theta power during flight
[~,tmp_idx] = max(smoothdata(Probe.theta_f,'movmedian',4));
Probe.tmp = 0*Probe.theta_f;    Probe.tmp(tmp_idx) = 1;
opt_int = intersect(1:numel(Probe.tmp),tmp_idx+[-5:5]);

%=== Get power spectrum at optimal channel range
LFP_o = LFP_interp(wBeats,opt_int); [PSD_o,f_PSD_o] = pwelch(LFP_o,[],[],[],Fs_imu);  PSD_o = PSD_o./sum(PSD_o,1);  f_PSD_smpl_o = mean(diff(f_PSD_o));
scaled_PSD_o = normalize(smoothdata(mean(PSD_o,2),'movmedian',1/f_PSD_smpl_o),'range');

%=== Initialize relevant matrices
NP_rms = NaN(row_idx(end),4);   NP_rms(sub2ind(size(NP_rms),row_idx,col_idx)) = rms(LFP,1);    
NP_avg = NaN(row_idx(end),4);   NP_avg(sub2ind(size(NP_avg),row_idx,col_idx)) = mean(LFP,1);
NP_rng = NaN(row_idx(end),4);   NP_rng(sub2ind(size(NP_rng),row_idx,col_idx)) = range(LFP,1);    
NP_cor = NaN(row_idx(end),4);   NP_cor(sub2ind(size(NP_cor),row_idx,col_idx)) = corr(a_abs_NP,LFP_interp);
NP_tta = NaN(row_idx(end),4);   NP_tta(sub2ind(size(NP_tta),row_idx,col_idx)) = Probe.theta;  
NP_ttf = NaN(row_idx(end),4);   NP_ttf(sub2ind(size(NP_ttf),row_idx,col_idx)) = Probe.theta_f;   
NP_tti = NaN(row_idx(end),4);   NP_tti(sub2ind(size(NP_tti),row_idx,col_idx)) = Probe.theta_i;   
NP_tmp = NaN(row_idx(end),4);   NP_tmp(sub2ind(size(NP_tmp),row_idx,col_idx)) = Probe.tmp;
NP_rpl = NaN(row_idx(end),4);   NP_rpl(sub2ind(size(NP_rpl),row_idx,col_idx)) = Probe.ripple;
NP_dif = NaN(row_idx(end),1);   NP_dif(unique(row_idx)) = med_abs_dev;
NP_spk = NaN(row_idx(end),1);   NP_spk(unique(row_idx)) = ch_spikes;
NP_mua = NaN(row_idx(end),1);   NP_mua(unique(row_idx)) = ch_mua;

%=== Probe layout for visualization only
NP_matrix = zeros(480,4);
NP_matrix(sub2ind(size(NP_matrix),Probe.row,Probe.col)) = 1;
Probe_tip = NaN(15,4);  Probe_tip(1:5,4) = 0;   Probe_tip(6:10,3:4) = 0;    Probe_tip(11:15,2:4) = 0;   
NP_matrix = [Probe_tip;NP_matrix];

%% Inspect LFP on probe channels

%=== Look at basic statistics across channels
figure('units','normalized','outerposition',[.3 .1 .4 .6]);
tiledlayout(1,11,'TileSpacing','compact');
nexttile;   imagesc(NP_avg,'AlphaData',~isnan(NP_avg));   set(gca,'YTickLabel',[],'XTickLabel',[],'YDir','normal'); axis('off');   title('LFP Avg.');       
nexttile;   imagesc(NP_rng,'AlphaData',~isnan(NP_rng));   set(gca,'YTickLabel',[],'XTickLabel',[],'YDir','normal'); axis('off');   title('LFP Range');   
nexttile;   imagesc(NP_rms,'AlphaData',~isnan(NP_rms));   set(gca,'YTickLabel',[],'XTickLabel',[],'YDir','normal'); axis('off');   title('LFP RMS');     
nexttile;   imagesc(NP_dif,'AlphaData',~isnan(NP_dif));   set(gca,'YTickLabel',[],'XTickLabel',[],'YDir','normal'); axis('off');   title('Coll. Diff');  
nexttile;   imagesc(NP_cor,'AlphaData',~isnan(NP_cor));   set(gca,'YTickLabel',[],'XTickLabel',[],'YDir','normal'); axis('off');   title('Corr with Acc.');  
nexttile;   imagesc(NP_tta,'AlphaData',~isnan(NP_tta));   set(gca,'YTickLabel',[],'XTickLabel',[],'YDir','normal'); axis('off');   title('Rel Theta (R) Pwr');   
nexttile;   imagesc(NP_ttf,'AlphaData',~isnan(NP_ttf));   set(gca,'YTickLabel',[],'XTickLabel',[],'YDir','normal'); axis('off');   title('Rel Theta (F) Pwr'); 
nexttile;   imagesc(NP_tti,'AlphaData',~isnan(NP_tti));   set(gca,'YTickLabel',[],'XTickLabel',[],'YDir','normal'); axis('off');   title('Rel Theta (I) Pwr');   
nexttile;   imagesc(NP_tmp,'AlphaData',~isnan(NP_tmp));   set(gca,'YTickLabel',[],'XTickLabel',[],'YDir','normal'); axis('off');   title('Max Ch');
nexttile;   imagesc(NP_rpl,'AlphaData',~isnan(NP_rpl));   set(gca,'YTickLabel',[],'XTickLabel',[],'YDir','normal'); axis('off');   title('Rel Ripple Pwr');  
nexttile;   imagesc(NP_spk,'AlphaData',~isnan(NP_spk));   set(gca,'YTickLabel',[],'XTickLabel',[],'YDir','normal'); axis('off');   title('Units'); 
cb = colorbar;   set(cb, 'YTick', []);  colormap('viridis'); 
sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3},', Probe ',num2str(probe_number)],'Interpreter','none');
fig_count = saveFig(analysis_directory,[cell2mat(unique_ID)],fig_count,options.savefigures);

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


%% LOOK AT PSD AND LFP

%=== Plot the PSDs
figure('units','normalized','outerposition',[0.1 0.3 0.25 0.3]);
tiledlayout(1,2,'TileSpacing','tight');
nexttile;   plot(f_PSD_r,scaled_PSD_r,'LineWidth',2);    hold on;     plot(f_PSD_o,scaled_PSD_o,'LineWidth',2); plot(f_PSD_i,scaled_PSD_i,'LineWidth',2);  plot(f_PSD_a,scaled_PSD_a,'LineWidth',2);    
rectangle('Position',[4 0 7 1],'FaceColor',[0 0 0 0.2],'EdgeColor','none');
xlim([0 20]);   xlabel('Frequency (Hz)');   ylabel('PSD (Norm.)');  legend('LFP (Rest)','LFP (Flight)','LFP (Immobile)','Accelerometer (Flight)');
nexttile;   plot(f_PSD_r,scaled_PSD_r,'LineWidth',2);    hold on;     plot(f_PSD_o,scaled_PSD_o,'LineWidth',2); plot(f_PSD_i,scaled_PSD_i,'LineWidth',2);  plot(f_PSD_a,scaled_PSD_a,'LineWidth',2);    
rectangle('Position',[4 0 7 1],'FaceColor',[0 0 0 0.2],'EdgeColor','none');
[~,tmp_idx] = max(scaled_PSD_a);    plot(f_PSD_a(tmp_idx)*[1 1],ylim,'k--');
xlim([3 12]);   xlabel('Frequency (Hz)');   ylabel('PSD (Norm.)');  set(gca, 'YScale', 'log');
%ylim([0 max(scaled_PSD_f(f_PSD_f>4 & f_PSD_f<11))]);   
sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3},', Probe ',num2str(probe_number)],'Interpreter','none');  
fig_count = saveFig(analysis_directory,[cell2mat(unique_ID)],fig_count,options.savefigures);

%=== Accumulate the LFP for each flight, aligned to takeoff and landing
LFP_tht = bandpass(mean(LFP_interp(:,opt_int),2),[4 11],Fs_imu);
LFP_ave = mean(LFP_interp(:,opt_int),2);
[~,tmp_idx] = sort(diff(f_smp,1,2));
int_t = [-3 6]; smht = [0.5 1]; sat_c = [2 98];
interval = [round(int_t(1)*Fs_imu) : round(int_t(2)*Fs_imu)];
LFP_htmp = zeros(f_num,numel(interval),2);
LTH_htmp = zeros(f_num,numel(interval),2);
ACC_htmp = zeros(f_num,numel(interval),2);
for i=1:f_num
    LFP_htmp(i,:,1) = LFP_ave(f_smp(tmp_idx(i),1)+interval)';
    LTH_htmp(i,:,1) = LFP_tht(f_smp(tmp_idx(i),1)+interval)';
    ACC_htmp(i,:,1) = a_abs_NP(f_smp(tmp_idx(i),1)+interval)';
    LFP_htmp(i,:,2) = LFP_ave(f_smp(tmp_idx(i),2)+interval)';
    LTH_htmp(i,:,2) = LFP_tht(f_smp(tmp_idx(i),2)+interval)';
    ACC_htmp(i,:,2) = a_abs_NP(f_smp(tmp_idx(i),2)+interval)';
end

%=== Plot the LFP for each flight, aligned to takeoff and landing
figure('units','normalized','outerposition',[0.4 0.3 0.35 0.5]);
tiledlayout(2,3,'TileSpacing','tight');
ax(1) = nexttile;   imagesc(int_t,[],imgaussfilt(ACC_htmp(:,:,1),smht),prctile(ACC_htmp,sat_c,'all')');  colormap(ax(1),gray);  hold on;      
plot([0 0],ylim,'k--');   xlabel('Time to takeoff (s)');  ylabel('Sorted Flight #'); title('Accelerometer');
ax(2) = nexttile;   imagesc(int_t,[],imgaussfilt(LFP_htmp(:,:,1),smht),prctile(LFP_htmp,sat_c,'all')');  colormap(ax(2),redblue); hold on;    
plot([0 0],ylim,'k--');   xlabel('Time to takeoff (s)');  ylabel('Sorted Flight #'); title('LFP (All Freq.)');
ax(3) = nexttile;   imagesc(int_t,[],imgaussfilt(LTH_htmp(:,:,1),smht),prctile(LTH_htmp,sat_c,'all')');  colormap(ax(3),redblue); hold on;      colorbar;  
plot([0 0],ylim,'k--');   xlabel('Time to takeoff (s)');  ylabel('Sorted Flight #'); title('LFP (4-11 Hz)');
ax(4) = nexttile;   imagesc(int_t,[],imgaussfilt(ACC_htmp(:,:,2),smht),prctile(ACC_htmp,sat_c,'all')');  colormap(ax(4),gray);  hold on;         
plot([0 0],ylim,'k--');   xlabel('Time to landing (s)');  ylabel('Sorted Flight #'); 
ax(5) = nexttile;   imagesc(int_t,[],imgaussfilt(LFP_htmp(:,:,2),smht),prctile(LFP_htmp,sat_c,'all')');  colormap(ax(5),redblue); hold on;    
plot([0 0],ylim,'k--');   xlabel('Time to landing (s)');  ylabel('Sorted Flight #'); 
ax(6) = nexttile;   imagesc(int_t,[],imgaussfilt(LTH_htmp(:,:,2),smht),prctile(LTH_htmp,sat_c,'all')');  colormap(ax(6),redblue); hold on;    
plot([0 0],ylim,'k--');   xlabel('Time to landing (s)');  ylabel('Sorted Flight #'); 
sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3},', Probe ',num2str(probe_number)],'Interpreter','none');  colorbar;
fig_count = saveFig(analysis_directory,[cell2mat(unique_ID)],fig_count,options.savefigures);

%=== Plot the accelerometer and the theta 
figure('units','normalized','outerposition',[0 0.3 1 0.4]);
tiledlayout(3,1,'TileSpacing','tight');
bx(1) = nexttile;       plot(NP_imu.t,a_flt_NP);    
bx(2) = nexttile;       plot(NP_imu.t,mean(LFP_interp(:,opt_int),2));  
bx(3) = nexttile;       plot(NP_imu.t,LFP_tht);
linkaxes(bx,'x');


end