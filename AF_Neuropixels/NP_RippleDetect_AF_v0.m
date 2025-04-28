function NP_RippleDetect_AF_v0(probe_number)
%% Function for detecting ripples from NP recordings. First draft by AF on May 2024

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
Fs = 1e3;                                       % Sampling frequency for working on the LFP  
Theta_band = [4 8];                             % Theta  frequency band
Ripple_band = [100 200];                        % Ripples frequency band
Ripple_th = 3;                                  % Threshold for zscored ripple power
int50ms = [round(Fs*-0.05):round(Fs*0.05)];     % Samples corresponding to a 50 ms interval
int100ms = [round(Fs*-0.1):round(Fs*0.1)];      % Samples corresponding to a 100 ms interval
options.savedata = 1;                           % If saving the results of analysis
options.savefigures = 1;                        % Save Figures  
fig_count = 1;                                  % Save Figures

%=== Create analysis folder for storing the results
if options.savedata
    analysis_directory=fullfile(pwd,['Evaluation_Probe_',num2str(probe_number)]);
    if ~exist(analysis_directory,'dir')
        mkdir(analysis_directory);
    end
end

%=== Transform the LFP into double and interpolate at 1kHz even sampling
t = linspace(red_out.t_ds(1),red_out.t_ds(end),(red_out.t_ds(end)-red_out.t_ds(1))*Fs);
LFP = interp1(red_out.t_ds,double(red_out.lfp).*red_out.voltage_scaling,t);
N_channels = size(LFP,2);
scale_f = 0.5*median(range(LFP,1));

%=== Calculate and normalize the PSDs of the LFP
[psd,f_psd] = pwelch(LFP,[],[],[],Fs);  psd = psd./sum(psd,1);

%=== Extract flight periods from accelerometer
a_abs_NP = vecnorm(NP_imu.acc,2,2);                         % Absolute acceleration
a_flt_NP = bandpass(a_abs_NP,[7 9],NP_imu.Fs);              % Filtered at the wing-beat frequency
[up,lo] = envelope(a_flt_NP,round(0.06*NP_imu.Fs),'peak');  % Upper and lower envelopes 
env = normalize(up - lo,'range');                           % Amplitude of the envelope
env_th = otsuthresh(histcounts(env));                       % Threshold (based on Otsu method). Can be set at 0.35
wBeats = movsum(env>env_th,2*NP_imu.Fs)>NP_imu.Fs/5;        % Euristic criterion for flight detection

%=== Create table for rapid assessment of the LFP config
Probe = table();
Probe.idx = [1:size(LFP,2)]';                                               % Row in the LFP matrix
Probe.ch = red_out.channelMap';                                             % Probe Channel 
Probe.xPos = red_out.channelPositions(:,1);                                 % Horizontal Position
Probe.yPos = red_out.channelPositions(:,2);                                 % Vertical Position (from the tip!)
Probe.col = ((Probe.xPos/8)+5)/2;                                           % Probe Column (1 to 4)
Probe.row = (Probe.yPos-red_out.yPos(1))/20+1;                              % Probe Row (from the tip!!) 
Probe.rms = rms(LFP,1)';                                                    % RMS of the LFP
Probe.RPL_power = sum(psd(f_psd>Ripple_band(1) & f_psd<Ripple_band(2),:))'; % Relative Ripple power
Probe.Tht_power = sum(psd(f_psd>Theta_band(1) & f_psd<Theta_band(2),:))';   % Relative Theta power

%=== Probe layout for visualization only
NP_matrix = zeros(480,4);
NP_matrix(sub2ind(size(NP_matrix),Probe.row,Probe.col)) = 1;
Probe_tip = NaN(15,4);  Probe_tip(1:5,4) = 0;   Probe_tip(6:10,3:4) = 0;    Probe_tip(11:15,2:4) = 0;   
NP_matrix = [Probe_tip;NP_matrix];

% %=== Visualize probe (with tip for clarity)
% figure('units','normalized','outerposition',[.3 .3 .1 .6]);
% imagesc(NP_matrix,'AlphaData',~isnan(NP_matrix));
% set(gca,'YTickLabel',[],'XTickLabel',[],'YDir','normal'); axis('off');
% sgtitle(['Probe ',num2str(probe_number)]);

%=== Select channels for processing SWR (when 2 are collinear, keep the best RMS)
Probe.used = logical(Probe.idx*0);          
for ii = unique(Probe.row)'
    if sum(Probe.row == ii)==2
        collinear = find(Probe.row == ii);
        if Probe.rms(collinear(1))>Probe.rms(collinear(2))
            Probe.used(collinear(1)) = 1;
        else
            Probe.used(collinear(2)) = 1;
        end
    else
        Probe.used(Probe.row == ii) = 1;
    end
end
LFP_red = LFP(:,Probe.used);                % Select the subset of channels to be inspected
num_good_channels = size(LFP_red,2);        % Number of good channels
ids_good_channels = Probe.ch(Probe.used);   % Ids of the good channels
row_good_channels = Probe.row(Probe.used);  % Rows of the good channels

%=== Define offset vector for plotting LFPs
offset_vec = (row_good_channels-row_good_channels(1))'*scale_f/5;

%=== Populate the s cell with the single units and calculate firing rate
NP_unit = table();
NP_unit = [NP_unit; out.good_units];
NP_unit.fr = cellfun(@(x) 1/mean(diff(x/1e6)),NP_unit.spikeTimes_usec);
NP_unit = NP_unit(NP_unit.fr<1,:);      % Exclude neurons with high-firing rate
n_cells = size(NP_unit,1);              % Number of units
s = cell(n_cells,1);                    % Single units' spikes
for nc = 1:n_cells
    s{nc,1} = NP_unit.spikeTimes_usec{nc,1}/1e6;                    % Each row contains the time (in s) of the spikes from a single unit
    too_close = find(diff(s{nc,1})<0.001);                          % Clean up duplicated spikes
    if ~isempty(too_close)
        s{nc,1}(too_close+1) = [];
    end
end

%%  Process LFP to detect SWRs

%=== Process LFP to generate Ripple power signal
disp('Detecting Ripples');
fly_vector = logical(interp1(NP_imu.t,double(wBeats),t));       % Vector of 1s (flight) and 0s (rest), sampled at the LFP sampling frequency 
ripple_signal = bandpass(LFP_red,Ripple_band,Fs);               % Bandpass signal at the Ripple Band    
ripple_power  = zscore(abs(hilbert(ripple_signal)),[],1);       % Calculate z-scored power as the magnitude of the Hilbert transform
ripple_power  = smoothdata(ripple_power,1,'gaussian',0.05*Fs);  % Smooth with a 50 ms Gaussian kernel
ripple_power(fly_vector,:) = 0;                                 % Ripple Power during fligth periods is forced to 0
ripple_power(1:Fs*1,:) = 0; ripple_power(end+[-Fs*1:0],:) = 0;  % Ripple Power at the tails of the session if forced to 0

%===Detect candidate ripples for each channel
candidate_RPL = [];
for i=1:num_good_channels
    [RPL_a,RPL_t,RPL_d] = findpeaks(ripple_power(:,i),t,'MinPeakHeight',Ripple_th,'MinPeakDistance',0.05,'MinPeakWidth',0.01);
    candidate_RPL = [candidate_RPL;[RPL_t',RPL_a,ids_good_channels(i)*ones(numel(RPL_t),1),RPL_d']];       % Accumulate events
end

%=== Convert the matrix to a table with custom column names
candidate_RPL = array2table(candidate_RPL, 'VariableNames', {'t', 'amp', 'ch','dur'});

%=== Sort by time, find burst initiators and keep largest amplitude events
RPL = table();
candidate_RPL = sortrows(candidate_RPL,'t','ascend');
candidate_RPL.brst = [0; diff(candidate_RPL.t)<0.05];
initiators = find(~candidate_RPL.brst);
for ii = 1:numel(initiators)-1
    subtable = candidate_RPL(initiators(ii):initiators(ii+1)-1,:);
    [~,tmp_idx] = max(subtable.amp);
    RPL = [RPL; subtable(tmp_idx,:)];
end
RPL = join(RPL, Probe(:,{'ch','row'}), 'Keys', 'ch');

%=== Count the ripples at each depth
for i=1:N_channels
    if Probe.used(i)
        Probe.n_rpl(i)  = sum(RPL.row == Probe.row(i));
    else
        Probe.n_rpl(i) = NaN;
    end
end

%=== Store all the putative ripples in a 3D-matrix
All_RPL_waveforms = zeros(numel(int50ms),size(LFP_red,2),size(RPL,1));
for i=1:size(RPL,1)
    smp = 1+round((RPL.t(i,1)-t(1))*Fs);
    interval = smp+int50ms;
    All_RPL_waveforms(:,:,i) = LFP_red(interval,:);
end

%=== Calculate correlation between putative ripples and average
reshaped_All_RPL_waveforms = reshape(All_RPL_waveforms, [numel(int50ms)*num_good_channels,size(RPL,1)]);
RPL.corr = corr(reshaped_All_RPL_waveforms,mean(reshaped_All_RPL_waveforms,2));

%% Plot some results

%=== For each depth look at the amplitude and the correlation value with the average SWR
RPL_summary = varfun(@mean,RPL,'GroupingVariables','row');
figure('units','normalized','outerposition',[.3 .2 .1 .3]);
tiledlayout(1,3,'TileSpacing','compact');
nexttile;   plot(normalize(RPL_summary.GroupCount,'range'),RPL_summary.row,'k.');   xlabel('Count');                        ylim('tight');   
nexttile;   plot(normalize(RPL_summary.mean_amp,'range'),RPL_summary.row,'r.');     xlabel('Amplitude');     yticks([]);    ylim('tight'); 
nexttile;   plot(normalize(RPL_summary.mean_corr,'range'),RPL_summary.row,'b.');    xlabel('Corr w Avg');    yticks([]);    ylim('tight');
sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3},', Probe ',num2str(probe_number)],'Interpreter','none');
fig_count = saveFig(analysis_directory,[cell2mat(unique_ID),'RPL'],fig_count,options.savefigures);

%=== Show average Ripple profile
figure('units','normalized','outerposition',[.3 .2 .1 .5]);
sgtitle('Pick up the putative CA1 channel (just above)');
Ripple_avg = zeros(numel(int100ms),size(LFP_red,2));
for i=1:size(RPL,1)
        smp = 1+round((RPL.t(i)-t(1))*Fs);
        interval = smp+int100ms;
        Ripple_avg = Ripple_avg+LFP_red(interval,:);
end
Ripple_avg = Ripple_avg./size(RPL,1);       
plot(t(interval)-t(smp),Ripple_avg+offset_vec/4,'k','LineWidth',2); hold on;
xlim('tight');  ylim('tight');  plot([0 0],ylim,'r--');         
plot(-mean(RPL.dur)*[.5 .5],ylim,'b--');    plot(+mean(RPL.dur)*[.5 .5],ylim,'b--');    
plot(-mean(RPL.dur)*[1 1],ylim,'b--');    plot(+mean(RPL.dur)*[1 1],ylim,'b--');    
hold off;
xlabel('Time (s)'); 
sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3},', Probe ',num2str(probe_number)],'Interpreter','none');
fig_count = saveFig(analysis_directory,[cell2mat(unique_ID),'RPL'],fig_count,options.savefigures);

%=== Plot the best ripples (correlations)
[~,best_rpl] = sort(RPL.corr,'descend');
figure('units','normalized','outerposition',[0 0 1 1]);
tiledlayout(4,16,'TileSpacing','compact');
Conc_Ripples = [];
for i=1:min(size(RPL,1),4*16)
    smp = 1+round((RPL.t(best_rpl(i),1)-t(1))*Fs);
    interval = smp+int100ms;
    nexttile;    plot(t(interval)-t(smp),LFP_red(interval,:)+offset_vec,'k'); hold on;
    rectangle('Position',[-RPL.dur(best_rpl(i))*0.5 0 RPL.dur(best_rpl(i)) offset_vec(end)],'FaceColor',[1 0 0 0.1],'EdgeColor','none'); 
    xlim('tight');  ylim('tight');  plot([0 0],ylim,'r--');  hold off;
    yticks([]);  xlabel('Time (s)');
    if i==1, ylabel('Sorted by correlation with mean'); end
    Conc_Ripples = [Conc_Ripples ;LFP_red(interval,:)];
end
sgtitle('Best Example Ripples');
sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3},', Probe ',num2str(probe_number)],'Interpreter','none');
fig_count = saveFig(analysis_directory,[cell2mat(unique_ID),'RPL'],fig_count,options.savefigures);

%=== Plot the best ripples (amplitude)
[~,best_rpl] = sort(RPL.amp,'descend');
figure('units','normalized','outerposition',[0 0 1 1]);
tiledlayout(4,16,'TileSpacing','compact');
for i=1:min(size(RPL,1),4*16)
    smp = 1+round((RPL.t(best_rpl(i),1)-t(1))*Fs);
    interval = smp+int100ms;
    nexttile;    plot(t(interval)-t(smp),LFP_red(interval,:)+offset_vec,'k'); hold on;
    rectangle('Position',[-RPL.dur(best_rpl(i))*0.5 0 RPL.dur(best_rpl(i)) offset_vec(end)],'FaceColor',[1 0 0 0.1],'EdgeColor','none'); 
    xlim('tight');  ylim('tight');  plot([0 0],ylim,'r--');  hold off;
    yticks([]);  xlabel('Time (s)');
    if i==1, ylabel('Sorted by amplitude'); end
end
sgtitle('Best Example Ripples');
sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3},', Probe ',num2str(probe_number)],'Interpreter','none');
fig_count = saveFig(analysis_directory,[cell2mat(unique_ID),'RPL'],fig_count,options.savefigures);

%=== Ripple triggered PSTH (triggered by high quality ripples)
all_s = sort(vertcat(s{:}));
figure('units','normalized','outerposition',[.3 .1 .1 .4]);
Raster_AF_v3(all_s,RPL.t(RPL.corr>0.3),[],[],[-0.25 0.25],[],1,1);
sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3},', Probe ',num2str(probe_number)],'Interpreter','none');
fig_count = saveFig(analysis_directory,[cell2mat(unique_ID),'RPL'],fig_count,options.savefigures);

%=== Find the best channel from the best ripples
[psd,f_psd] = pwelch(Conc_Ripples,[],[],[],Fs);  psd = psd./sum(psd,1);
[~,best_red] = max(sum(psd(f_psd>Ripple_band(1) & f_psd<Ripple_band(2),:)));
best_ch = ids_good_channels(best_red);
inf_red = max(1,best_red-0);
sup_red = min(best_red+0,num_good_channels);
LFP_depth_prof = LFP_red(:,inf_red:sup_red);

%=== Plot the LFP trace and the detected ripples (only good)
ripple_logical = zeros(size(t));
for i=find(RPL.corr>0.1)'
    smp1 = 1+round((RPL.t(i)-RPL.dur(i)/2-t(1))*Fs);
    smp2 = 1+round((RPL.t(i)+RPL.dur(i)/2-t(1))*Fs);
    ripple_logical(smp1:smp2) = 1;
end
figure('units','normalized','outerposition',[0 .3 1 .4]);
plot(t,LFP_depth_prof+offset_vec(inf_red:sup_red),'k','LineWidth',0.7);       xlim('tight'); plot_lim = ylim;   
hold on;  area(t,ripple_logical*plot_lim(2),'FaceAlpha',0.5,'LineStyle','none','FaceColor',[0 0 1]);    ylim(plot_lim);
sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3},', Probe ',num2str(probe_number)],'Interpreter','none');
fig_count = saveFig(analysis_directory,[cell2mat(unique_ID),'RPL'],fig_count,options.savefigures);

%% Save data

%=== Create a structure for saving the relevant variables
RPL_out.probe = probe_number;                                  % Probe number
RPL_out.Fs = Fs;                                               % Sampling Frequency
RPL_out.th = Ripple_th;                                        % z-score threshold for detection
RPL_out.table = RPL;                                           % Table with all detected ripples
RPL_out.probe_info = Probe;                                    % Table with the Probe infos
RPL_out.template = Ripple_avg;                                 % Average Ripple Template
RPL_out.waveforms = All_RPL_waveforms;                         % All Ripple waveforms
RPL_out.LFP_trace = LFP_red(:,best_red);                       % LFP trace from the best channel
RPL_out.ch = best_ch;                                          % Channel corresponding to putative CA1 pyramidal layer
RPL_out.time_vector = t;                                       % Time vector for plotting the LFP

save([cd,'\RPL_probe',num2str(probe_number),'.mat'],'RPL_out');

end