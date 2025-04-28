%% Script for the analysis of LFP recorded with Neuropixels1.0
% First draft by A.F. on March 2024

%% Load data and extract relevant variables

probe_number = 2;

%=== Load data
load(['LFP_probe',num2str(probe_number),'.mat']);     % LFP data, saved by loadTrodesLFP_AF_v0
load('IMU_data.mat');                                 % IMU data, saved by loadTrodesAnalog_AF_v0

%=== Parameters
Fs = 1e3;                               % Sampling frequency for working on the LFP  
Fs_imu = 1/mean(diff(NP_imu.t));        % Sampling frequency of the IMU data
Ripple_band = [100 200];                % Ripples frequency band (From MY: [80 160])
Ripple_th = 3;                          % Threshold on zscored Ripple power (3.5)
Ripple_th_ave = 5;                      % Threshold on zscored Ripple power for calculating the average template
use_CA1_bank = 0;                       % If using only the upper bank of channels to calculate the correlation with the template

%=== Transform the LFP into double and interpolate at 1kHz even sampling
t = linspace(red_out.t_ds(1),red_out.t_ds(end),(red_out.t_ds(end)-red_out.t_ds(1))*Fs);
LFP = interp1(red_out.t_ds,double(red_out.lfp).*red_out.voltage_scaling,t);
N_channels = size(LFP,2);

%=== Calculate and normalize the PSDs of the LFP
[psd,f_psd] = pwelch(LFP,[],[],[],Fs);  psd = psd./sum(psd,1);

%% Process IMU data to extract flights, look at average LFP

%=== Extract flight periods from accelerometer
a_abs_NP = vecnorm(NP_imu.acc,2,2);                         % Absolute acceleration
a_flt_NP = bandpass(a_abs_NP,[7 9],Fs_imu);                 % Filtered at the wing-beat frequency
[up,lo] = envelope(a_flt_NP,round(0.06*Fs_imu),'peak');     % Upper and lower envelopes 
env = normalize(up - lo,'range');                           % Amplitude of the envelope
env_th = otsuthresh(histcounts(env));                       % Threshold (based on Otsu method). Can be set at 0.35
wBeats = movsum(env>env_th,2*Fs_imu)>Fs_imu/5;              % Euristic criterion for flight detection
LFP_interp = interp1(t,LFP,NP_imu.t,'linear','extrap');     % Interpolate LFP at the accelerometer time samples (useful later)

%=== Plot wing-beats and average LFP
figure('units','normalized','outerposition',[0 .2 1 .5]);
tiledlayout(2,1,'TileSpacing','compact');
ax(1) = nexttile;       area(NP_imu.t,wBeats*1,'FaceAlpha',0.3,'LineStyle','none');  hold on;
plot(NP_imu.t,normalize(a_flt_NP,'range',[-1 1]));  plot(NP_imu.t,normalize(env,'range'));  hold off;   ylabel('Accelerometer');        xticks([]);
ax(2) = nexttile;    area(NP_imu.t,wBeats*1,'FaceAlpha',0.3,'LineStyle','none');  hold on;
plot(t,normalize(mean(LFP,2),'range'));                                                                 ylabel('Average LFP (norm)');   xticks([]);
linkaxes(ax,'x');   xlim('tight');  xlabel('Time (s)');

%=== Calculate the cross-correlation between ave LFP and Accelerometer
figure('units','normalized','outerposition',[.1 .2 .2 .3]);
[acc_r,acc_lags] = xcorr(mean(LFP_interp,2),a_abs_NP,round(0.3*Fs_imu));
plot(acc_lags/Fs_imu,acc_r);    hold on;    plot([0 0],ylim);   xlabel('Lag(s)');   ylabel('Cross-correlation');
[~,max_corr_loc] = max(acc_r);  title(['Delay: ',num2str(acc_lags(max_corr_loc)/Fs,1),'ms']);

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
t_10Hz = [t(1):0.1:t(end)];
LFP_2look = interp1(t,LFP,t_10Hz,'linear','extrap');
scale_f = 0.5*median(range(LFP,1));
figure('units','normalized','outerposition',[0 0 .7 1]);
plot(t_10Hz,LFP_2look+flip([1:N_channels]*scale_f),'k');    hold on;
area(t_10Hz,logical(interp1(NP_imu.t,double(wBeats),t_10Hz))*N_channels*scale_f,'FaceAlpha',0.5,'LineStyle','none','FaceColor',[0 0 1]);
xlim('tight');  ylim('tight');  xlabel('Time (s)'); ylabel('LFP');
Y_labels = num2str(flip([row_idx,col_idx,red_out.channelMap'],1));                    
yticks([1:N_channels]*scale_f);   yticklabels(Y_labels);

%% Ripple detection

%=== Keep only one of the two collinear channels (the one with higher RMS)
tmp_diff = diff(rms(LFP,1));    good_ch = [1:2:N_channels]+(tmp_diff(1:2:end)>0);
LFP_red = LFP(:,good_ch);       % [row_idx(good_ch),col_idx(good_ch)]

%=== Process LFP to generate Ripple power signal
disp('Detecting Ripples');
fly_vector = logical(interp1(NP_imu.t,double(wBeats),t));       % Vector of 1s (flight) and 0s (rest), sampled at the LFP sampling frequency 
ripple_signal = bandpass(LFP_red,Ripple_band,Fs);               % Bandpass signal at the Ripple Band    
ripple_power  = zscore(abs(hilbert(ripple_signal)),[],1);       % Calculate z-scored power as the magnitude of the Hilbert transform
ripple_power  = smoothdata(ripple_power,1,'gaussian',0.05*Fs);  % Smooth with a 50 ms Gaussian kernel
ripple_power(fly_vector,:) = 0;                                 % Ripple Power during fligth periods is forced to 0

%===Detect candidate ripples for each channel
candidate_RPL = [];
for i=1:numel(good_ch)
    
    [RPL_a,RPL_t,RPL_d] = findpeaks(ripple_power(:,i),t,'MinPeakHeight',Ripple_th,'MinPeakDistance',0.1,'MinPeakWidth',0.01);
    
    %===Clean up times too close to start and stop of the recording
    valid = RPL_t>1 & RPL_t<(t(end)-1);
    RPL_a = RPL_a(valid);
    RPL_t = RPL_t(valid);
    RPL_d = RPL_d(valid);
    
    candidate_RPL = [candidate_RPL;[RPL_t',RPL_a,i*ones(numel(RPL_t),1),RPL_d']];       % Accumulate events
end

%=== Eliminate temporally-adjacent events across channels
candidate_RPL = sortrows(candidate_RPL);
parsed_seq = parse_into_bursts_AF_v0(candidate_RPL(:,1),0.10);                                            % Parse candidate events into bursts
parsed_amplitudes = mat2cell(candidate_RPL(:,2),cellfun(@(x) size(x,1),parsed_seq),1);                    % Parse the corresponding amplitudes
parsed_channels = mat2cell(candidate_RPL(:,3),cellfun(@(x) size(x,1),parsed_seq),1);                      % Parse the corresponding channels
parsed_durations = mat2cell(candidate_RPL(:,4),cellfun(@(x) size(x,1),parsed_seq),1);                     % Parse the corresponding durations
RPL_cell = cellfun(@(x,y,z,d) [x,y,z,d],parsed_seq,parsed_amplitudes,parsed_channels,parsed_durations,'UniformOutput',false);  % Generate a cell containing all the information
put_RPL_t = cellfun(@(x,y) x(y,1),RPL_cell,num2cell(cellfun(@(x) find(x(:,3)== max(x(:,3))),RPL_cell)));  % Within a burst, keep the event associated with the largest amplitude
put_RPL_a = cellfun(@(x,y) x(y,3),RPL_cell,num2cell(cellfun(@(x) find(x(:,3)== max(x(:,3))),RPL_cell)));  % Store the amplitudes of the events
put_RPL_c = cellfun(@(x,y) x(y,4),RPL_cell,num2cell(cellfun(@(x) find(x(:,3)== max(x(:,3))),RPL_cell)));  % Store the channels of the events
put_RPL_d = cellfun(@(x,y) x(y,5),RPL_cell,num2cell(cellfun(@(x) find(x(:,3)== max(x(:,3))),RPL_cell)));  % Store the durations of the events

%=== Plot the time differences between putative 'same ripples'

%=== Plot the distribution of channels where the putative ripples were found
%=== and the time difference between pooled same-ripples from different channels
figure('units','normalized','outerposition',[.3 .2 .3 .7]);
tiledlayout(2,2,'TileSpacing','compact');
rpl_counts = histcounts(candidate_RPL(:,3),0.5+[0:N_channels/2]);
nexttile(1,[2,1]);   barh(red_out.channelPositions(good_ch,2),rpl_counts);   set(gca, 'YDir', 'reverse');    
ylabel('Depth');    xlabel('Count');    title('Putative Ripples');
tmp_isi = cellfun(@(x) diff(x(:,1)),parsed_seq,'UniformOutput',false);
nexttile;   histogram(vertcat(tmp_isi{:}),[0:1e-3:0.1]);    xlabel('Time (s)'); ylabel('Counts');   title('ISI between same-ripples');    
nexttile;   histogram(diff(put_RPL_t),[0:1e-2:1]);          xlabel('Time (s)'); ylabel('Counts');   title('Ripple ISI');

%=== Store all the putative ripples in a 3D-matrix
int50ms = [round(Fs*-0.05):round(Fs*0.05)];
int100ms = [round(Fs*-0.1):round(Fs*0.1)];
All_RPL_waveforms = zeros(numel(int50ms),size(LFP_red,2),size(put_RPL_t,1));
for i=1:size(put_RPL_t,1)
    smp = 1+round((put_RPL_t(i,1)-t(1))*Fs);
    interval = smp+int50ms;
    All_RPL_waveforms(:,:,i) = LFP_red(interval,:);
end

%=== Calculate correlation between putative ripples and average
reshaped_All_RPL_waveforms = reshape(All_RPL_waveforms, [numel(int50ms)*N_channels/2,size(put_RPL_t,1)]);
CA1_bank = numel(int50ms)*find(diff(red_out.channelID)>1)/2;
if use_CA1_bank
    put_RPL_corr = corr(reshaped_All_RPL_waveforms(CA1_bank,:),mean(reshaped_All_RPL_waveforms(CA1_bank,put_RPL_a>=Ripple_th_ave),2));
else
    put_RPL_corr = corr(reshaped_All_RPL_waveforms,mean(reshaped_All_RPL_waveforms(:,put_RPL_a>=Ripple_th_ave),2));
end
[~,best_rpl] = sort(put_RPL_corr,'descend');

%=== Plot the average ripple waveform across channels and the distribution of correlation values
%=== Manually select the channel corresponding to just above the putative CA1 pyramidal layer
figure('units','normalized','outerposition',[.3 .2 .4 .7]);
tiledlayout(1,2,'TileSpacing','compact');
sgtitle('Pick up the putative CA1 channel (just above)');
Ripple_avg = zeros(numel(int100ms),size(LFP_red,2));
for i=1:size(put_RPL_t,1)
    if put_RPL_a(i)>=Ripple_th_ave
        smp = 1+round((put_RPL_t(i,1)-t(1))*Fs);
        interval = smp+int100ms;
        Ripple_avg = Ripple_avg+LFP_red(interval,:);
    end
end
Ripple_avg = Ripple_avg./size(put_RPL_t,1);       
offset_vec = ((red_out.channelPositions(good_ch,2)-red_out.channelPositions(1,2))./20)';
%nexttile;   plot(t(interval)-t(smp),1+(-Ripple_avg+Ripple_avg(1,:))/std(Ripple_avg,[],'all')+offset_vec,'k','LineWidth',2); hold on;
nexttile;   plot(t(interval)-t(smp),1+(+Ripple_avg-Ripple_avg(1,:))/std(Ripple_avg,[],'all')+offset_vec,'k','LineWidth',2); hold on;
xlim('tight');  ylim('tight');  plot([0 0],ylim,'r--');         
plot(-mean(put_RPL_d)*[.5 .5],ylim,'b--');    plot(+mean(put_RPL_d)*[.5 .5],ylim,'b--');    
plot(-mean(put_RPL_d)*[1 1],ylim,'b--');    plot(+mean(put_RPL_d)*[1 1],ylim,'b--');    
hold off;
xlabel('Time (s)'); title(['Average of ', num2str(sum(put_RPL_a>=Ripple_th_ave)) ,' events']);
set(gca, 'YDir','reverse');
[~,best_ch_RPL] = ginput(1);    best_ch_RPL = knnsearch(offset_vec',best_ch_RPL);
nexttile;   barh(sort(put_RPL_corr,'descend'));    hold on;    plot(.3*[1 1],ylim);    
xlabel('Correlation with template');    ylabel('Count');
sgtitle(['Putative CA1 channel (just above):',num2str(best_ch_RPL)]);

%=== Plot the LFP trace and the detected ripples (only good)
ripple_logical = zeros(size(t));
for i=find(put_RPL_c>0.6 & put_RPL_a>5)'
    smp1 = 1+round((put_RPL_t(i)-put_RPL_d(i)/2-t(1))*Fs);
    smp2 = 1+round((put_RPL_t(i)+put_RPL_d(i)/2-t(1))*Fs);
    ripple_logical(smp1:smp2) = 1;
end
figure('units','normalized','outerposition',[0 .3 1 .5]);
LFP_depth_prof = -LFP_red(:,best_ch_RPL+[-5:5]);
offset_vec = flip([1:size(LFP_depth_prof,2)])*0.07;
plot(t,normalize(LFP_depth_prof,1,'range')+offset_vec,'k','LineWidth',0.7);       xlim('tight');   
hold on;  area(t,ripple_logical*100,'FaceAlpha',0.5,'LineStyle','none','FaceColor',[1 0 1]);    ylim([0 2]);       

%=== Plot the best ripples
figure('units','normalized','outerposition',[0 0 1 1]);
tiledlayout(4,16,'TileSpacing','compact');
for i=1:4*16
    %smp = 1+round((put_RPL(datasample(1:size(put_RPL,1),1),1)-t(1))*Fs);
    smp = 1+round((put_RPL_t(best_rpl(i),1)-t(1))*Fs);
    interval = smp+int100ms;
    offset_vec = ((red_out.channelPositions(good_ch,2)-red_out.channelPositions(1,2))./20)'*scale_f/5;
    nexttile;    plot(t(interval)-t(smp),LFP_red(interval,:)+offset_vec,'k'); hold on;
    %plot(-put_RPL_d(best_rpl(i))*[.5 .5],ylim,'b--');    plot(+put_RPL_d(best_rpl(i))*[.5 .5],ylim,'b--'); 
    rectangle('Position',[-put_RPL_d(best_rpl(i))*0.5 0 put_RPL_d(best_rpl(i)) offset_vec(end)],'FaceColor',[1 0 0 0.1],'EdgeColor','none'); 
    xlim('tight');  ylim('tight');  plot([0 0],ylim,'r--');  hold off;
    yticks([]); set(gca, 'YDir','reverse'); 
    xlabel('Time (s)'); 
end
sgtitle('Best Example Ripples');

%=== Plot some statistics of the events
corr_th = 0.3;
figure('units','normalized','outerposition',[.3 .2 .3 .7]);
tiledlayout(2,2,'TileSpacing','compact');
rpl_counts = histcounts(put_RPL_c(put_RPL_corr>corr_th),0.5+[0:N_channels/2]);
nexttile(1,[2,1]);   barh(red_out.channelPositions(good_ch,2),rpl_counts);   set(gca, 'YDir', 'reverse');    
ylabel('Depth');    xlabel('Count');    title('Putative Ripples');
tmp_isi = cellfun(@(x) diff(x(:,1)),parsed_seq,'UniformOutput',false);
nexttile;   histogram(put_RPL_d(put_RPL_corr>corr_th),[0:1e-2:0.2]);        xlabel('Time (s)'); ylabel('Counts');   title('Ripple Duration');    
nexttile;   histogram(diff(put_RPL_t(put_RPL_corr>corr_th)),[0:1e-2:1]);    xlabel('Time (s)'); ylabel('Counts');   title('Ripple ISI');
sgtitle('Good Ripples');


%% Save the detected Ripples 

%=== Create a structure for saving the relevant variables
RPL.probe = probe_number;                                  % Probe number
RPL.Fs = Fs;                                               % Sampling Frequency
RPL.th = Ripple_th;                                        % z-score threshold for detection
RPL.th_ave = Ripple_th_ave;                                % z-score threshold for calculating the average
RPL.t = put_RPL_t;                                         % Time of the Ripple
RPL.ch = put_RPL_c;                                        % Channel with maximal amplitude
RPL.a = put_RPL_a;                                         % Amplitude of the event
RPL.dur = put_RPL_d;                                       % Duration of the event
RPL.corr = put_RPL_corr;                                   % Correlation with template
RPL.ch_geometry = red_out.channelPositions(good_ch,:);     % Horz and Vert position of each channel
RPL.template = Ripple_avg;                                 % Average Ripple Template
RPL.waveforms = All_RPL_waveforms;                         % All Ripple waveforms
RPL.CA1_ch_ds = best_ch_RPL;                               % Channel (downsampled) corresponding to just above the putative CA1 pyramidal layer
RPL.CA1_ch_real = good_ch(best_ch_RPL);                    % Channel (real) corresponding to just above the putative CA1 pyramidal layer                                       
RPL.CA1_LFP = -LFP_red(:,best_ch_RPL);                     % LFP from the putative CA1 pyramidal layer
RPL.CA1_LFP_nb = -LFP_red(:,best_ch_RPL-5:best_ch_RPL+5);  % LFP from the putative CA1 pyramidal layer (across depth)
RPL.time_vector = t;                                       % Time vector for plotting the LFP

save([cd,'\Ripples_probe',num2str(probe_number),'.mat'],'RPL');
