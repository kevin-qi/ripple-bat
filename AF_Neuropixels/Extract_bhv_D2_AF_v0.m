%% Brief Function to Extract positional data for Ciholas recordings at the Field Station 
clr;

%=== Load data and session name
extracted_CDPfile = dir(fullfile(cd, '*extracted_*cdp_1*'));    load(extracted_CDPfile.name);    
bat_filename = strsplit(extracted_CDPfile.name,'_');            bat_filename = bat_filename{1,2};
batdate = bat_filename(1:6);

%===Room dimensions and params
r_lim = [-4.6 5.5; -2.3 1.8; 0 2.8]*1.1;           % Room boundaries (Field Station)
Fs = 120;                                          % Acquisition Frame Rate (Original was 100 Hz)
v_th = 0.5;                                        % Velocity threshold (m/s) for flight segmentation
tag_id = 2;                                        % Identity of the tag associated to the bat

options.use_r_corr = 1;                            % If using r corrected
options.unique_ID = {'Dataset_2','14445',batdate}; % Save the session identifier
options.savedata = 1;                              % If saving the data
options.savefigures = 1;                           % Save Figures
fig_count = 1;                                     % Save Figures

%=== Analysis folder for storing the results
if options.savedata
    analysis_directory=fullfile(pwd,['Ext_Behavior_',datestr(now, 'yymmdd_HHMM')]);
    if ~exist(analysis_directory,'dir')
        mkdir(analysis_directory);
    end
end

%% Calculate evenly sampled kinematic variables (r,v,a) and wing-beats epochs

%=== Remove duplicate samples 
[~,ia,~] = unique(tag_data{1,tag_id}(:,8),'stable');         tag_data{1,tag_id} = tag_data{1,tag_id}(ia,:);              
[~,ia,~] = unique(tag_data_filt{1,tag_id}(:,8),'stable');    tag_data_filt{1,tag_id} = tag_data_filt{1,tag_id}(ia,:);    
[~,ia,~] = unique(tag_ac_data{1,tag_id}(:,8),'stable');      tag_ac_data{1,tag_id} = tag_ac_data{1,tag_id}(ia,:); 

%=== Define t vector
t = [CDPmtdata.TTL_times(1):1/Fs:CDPmtdata.TTL_times(end)]';       %Evenly sampled time vector from first to last TTL.
T = length(t);                                                     %Number of time samples

%=== Interpolate position at evenly spaced time points
r =  interp1(tag_data_filt{1,tag_id}(:,8), tag_data_filt{1,tag_id}(:,[3:5]),t,'linear','extrap'); %DO NOT use spline interpolation!

%=== Ad hoc corrections
if strcmp(extracted_CDPfile.name  ,'extracted_231210_cdp_2.mat')
    for i=1:3
        r(:,i) = interp1(tag_data_filt{1,tag_id}(:,8), tag_data_filt{1,tag_id}(:,2+i),t,'linear',tag_data_filt{1,tag_id}(1,2+i));
    end
end

%=== Interpolate acceleration at evenly spaced time points
a =  interp1(tag_ac_data{1,tag_id}(:,8), tag_ac_data{1,tag_id}(:,[3:5]),t,'linear','extrap'); %DO NOT use spline interpolation!

%=== Ad hoc corrections
if strcmp(extracted_CDPfile.name  ,'extracted_231210_cdp_2.mat')
    for i=1:3
        a(:,i) = interp1(tag_ac_data{1,tag_id}(:,8), tag_ac_data{1,tag_id}(:,2+i),t,'linear',tag_ac_data{1,tag_id}(1,2+i));
    end
end

a_abs = squeeze(vecnorm(a,2,2));    %Modulus
a_flt = bandpass(a_abs,[7 9],100);  %Filtered at the wing-beat frequency

%=== Calculate velocity and 2D-direction of motion
v = diff(r,1,1)./diff(t); v=cat(1,zeros(1,3),v);   v_abs = squeeze(vecnorm(v,2,2));

%=== Detect flight epochs based on wing-beat signal
[up,lo] = envelope(a_flt,10,'peak');                                        %Envelope of the acceleration signal (noise addedd to avoid problems with splining)
env = normalize(up - lo,'range');                                           %Amplitude of the envelope
env_th = otsuthresh(histcounts(env));                                       %Threshold (based on Otsu method). Can be set at 0.35
wBeats = movsum(env>env_th,2*Fs)>Fs/5;                                      %Euristic criterion for flight detection
   
%=== Sanity check on sampling intervals
figure('units','normalized','outerposition',[.2 .4 .15 .3]);
histogram(diff(t)*1e3); hold on;    plot(1e3/Fs*[1 1],ylim);  hold off;   xlim(10*[-1,1]+1e3/Fs); xlabel('Time between samples (ms)');  ylabel('Fraction');
sgtitle(['Session ', batdate],'Interpreter','none');
fig_count = saveFig(analysis_directory,batdate,fig_count,options.savefigures);

% % Uncomment to control 
% area(t,wBeats*1,'FaceAlpha',0.3,'LineStyle','none');  hold on;
% plot(t,normalize(v_abs,'range'));
% plot(t,r(:,1),t,r(:,2));
% plot(t,normalize(a_flt,'range',[-1 1]));
% plot(t,normalize(movsum(env>env_th,2*Fs),'range'));
% hold off;

% %=== Check interpolation quality, Uncomment to control
% figure('units','normalized','outerposition',[0 0 1 1]);
% tiledlayout(1,3,'TileSpacing','tight');
% ax(1) = subplot(311);  plot(tag_data{1,tag_id}(:,8),tag_data{1,tag_id}(:,3)',t,r(:,1));    legend('Raw','Processed')
% ax(2) = subplot(312);  plot(tag_data{1,tag_id}(:,8),tag_data{1,tag_id}(:,4)',t,r(:,2));
% ax(3) = subplot(313);  plot(tag_data{1,tag_id}(:,8),tag_data{1,tag_id}(:,5)',t,r(:,3));
% linkaxes(ax,'x');
% sgtitle(['Session ', batdate],'Interpreter','none');

%% Correct position by median filtering when the bat is not flapping its wings

%=== Find stationary epochs
stat_periods = repmat(~wBeats,1,1,3);    stat_periods = permute(stat_periods,[1 3 2]);
r_stat = r;    r_stat(~stat_periods) = nan;

%=== Filter position during stationary epochs
r_stat =  smoothdata(r_stat,1,'movmedian',Fs*5,'omitnan');
r_stat(~stat_periods) = nan;

%=== Substitute median filtered data when bat is not flapping its wings
r_stat = fillmissing(r_stat,'constant',0,1);
r_corr = r_stat.*stat_periods+r.*(~stat_periods);

%=== Connect segments of filtered and non-filtered data
r_corr = smoothdata(r_corr,1,'loess',round(Fs*0.7));     %DO NOT USE lowess!!
%r_corr = smoothdata(r_corr,1,'movmedian',round(Fs*1));     %DO NOT USE lowess!!

%=== Check filtering quality
figure('units','normalized','outerposition',[0 0 1 1]);
tiledlayout(1,3,'TileSpacing','tight');
ax(1) = subplot(311);  plot(tag_data{1,tag_id}(:,8),tag_data{1,tag_id}(:,3)',t,r_corr(:,1));    legend('Raw','Processed')
ax(2) = subplot(312);  plot(tag_data{1,tag_id}(:,8),tag_data{1,tag_id}(:,4)',t,r_corr(:,2));
ax(3) = subplot(313);  plot(tag_data{1,tag_id}(:,8),tag_data{1,tag_id}(:,5)',t,r_corr(:,3));
linkaxes(ax,'x');
sgtitle(['Session ', batdate],'Interpreter','none');
fig_count = saveFig(analysis_directory,batdate,fig_count,options.savefigures);
  
%% Correct position, recalculate velocity and heading angle

if options.use_r_corr
    r_old = r;  v_old = v_abs;  %Store old versions, in case you need them
    r = r_corr; 
    v = diff(r,1,1)./diff(t); v=cat(1,zeros(1,3),v);   v_abs = squeeze(vecnorm(v,2,2)); 
end

%=== FIGURE: Raw Velocity VS Corrected Velocity
figure('units','normalized','outerposition',[0 .3 1 .3]);
plot(t,v_old,'.');    hold on;     sgtitle('v');
plot(t,v_abs);    hold off;    legend('raw','corrected');      ylabel('m/s');
xlabel('Time(s)'); xlim([0,t(T)]);
sgtitle(['Session ', batdate],'Interpreter','none');
fig_count = saveFig(analysis_directory,batdate,fig_count,options.savefigures);

%=== Visualize raw and processed data
figure('units','normalized','outerposition',[.2 .3 .5 .5]);
tiledlayout(1,3,'TileSpacing','tight');
sgtitle('C3D Raw Data');
nexttile;   plot(tag_data{1,tag_id}(:,3),tag_data{1,tag_id}(:,4),'.r');   hold on;  plot(r(:,1),r(:,2),'k-');   title('TOP');     axis equal; xticks([]); yticks([]); xlim(r_lim(1,:)*1.1); ylim(r_lim(2,:)*1.1);
nexttile;   plot(tag_data{1,tag_id}(:,3),tag_data{1,tag_id}(:,5),'.r');   hold on;  plot(r(:,1),r(:,3),'k-');   axis equal; xticks([]); yticks([]);   xlim(r_lim(1,:)*1.1); ylim(r_lim(3,:)*1.1); title('SIDE1');
nexttile;   plot(tag_data{1,tag_id}(:,4),tag_data{1,tag_id}(:,5),'.r');   hold on;  plot(r(:,2),r(:,3),'k-');   axis equal; xticks([]); yticks([]);   xlim(r_lim(2,:)*1.1); ylim(r_lim(3,:)*1.1); title('SIDE2');
sgtitle(['Session ', batdate],'Interpreter','none');
fig_count = saveFig(analysis_directory,batdate,fig_count,options.savefigures);

%% Flight segmentation and clustering

%=== Flight segmentation
[bflying,f_num,f_smp(:,1),f_smp(:,2)] = FlightSegm_AF_v0(v_abs,v_th,Fs);

%=== FIGURE: Velocity and flight segmentation
figure('units','normalized','outerposition',[0 0 1 .5]);
area(t,bflying*5,'FaceAlpha',0.3,'LineStyle','none');  hold on;
plot(t,v_abs,'r');     plot(t,r(:,1),'k--');    plot(t,r(:,2),'b--');  ylabel('Velocity (m/s)');     hold off;
legend('Fly','Vel','x(m)','y(m)'); %ylim([-3 6]);
title([num2str(f_num) ' flights']);
xlabel('Time(s)');  xlim('tight');
sgtitle(['Session ', batdate],'Interpreter','none');
fig_count = saveFig(analysis_directory,batdate,fig_count,options.savefigures);

%=== Cluster flights
alpha_clus = 1.3;         %Parameter for flight clustering
f_clus = FlightClus_AF_v3(r,bflying,Fs,'Alpha',alpha_clus,'Frechet',1,'Points',7);
sgtitle(['Session ', batdate],'Interpreter','none');
fig_count = saveFig(analysis_directory,batdate,fig_count,options.savefigures);

%% Save Relevant Data

% T:        total number of samples covered by the tracking system
% t:        time relative to the first M-9 TTL and sampled at Fs
% bflying:  T x 1 column vector with flight status (1 flying, 0 resting)
% f_num:    number of flights
% f_smp:    f_num x 2 matrix with the samples of takeoff and landing    
% Fs:       Sampling Frequency (default 120 Hz, inherited from Cortex)
% options:  options and metadata
% r:        T x 3 matrix with the average x,y,z position of the cap
% r_lim:    Room boundaries
% v:        T x 3 matrix with the average x,y,z velocity of the cap
% v_abs:    T x 1 matrix with the average absolute velocity of the cap
% v_th:     Velocity threshold used for flight segmentation

%=== Analysis folder for storing the results
if options.savedata
    save([analysis_directory,'/Extracted_Behavior_', batdate, '.mat'],'bflying','f_num','f_smp','Fs','options','r','r_lim','t','T','v','v_abs','v_th');
end
