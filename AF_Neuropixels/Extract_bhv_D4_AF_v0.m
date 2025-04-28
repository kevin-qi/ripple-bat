%% Brief Function to Extract c3d Marker data for Standard Cortex based recordings in the flight room 

%=== Load Data
clr;
c3d_clus_file = dir(fullfile(cd, '*_tracking_1.mat*'));                 c3d_data = load(c3d_clus_file.name);       % Cortex
bat_filename = strsplit(c3d_clus_file.name,'_');                    
batdate = bat_filename{1, 2};
batname = bat_filename{1, 1};
check_TTLs = 1;

%===Room dimensions and params
r_lim = [-2.9 2.9; -2.6 2.6; 0 2.30];              % Room boundaries
Fs = c3d_data.AnalogFrameRate;                     % Acquisition Frame Rate
v_th = 0.5;                                        % Velocity threshold (m/s) for flight segmentation
T = numel(c3d_data.AnalogSignals(:,3));            % Total number of c3d samples
correct_artifacts = 0;                             % If attempting to correct tracking artifacts
options.unique_ID = {'Dataset_4',batname,batdate}; % Save the session identifier
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

%% Extract position and time vectors

%=== Convert c3d position into m and set missed points to NaN
c3d_data.Markers = c3d_data.Markers(:,1:3,:)/1e3;   c3d_data.Markers(c3d_data.Markers==0) = nan;

%=== Ad hoc corrections
if isequal(options.unique_ID,[{'Dataset_4'},{'14640'},{'250324'}])
    artifact = c3d_data.AnalogSignals(:,end)<0;
    c3d_data.AnalogSignals(artifact,end) = 0;
    c3d_data.AnalogSignals(1:5e3,end) = 0;
end

%=== Set time 0 to the rising edge of the first TTL, then interpolate by assuming each TTL is separated by 3s
[~,smp_rising] = pulsewidth(c3d_data.AnalogSignals(:,3));

%=== Ad hoc corrections
if isequal(options.unique_ID,[{'Dataset_4'},{'14640'},{'250324'}]) || isequal(options.unique_ID,[{'Dataset_4'},{'14640'},{'250321'}])
    smp_rising = [smp_rising(1)-mean(diff(smp_rising)) ;   smp_rising];
end

TTL_times = linspace(0,3*(numel(smp_rising)-1),numel(smp_rising));
t = interp1(smp_rising,TTL_times,1:T,'linear','extrap');

if check_TTLs
    load('SU_kilosort4_outdir_probe1.mat')
    if numel(out.local_ttl_timestamps_usec)~=numel(smp_rising)
        disp('Something is wrong with TTLs');
        disp([num2str(numel(out.local_ttl_timestamps_usec)-numel(smp_rising)),' TTL lost by Cortex']);
    else
        disp('All good with TTLs');
    end
    
end

%=== Sanity check on Cortex sampling
figure('units','normalized','outerposition',[.2 .4 .3 .3]);
tiledlayout(1,2,'TileSpacing','tight');
nexttile;   histogram(diff(t)*1e3); hold on;    plot(1e3/Fs*[1 1],ylim);  hold off;   xlim(10*[-1,1]+1e3/Fs); xlabel('Time between samples (ms)');  ylabel('Fraction');
nexttile;   histogram(diff(smp_rising./c3d_data.VideoFrameRate));   xlabel('Time between detected TTLs (ms)');     
sgtitle(['Session ', batdate],'Interpreter','none');
fig_count = saveFig(analysis_directory,batdate,fig_count,options.savefigures);

%% Process Cortex Data

%=== Get the mean marker position
r_mn = squeeze(median(c3d_data.Markers,2,'omitnan'));

%=== Ad hoc corrections
if isequal(options.unique_ID,[{'Dataset_3'},{'14543'},{'240422'}])

    inv_x = r_mn(:,1)>1.65 & r_mn(:,1)<1.67;
    inv_y = r_mn(:,2)>1.33 & r_mn(:,2)<1.35;
    inv_z = r_mn(:,3)>1.74 & r_mn(:,3)<1.76;
    
    r_mn(inv_x & inv_y & inv_z,:) = NaN;
    
end


%=== Delete teleported data (too fast)
v_tlp_th = 6;
v_tlp = v_tlp_th+1;
while sum(v_tlp>v_tlp_th)
    cap_idx = ~isnan(r_mn(:,1));    nonnan_markers = find(cap_idx);
    v_tlp = vecnorm(diff(r_mn(cap_idx,:))./diff(t(cap_idx)'),2,2);
    r_mn(nonnan_markers(find(v_tlp>v_tlp_th)+1),:) = NaN(numel(find(v_tlp>v_tlp_th)),3);
end

%=== Fill in missing data and smooth
cap_idx = ~isnan(r_mn(:,1));
r_fl = interp1(t(cap_idx),r_mn(cap_idx,:),t,'previous','extrap');
r_ft = smoothdata(r_fl,1,'lowess',round(Fs*0.5));

%=== Define final vectors
r = r_ft;
v = diff(r,1,1)./diff(t)'; v=cat(1,zeros(1,3),v);   v_abs = vecnorm(v,2,2);

%=== Plot at different steps of the pre-processing
figure('units','normalized','outerposition',[0 0 1 1]);
tiledlayout(4,1,'TileSpacing','tight');
for i=1:3
    ax(i) = nexttile;   plot(t,r_mn(:,i),'r.');    hold on; plot(t,r_fl(:,i),'g-'); plot(t,r_ft(:,i),'k-'); ylim(r_lim(i,:));
    legend('Raw','Cleaned','Smoothed'); ylabel(['x',num2str(i),' (m)']);
end
ax(4) = nexttile;   plot(t,v_abs,'b'); ylabel('Velocity (m/s)');
linkaxes(ax,'x');   xlabel('Time (s)');  xlim('tight');
sgtitle(['Session ', batdate],'Interpreter','none');
fig_count = saveFig(analysis_directory,batdate,fig_count,options.savefigures);

%% Correct remaining tracking errors

if correct_artifacts
    
    %=== Flight segmentation
    f_smp = []; [bflying,f_num,f_smp(:,1),f_smp(:,2)] = FlightSegm_AF_v0(v_abs,v_th,Fs);
    
    %=== Find periods of abnormal acceleration
    a_diff = abs([0; diff(v_abs)]);
    bartifact = double(a_diff>0.12);
    bartifact = movmax(bartifact,round(Fs*0.1));
    
    %=== Force regions of high acceleration (and no flights) to NaN and fill missing entries
    r_corr = r;
    [~,start_artifact,stop_artifact] = pulsewidth(bartifact);
    for i=1:numel(start_artifact)-1
        candidate_interval = floor(start_artifact(i)):ceil(start_artifact(i+1));
        if ~any(bflying(candidate_interval)) || numel(candidate_interval)<Fs*1
            r_corr(floor(start_artifact(i)):ceil(stop_artifact(i+1)),:) = NaN;
        else
            r_corr(floor(start_artifact(i)):ceil(stop_artifact(i)),:) = NaN;
        end
    end
    r_corr = fillmissing(r_corr,'previous');
    
    %=== Plot at different steps of the pre-processing
    figure('units','normalized','outerposition',[0 0 1 1]);
    tiledlayout(4,1,'TileSpacing','tight');
    for i=1:3
        ax(i) = nexttile;   plot(t,r_corr(:,i),'r.');    hold on; plot(t,r(:,i),'g-');
        ylim(r_lim(i,:));   ylabel(['x',num2str(i),' (m)']);
    end
    ax(4) = nexttile;   plot(t,v_abs,'b'); ylabel('Velocity (m/s)');
    linkaxes(ax,'x');   xlabel('Time (s)');  xlim('tight');
    sgtitle(['Session ', batdate],'Interpreter','none');
    fig_count = saveFig(analysis_directory,batdate,fig_count,options.savefigures);
    
    %=== Recalculate r and v
    r = smoothdata(r_corr,1,'lowess',round(Fs*0.3));
    v = diff(r,1,1)./diff(t)'; v=cat(1,zeros(1,3),v);   v_abs = vecnorm(v,2,2);
end

%% Visualize raw and processed data

figure('units','normalized','outerposition',[.2 .3 .5 .5]);
tiledlayout(1,3,'TileSpacing','tight');
sgtitle('C3D Raw Data');
nexttile;   plot(c3d_data.Markers(:,:,1),c3d_data.Markers(:,:,2),'.r');   hold on;  plot(r_ft(:,1),r_ft(:,2),'k-');   title('TOP');     axis equal; xticks([]); yticks([]); xlim(r_lim(1,:)*1.1); ylim(r_lim(2,:)*1.1);
nexttile;   plot(c3d_data.Markers(:,:,1),c3d_data.Markers(:,:,3),'.r');   hold on;  plot(r_ft(:,1),r_ft(:,3),'k-');   axis equal; xticks([]); yticks([]);   xlim(r_lim(1,:)*1.1); ylim(r_lim(3,:)*1.1); title('SIDE1');
nexttile;   plot(c3d_data.Markers(:,:,2),c3d_data.Markers(:,:,3),'.r');   hold on;  plot(r_ft(:,2),r_ft(:,3),'k-');   axis equal; xticks([]); yticks([]);   xlim(r_lim(2,:)*1.1); ylim(r_lim(3,:)*1.1); title('SIDE2');
sgtitle(['Session ', batdate],'Interpreter','none');
fig_count = saveFig(analysis_directory,batdate,fig_count,options.savefigures);

%% Flight segmentation and clustering

%=== Flight segmentation
f_smp = []; [bflying,f_num,f_smp(:,1),f_smp(:,2)] = FlightSegm_AF_v0(v_abs,v_th,Fs);

%=== FIGURE: Velocity and flight segmentation
figure('units','normalized','outerposition',[0 0 1 .5]);
area(t,bflying*5,'FaceAlpha',0.3,'LineStyle','none');  hold on;
plot(t,v_abs,'r');     plot(t,r(:,1),'k--');  ylabel('Velocity (m/s)');     hold off;
legend('Fly','Vel','x(m)'); ylim([-3 6]);
title([num2str(f_num) ' flights']);
xlabel('Time(s)');  xlim('tight');
sgtitle(['Session ', batdate],'Interpreter','none');
fig_count = saveFig(analysis_directory,batdate,fig_count,options.savefigures);

%=== Cluster flights
alpha_clus = 1.0;         %Parameter for flight clustering
f_clus = FlightClus_AF_v3(r,bflying,Fs,'Alpha',alpha_clus,'Frechet',1,'Points',7);
sgtitle(['Session ', batdate],'Interpreter','none');
fig_count = saveFig(analysis_directory,batdate,fig_count,options.savefigures);

%=== Plot only flight data
figure('units','normalized','outerposition',[.2 .3 .5 .5]);
tiledlayout(1,3,'TileSpacing','tight');
nexttile;   plot(r(logical(bflying),1),r(logical(bflying),2),'k.');   title('TOP');     axis equal; xticks([]); yticks([]); xlim(r_lim(1,:)*1.1); ylim(r_lim(2,:)*1.1);
nexttile;   plot(r(logical(bflying),1),r(logical(bflying),3),'k.');   axis equal; xticks([]); yticks([]);   xlim(r_lim(1,:)*1.1); ylim(r_lim(3,:)*1.1); title('SIDE1');
nexttile;   plot(r(logical(bflying),2),r(logical(bflying),3),'k.');   axis equal; xticks([]); yticks([]);   xlim(r_lim(2,:)*1.1); ylim(r_lim(3,:)*1.1); title('SIDE2');
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
