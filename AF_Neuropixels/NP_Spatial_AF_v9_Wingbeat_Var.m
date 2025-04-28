function NP_Spatial_AF_v9_Wingbeat_Var(folder_name)
%% SCRIPT FOR LOOKING AT NP RECORDINGS WITH MULTIPLE PROBES IN THE HIPPOCAMPUS DURING NAVIGATION
% Updated by A.F. on May 2024, based on previous versions
% BRIEF DESCRIPTION (TBD)

fig_visibility = 'on';

%% LOAD DATA and SAVING OPTIONS
for hide=1
    
    %=== Load data
    disp('Loading data...');
    bhvFile = dir('Extracted_Behavior_*');  load([bhvFile.folder '/' bhvFile.name]);    % Behavioral file
    imuFile = dir('IMU_data.mat');          load([imuFile.folder '/' imuFile.name]);    % IMU data
    unique_ID = options.unique_ID;                                                      % Session identifier
   
    %=== Echolocation Clicks
    if ~isempty(dir(fullfile('Detected_Clicks_*')))
        Clck_file = dir(fullfile('Detected_Clicks_*'));                     % File with times of the detected echolocation clicks
        load([Clck_file.folder,'\',Clck_file.name]);
    else
        Detected_Clicks =[];
        Detected_Clicks.times = [];
    end
    
    %=== Create analysis folder for storing the results
    options.savedata = 1;                                   % If creating folder for saving the data or not
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
    
end

%% MAKE SURE BEHAVIOR AND EPHYS HAVE THE SAME DURATION
for hide=1
    
    % Make sure t is a column vector
    if isrow(t),t=t';end
    
    %=== Show difference between recording times
    dt_str = NP_imu.t(1)-t(1);
    dt_end = NP_imu.t(end)-t(T);
    disp(['NP started ',num2str(-dt_str),' s before behavior']);
    disp(['NP ended '  ,num2str(dt_end),' s after behavior']);
    
    %=== Cut the end of recordings accordingly
    if dt_end<0
        %=== Find end of recording
        t_max = NP_imu.t(end);
        T_max = find(diff(t<t_max)==-1);
        
        %=== Cut behavioral variables
        T = T_max;
        bflying = bflying(1:T_max,1);
        r = r(1:T_max,:);
        v = v(1:T_max,:);
        t = t(1:T_max,1);
        v_abs = v_abs(1:T_max,1);
        f_smp = f_smp(all(f_smp<T_max,2),:);
        f_num = size(f_smp,1);
    end
    
end

%% BASIC PARAMS, ASSIGNMENTS
for hide=1
    
    %=== Recluster Flights (default is .8)
    alpha_clus = .8;
    if strcmp(unique_ID{1},'Dataset_2'),alpha_clus = 1.2;end    % More liberal clustering for Field Station Dataset
    clear f_clus;
    f_clus = FlightClus_AF_v3(r,bflying,Fs,'Alpha',alpha_clus,'Frechet',1,'Points',7);
    while numel(unique(f_clus.id))==1 && alpha_clus<1.5         % Try increasing the clustering param if only 1 cluster
        alpha_clus = alpha_clus+0.1;
        clear f_clus;   close all;
        f_clus = FlightClus_AF_v3(r,bflying,Fs,'Alpha',alpha_clus,'Frechet',1,'Points',7);
    end
    sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3}],'Interpreter','none');
    fig_count = saveFig(analysis_directory,[cell2mat(unique_ID)],fig_count,options.savefigures);
    
    %=== Params
    if isrow(t),t=t';end                                    % Make sure t is a column vector
    bflying = logical(bflying);                             % Convert to logical
    imp_bat_clusters =  unique(f_clus.id);                  % Surviving fligth clusters
    n_surv_clusters = numel(imp_bat_clusters);              % Number of surviving flight clusters
    col_clus = hsv(n_surv_clusters);                        % Colors for clusters
    t_Hi = NP_imu.t;                                        % Time for fast sampling (500 Hz, IMU and LFP)
    bflying_Hi = logical(interp1(t,double(bflying),t_Hi,'nearest',0));
    Fs_Hi = NP_imu.Fs;                                      % Sampling frequency for fast sampling (500 Hz, IMU and LFP)
    
    %=== Extract flight periods from accelerometer
    a_abs_NP = vecnorm(NP_imu.acc,2,2);                     % Absolute acceleration
    a_flt_NP = bandpass(a_abs_NP,[6 10],Fs_Hi);             % Filtered at the wing-beat frequency
    v_abs_Hi = interp1(t,v_abs,t_Hi,'linear',0);            % Velocity at the IMU sampling
    [up,lo] = envelope(a_flt_NP,round(0.06*Fs_Hi),'peak');  % Upper and lower envelopes
    env = normalize(up - lo,'range');                       % Amplitude of the envelope
    env_th = otsuthresh(histcounts(env));                   % Threshold (based on Otsu method). Can be set at 0.35
    wBeats = movsum(env>env_th,2*Fs_Hi)>Fs_Hi/5;            % Euristic criterion for flight detection
    a_mdg = interp1(t_Hi,movmean(abs(a_abs_NP-1),Fs_Hi*0.5),t,'linear',0);  % mean deviation from g in a 0.5s integration window
    immobile = a_mdg<0.04;                                  % Immobility vector
    wb_phase = interp1(t_Hi,angle(hilbert(a_flt_NP)),t,'linear',0);  % Wingbeat phase at t sampling
    %wb_phase = angle(hilbert(a_flt_NP));
    
    %=== Calculate characteristic frequency of wingbeat signal
    n = 2^nextpow2(numel(a_flt_NP(bflying)));           % Optimal number of points for FFT
    Y = fft(a_flt_NP(bflying),n);                       % Calculate fft
    P = abs(Y/n).^2;                                    % Power density at all frequences (positive and negative)
    f = Fs_Hi*(0:n/2)/n;                                % Fs is the sampling frequency
    PSD = P(1:n/2+1);                                   % Power spectral density
    sm_PSD = smoothdata(PSD,'movmedian',n/Fs_Hi*3);     % Smoothed
    [~,max_loc] = max(sm_PSD); f_wBeats = f(max_loc);   % Characteristic frequency of wingbeat signal
    %     plot(f,sm_PSD,'LineWidth',2,'Color','r'); set(gca, 'YScale', 'log');
    
    %=== Cluster Positional Data and find feeders
    r_fd = [2.77,0.82,1.75; 2.79,-0.99,1.64];                                       % Feeder coordinates
    if strcmp(unique_ID{1,1},'Dataset_2')
        r_fd = [-3,0.82,2; 100,100,100];                                            % Correct coordinates for Field Station
    end
    X = r(~bflying & v_abs<0.5,:);                                                  % Consider only stationary epochs
    [~,centroid,~,~,~] = Cluster3D_AF_v1(X,30*Fs,0.3,100);                          % Perform clustering
    for i=1:2                                                                       % Correct feeder position(s)
        [feeder_distance,fd_idx] = min(vecnorm(r_fd(i,:)-centroid,2,2));
        if feeder_distance<0.2,r_fd(i,:) =  centroid(fd_idx,:);end                  % Do not correct if further than 20 cm
    end
    bat_on_feeder = vecnorm(r_fd(1,:)-r,2,2)<0.3 | vecnorm(r_fd(2,:)-r,2,2)<0.3;    % Bat on feeder
    
    %=== Transform start/stop samples of flight to the IMU sampling
    f_smp_Hi = f_smp;
    for i=1:f_num
        if t(f_smp(i,1))>=NP_imu.t(1) && t(f_smp(i,2))<=NP_imu.t(end)
            f_smp_Hi(i,1) = knnsearch_fast_AF_v0(NP_imu.t',t(f_smp(i,1)),0.1);
            f_smp_Hi(i,2) = knnsearch_fast_AF_v0(NP_imu.t',t(f_smp(i,2)),0.1);
        else
            f_smp_Hi(i,1) = NaN;
            f_smp_Hi(i,2) = NaN;
        end
    end
    f_smp_Hi = f_smp_Hi(~isnan(f_smp_Hi(:,1)),:);    f_num_Hi = size(f_smp_Hi,1);
    
    %=== Extract relevant phases and amplitudes
    wbt_phase_tmp = wrapTo2Pi(angle(hilbert(a_flt_NP)))-pi;   wbt_power_tmp = abs(hilbert(a_flt_NP)).^2;    wbt_freq_tmp = instfreq(a_flt_NP,Fs_Hi,'Method','hilbert'); 
    
    %=== Invalidate wingbeat outside of flight
    wbt_phase = NaN(size(wbt_phase_tmp));
    wbt_power = zeros(size(wbt_phase_tmp));
    wbt_freq = NaN(size(wbt_phase_tmp));
    for i=1:f_num_Hi
        wbt_phase(f_smp_Hi(i,1):f_smp_Hi(i,2))= wbt_phase_tmp(f_smp_Hi(i,1):f_smp_Hi(i,2));
        wbt_power(f_smp_Hi(i,1):f_smp_Hi(i,2))= wbt_power_tmp(f_smp_Hi(i,1):f_smp_Hi(i,2));
        wbt_freq(f_smp_Hi(i,1):f_smp_Hi(i,2))=  wbt_freq_tmp(f_smp_Hi(i,1):f_smp_Hi(i,2));
    end
    
    %=== Create structure for storing miscellanea variables
    MSC = struct();
    MSC.unique_ID = unique_ID;
    MSC.wbf_all = f_wBeats;
    
end

%% ADD A FEW FEATURES TO THE FCLUS STRUCTURE
for hide=1
    
    for jj=1:f_clus.N  
    
        %=== Extract flight 3D path
        flight_pos = f_clus.pos(:,:,jj);
        flight_pos = flight_pos(:,~isnan(flight_pos(1,:)))';
    
        %=== Calculate curvature
        norm_grad3D = diff(flight_pos,1,1)./vecnorm(diff(flight_pos,1,1),2,2);
        curv3D_1 = [0; vecnorm(diff(norm_grad3D,1,1),2,2); 0];
        [~,R_tmp,~] = curvature(flight_pos); % from Are Mjaavatten
        curv3D_2 = 1./R_tmp;
        f_clus.curv_1(1,jj) = {curv3D_1};
        f_clus.curv_2(1,jj) = {curv3D_2};
        middle_int = round(numel(curv3D_1)/2)+[-f_clus.Fs/2:f_clus.Fs/2];
        f_clus.med_curv1(1,jj) = median(curv3D_1(middle_int),'omitnan');
        f_clus.med_curv2(1,jj) = median(curv3D_2(middle_int),'omitnan');

        %=== Fit with Gaussian velocity profile
        if f_clus.id(jj)>1
            mean_v1D = diff(f_clus.lin_tr{1, jj}).*f_clus.Fs*f_clus.length(jj);
            mean_t1D = linspace(0,numel(mean_v1D)/f_clus.Fs,numel(mean_v1D));
            [fg1,gof1] = fit(mean_t1D',mean_v1D,'gauss1','StartPoint',[6 1 1]);
            [fg2,gof2] = fit(mean_t1D',mean_v1D,'gauss2','StartPoint',[6 1 1 6 3 1]);
            f_clus.Rsq1(1,jj) = gof1.adjrsquare;
            f_clus.Rsq2(1,jj) = gof2.adjrsquare;
        else
            f_clus.Rsq1(1,jj) = NaN;
            f_clus.Rsq2(1,jj) = NaN;
        end
         
        %=== Add accelerometer profile and frequencies
        if t(f_clus.strt_frame(jj))>t_Hi(1) && t(f_clus.stop_frame(jj))<t_Hi(end)
            strt_frame_hi = knnsearch_fast_AF_v0(t_Hi',t(f_clus.strt_frame(jj)),0.1);
            stop_frame_hi = knnsearch_fast_AF_v0(t_Hi',t(f_clus.stop_frame(jj)),0.1);
            f_clus.a_abs_NP(1,jj) = {a_abs_NP(strt_frame_hi:stop_frame_hi)};
            f_clus.a_flt_NP(1,jj) = {a_flt_NP(strt_frame_hi:stop_frame_hi)};
            f_clus.wbt_freq(1,jj) = {wbt_freq_tmp(strt_frame_hi:stop_frame_hi)};
            f_clus.v_abs_Hi(1,jj) = {v_abs_Hi(strt_frame_hi:stop_frame_hi)};
            
            x = a_abs_NP(strt_frame_hi:stop_frame_hi);
            n = length(x);                      % Length of the signal
            X = fft(x);                         % Perform FFT
            f = (0:n-1)*(Fs_Hi/n);              % Frequency vector
            magnitude = abs(X);                 % Magnitude of FFT
            half_n = floor(n/2) + 1;            % Keep only half
            f = f(1:half_n);                    % Half frequency
            magnitude = magnitude(1:half_n);    % Half spectrum
            M = magnitude(f>6 & f<10);          % Zoom in
            F = f(f>6 & f<10);                  % Zoom in
            [~, idx_peak] = max(M);             % Find the peak frequency (index of the maximum magnitude)
            f_clus.wb_frequency_FFT(1,jj) = F(idx_peak);                                                 % Peak frequency from FFT
            f_clus.wb_frequency_HLB(1,jj) = median(wbt_freq_tmp(strt_frame_hi:stop_frame_hi),'omitnan'); % Peak frequeny from Hilbert transform
            
        else
            f_clus.a_abs_NP(1,jj) = {NaN};
            f_clus.a_flt_NP(1,jj) = {NaN};
            f_clus.wbt_freq(1,jj) = {NaN};
            f_clus.v_abs_Hi(1,jj) = {NaN};
            f_clus.wb_frequency_FFT(1,jj) = NaN;
            f_clus.wb_frequency_HLB(1,jj) = NaN;
        end
    end
    
    MSC.f_clus = f_clus;
end

%% SAVE THE DATA
for hide=1
    if options.savedata
        save([analysis_directory,'/Analyzed_NPs_', unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3},'.mat'],'MSC');
    end
end

