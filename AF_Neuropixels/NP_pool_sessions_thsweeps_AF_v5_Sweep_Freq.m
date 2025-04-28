%% POOL DATA RELATED TO LFP AND THETA SWEEPS
%=== Load data
for hide=1
    clr;
    
    %=== Default settings
    set(groot,'defaultAxesFontSize',12);
    set(groot,'defaultHistogramEdgeColor','none','defaultHistogramFaceAlpha',0.5);
    
    %=== Sessions to load (comment sessions you don't want to include)
    sessions2include = {...
        %'Dataset_1','32622','231006';...
        'Dataset_1','32622','231007';...
        'Dataset_1','32622','231008';...
        'Dataset_1','32622','231009';...
        'Dataset_1','32622','231010';...
        'Dataset_2','14445','231208';...
        'Dataset_2','14445','231209';...
        %'Dataset_2','14445','231210';...
        %'Dataset_2','14445','231211';...
        %'Dataset_2','14445','231212';...
        'Dataset_2','14611','231208';...
        'Dataset_2','14611','231209';...
        'Dataset_2','14611','231210';...
        'Dataset_2','14611','231211';...
        %'Dataset_2','14611','231212';...
        'Dataset_3','00000','240402';...
        'Dataset_3','00000','240403';...
        'Dataset_3','00000','240404';...
        'Dataset_3','00000','240405';...
        'Dataset_3','14543','240419';...
        'Dataset_3','14543','240420';...
        'Dataset_3','14543','240421';...
        'Dataset_3','14543','240422';...
        'Dataset_3','29959','240402';...
        'Dataset_3','29959','240403';...
        'Dataset_3','29959','240404';...
        'Dataset_3','29959','240405';...
        ...
        };
    
    %=== Load data and aggregate them in the Multi_Day structure
    Folder = cd;    FileList = dir(fullfile(Folder, '**', 'Analyzed_NPs_*'));   LFP = [];   NP_units = [];  SWPs = []; ids = [];
    for nc = 1:length(FileList)
        cd(FileList(nc).folder);
        load(FileList(nc).name);
        if ~isempty(LFP_PSD) && any(cellfun(@(row) isequal(row, LFP_PSD.unique_ID(1,:)), num2cell(sessions2include, 2)))
            LFP = [LFP; LFP_PSD];
            NP_units = [NP_units;   NP_unit];
            SWPs = [SWPs;   SWP_table];
        end
        disp([num2str(length(FileList)-nc),' remaining sessions to load...']);
    end
    cd(Folder); clearvars -except LFP NP_units SWPs
    
    NP_units = struct2table(NP_units);
    
end

%% ECHOLOCATION PHASE
for hide=1
    
    %=== Look at echolocation phase
    N_flights = size(SWPs,1);
    click_phase = [];
    for zz=1:N_flights
        cond = [SWPs.pos_real{zz,1}/SWPs.pos_real{zz,1}(end)>0.2 & SWPs.pos_real{zz,1}/SWPs.pos_real{zz,1}(end)<.8 & SWPs.clk{zz,1}];
        click_phase  = [click_phase; SWPs.wbt_phase{zz,1}(cond)];
    end
    
    %=== Plot phase
    figure('units','normalized','outerposition',[.3 .3 .2 .25]);
    tiledlayout(1,2,'TileSpacing','compact');
    nexttile;   polarhistogram(click_phase,linspace(-pi,pi,70),'Normalization','probability','facealpha',.9,'edgecolor','none','FaceColor','k');
    nexttile;   histogram([click_phase;click_phase+2*pi],unique([linspace(-pi,pi,40),linspace(pi,3*pi,40)]),'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','k');
    yticks([]); xticks(-pi:pi:3*pi);    xticklabels({'-180', '0', '180', '360', '540'});   xlabel('Wingbeat');
    
end

%% THETA SWEEPS ANALAYSIS (TEMPLATE MATCHING)
for hide=1
    
    %=== Select subtable of good flights
    rng(1);
    SWPs_sst = SWPs(SWPs.rmsDec_error<1.3 & SWPs.prc_decoded>0.7,:);
    [~,~,SWPs_sst.groupID] =  unique([string(SWPs_sst.unique_ID),string(SWPs_sst.flight_id)],'rows');                                  % Id for each cluster from each session
    [~,~,SWPs_sst.sessionID] =  unique([string(SWPs_sst.unique_ID)],'rows');                                                           % Id for each cluster from each session
    N_flights = size(SWPs_sst,1);
    smooth_f = [1 .3];                                                  % [space bin,time bin]
    n_bins_phase = 11;
    bin_size_1D = 0.15;
    Fs_Hi = 1/mean(diff(SWPs_sst.bin_time{1,1}));                       % Sampling frequency (200 Hz)
    
    %=== Create template for finding candidate theta sweeps via template matching
    template_samples = 12;      %12
    sweep_binspan = 6;          %6
    conv_threshold = 2;         %2
    template = normalize(interp1(normpdf(-3:0.1:3),linspace(1,61,template_samples)),'range')*bin_size_1D*sweep_binspan;
    %plot(template)
    
    %=== Find theta sweeps by template matching
    SWPs_posterior = [];  SWPs_wbtphase = [];   SWPs_tmrphase = [];  SWPs_amp = [];   SWPs_wbtphase_ctrl = [];    SWPs_f_id = []; SWPs_s_id = [];   SWPs_flightID = []; SWPs_u_id = {}; SWPs_tmr = [];
    SWPs_acc = {}; SWPs_sample = []; wbt_phase_all = []; tmr_phase_all = []; tmr_power_all = []; fract_pos_all = [];    SWPs_phs = [];  PSD_sweeps = [];    WFR_sweeps = [];    SWPs_clk = [];
    warning('off');
    for zz=1:N_flights
        
        %=== Number of spatial bins
        N_spatial_bins = size(SWPs_sst.p_dec_flight{zz,1},1);
        spt_bin_ids = [1:N_spatial_bins]';
        
        %=== Shifted posterior and wingbeat phase
        sft_posterior = imgaussfilt(SWPs_sst.p_dec_shifted{zz,1},smooth_f);
        wbt_phase = SWPs_sst.wbt_phase{zz,1};
        tmr_phase = SWPs_sst.tmr_phase{zz,1};
        fract_pos = SWPs_sst.pos_real{zz,1}/SWPs_sst.pos_real{zz,1}(end);
        
        wbt_phase_all = [wbt_phase_all; wbt_phase];
        tmr_phase_all = [tmr_phase_all; tmr_phase];
        fract_pos_all = [fract_pos_all; fract_pos];
        tmr_power_all = [tmr_power_all; SWPs_sst.tmr_power{zz,1};];
        
        %=== Extract the center of mass and convolve with template
        %cnt_mass = (spt_bin_ids'*imgaussfilt(SWPs_sst.p_dec_shifted{zz,1},smooth_f)-N_spatial_bins/2)*bin_size_1D;
        [~,max_loc] = max(imgaussfilt(SWPs_sst.p_dec_shifted{zz,1},smooth_f),[],1);
        cnt_mass = (max_loc-N_spatial_bins/2)*bin_size_1D;
        %sample_trace = cnt_mass;
        sample_trace = detrend(cnt_mass,15,'omitnan');
        convolutionResult = conv(sample_trace, flip(template), 'same');
        
        %         figure('units','normalized','outerposition',[.3 .3 .2 .3]);
        %         plot(sample_trace); hold on; plot(normalize(convolutionResult));
        
        %=== Find candidate sweeps
        [candidate_swp_amp,candidate_swp] = findpeaks(convolutionResult,'MinPeakDistance',template_samples,'MinPeakHeight',conv_threshold);
        
        %=== Extract the power spectrum of the center of mass trace
        n = 2^nextpow2(numel(cnt_mass));                    % Optimal number of points for FFT
        Y = fft(cnt_mass,n);                                % Calculate fft
        P = abs(Y/n).^2;                                    % Power density at all frequences (positive and negative)
        f = Fs_Hi*(0:n/2)/n;                                % Fs is the sampling frequency
        PSD = P(1:n/2+1);                                   % Power spectral density
        sm_PSD = smoothdata(PSD,'movmedian',n/Fs_Hi*2);     % Smoothed
        %[~,max_loc] = max(sm_PSD); f_wBeats = f(max_loc);  % Characteristic frequency of wingbeat signal
        PSD_sweeps = [PSD_sweeps; interp1(f,sm_PSD,linspace(0,100,200))];
        WFR_sweeps = [WFR_sweeps; mean(instfreq(SWPs_sst.wbt{zz,1},Fs_Hi,'Method','hilbert'))];
        
        %=== For every sweep, store the posterior probability and the wingbeat phase
        for ss=1:numel(candidate_swp)
            if candidate_swp(ss)-template_samples>0 && candidate_swp(ss)+template_samples< size(sft_posterior,2) && round(N_spatial_bins/2)-1*sweep_binspan>0 && round(N_spatial_bins/2)+3*sweep_binspan<N_spatial_bins
                tmp_pst = sft_posterior(round(N_spatial_bins/2)+[-1*sweep_binspan:2*sweep_binspan],candidate_swp(ss)+[-template_samples:template_samples]);
                tmp_pst(isnan(tmp_pst))=0;
                SWPs_posterior = cat(3,SWPs_posterior,tmp_pst);
                SWPs_wbtphase = [SWPs_wbtphase; wbt_phase(candidate_swp(ss))];
                SWPs_tmrphase = [SWPs_tmrphase; tmr_phase(candidate_swp(ss))];
                
                %=== Control random phase
                available_samples = find(fract_pos>0.0 & fract_pos<1);
                SWPs_wbtphase_ctrl = [SWPs_wbtphase_ctrl; wbt_phase(available_samples(randi(numel(available_samples),1,100)))];
                SWPs_amp = [SWPs_amp; candidate_swp_amp(ss)];
                
                SWPs_f_id = [SWPs_f_id; SWPs_sst.groupID];
                SWPs_s_id = [SWPs_s_id; SWPs_sst.sessionID];
                SWPs_u_id = [SWPs_u_id; SWPs_sst.unique_ID(zz,:)];
                SWPs_flightID = [SWPs_flightID; zz];
                SWPs_sample = [SWPs_sample; candidate_swp(ss)];
                SWPs_acc = [SWPs_acc; {SWPs_sst.wbt{zz,1}}];
                SWPs_phs = [SWPs_phs; {SWPs_sst.wbt_phase{zz,1}}];
                SWPs_clk = [SWPs_clk; {SWPs_sst.clk{zz,1}}];
                SWPs_tmr = [SWPs_tmr; {SWPs_sst.LFP{zz,1}}];
                
            end
        end
        %SWPs_sst.pos_real{zz,1}/SWPs_sst.pos_real{zz,1}(end);
        
    end
    warning('on');
    
    %=== Train the network or load the pretrained network if needed
    if ~isfile('trainedNet.mat')
        theta_sweep_manual_class_AF_v0
    else
        load('trainedNet.mat');
        load('man_class_SWSs.mat');
    end
    
    %=== Prepare the images for prediction
    numImages = size(SWPs_posterior, 3); % Number of images
    inputSize = [227, 227, 3];           % AlexNet input size
    test_data_sst = SWPs_posterior;      % Copy data
    
    %=== Resize images for AlexNet
    SWPs_resized = zeros([inputSize, numImages], 'single');
    for i = 1:numImages
        test_data_sst(:, :, i) = SWPs_posterior(:, :, i)./max(SWPs_posterior(:, :, i),[],'all');
        SWPs_resized(:, :, :, i) = imresize(repmat(test_data_sst(:, :, i), [1, 1, 3]), inputSize(1:2));
    end
    
    %=== Predict the labels using the trained network
    testDs = augmentedImageDatastore(inputSize, SWPs_resized);
    predictedLabels = classify(trainedNet, testDs);
    good_ids = find(predictedLabels == '1');  N_good = numel(good_ids);
    baad_ids = find(predictedLabels == '0');  N_baad = numel(baad_ids);
    
    %=== Show some examples
    rng(3)
    figure('units','normalized','outerposition',[0 0 .5 1]);
    tiledlayout(10,10,'TileSpacing','tight');
    for i=1:min(N_good,100)
        tmp_pst = squeeze(SWPs_posterior(:,:,good_ids(randi(N_good))));
        nexttile;   imagesc(tmp_pst,prctile(tmp_pst, [1 99],'all')');
        colormap(flipud(gray)); axis off;  set(gca,'YDir','normal');
    end
    sgtitle('Good Sweeps');
    figure('units','normalized','outerposition',[.5 0 .5 1]);
    tiledlayout(10,10,'TileSpacing','compact');
    for i=1:min(N_baad,100)
        tmp_pst = squeeze(SWPs_posterior(:,:,baad_ids(randi(N_baad))));
        nexttile;   imagesc(tmp_pst,prctile(tmp_pst, [1 99],'all')');
        colormap(flipud(gray)); axis off;  set(gca,'YDir','normal');
    end
    sgtitle('Bad Sweeps');
    
    %% ==================================================== ANALYSIS OF SWEEP DT VS WINGBEAT FREQUENCY =========================================================================
    
    S2W = table();                                          % Sweep to wingbeat table
    S2W.flightID = SWPs_flightID(good_ids);                 % Id of the flights
    S2W.acc = SWPs_acc(good_ids);                           % Accelerometer
    S2W.phs = SWPs_phs(good_ids);                           % Wingbeat Phase
    S2W.sample = SWPs_sample(good_ids);                     % Sample of the sweep
    S2W.u_id = SWPs_u_id(good_ids,:);                       % Unique identifier
    [~,~,S2W.sessionID] = unique([string(S2W.u_id)],'rows');% Identifier
    S2W.unwrp_phs = cellfun(@(x) unwrap(x(~isnan(x))),S2W.phs,'UniformOutput',false);
    S2W.wbt_cycles = cellfun(@(x) (x(end)-x(1))./(2*pi),S2W.unwrp_phs);
    freq_Hbins = linspace(0,20,20);
    
    %=== Count number of sweeps per flight
    [sweep_count,sweep_flight] = groupcounts(S2W.flightID);
    good_flights = sweep_flight(sweep_count>1);
    
    S2W_summary = groupsummary(S2W,'flightID','mean',{'wbt_cycles'});
    S2W_summary.ratio = S2W_summary.GroupCount./S2W_summary.mean_wbt_cycles;
    
    mean_and_sem_AF_v0(S2W_summary.GroupCount)
    mean_and_sem_AF_v0(S2W_summary.ratio)
    
    %=== Accumulate a few variables
    dt_S2W = [];    f_est1_S2W = [];    f_est2_S2W = []; n_cyc = [];    id_tmp = [];
    for i=good_flights'
        S2W_sst = S2W(S2W.flightID == i,:);
        for j = 1:size(S2W_sst,1)-1
            wbt_signal = S2W_sst.acc{j};
            wbt_phase = unwrap(S2W_sst.phs{j});
            wbt_signal_freq = instfreq(wbt_signal,Fs_Hi,'Method','hilbert');
            dt_tmp = (S2W_sst.sample(j+1)-S2W_sst.sample(j))/Fs_Hi;
            
            dt_S2W = [dt_S2W; dt_tmp];
            f_est1_S2W = [f_est1_S2W; mean(wbt_signal_freq(S2W_sst.sample(j):S2W_sst.sample(j+1)))];
            f_est2_S2W = [f_est2_S2W; (wbt_phase(S2W_sst.sample(j+1))-wbt_phase(S2W_sst.sample(j)))/(2*pi*dt_tmp)];
            n_cyc = [n_cyc; (wbt_phase(S2W_sst.sample(j+1))-wbt_phase(S2W_sst.sample(j)))/(2*pi)];
            id_tmp = [id_tmp; S2W_sst.u_id(j,:)];
            
        end
    end
    
    %=== Define max time between consecutive sweeps
    consecutive_dt = 1.5/mean(f_est1_S2W);
    
    %=== Calculate frequency of consecutive sweeps and corresponding wingbeat frequency
    swp_freq = 1./dt_S2W(dt_S2W<consecutive_dt);
    wbt_freq = f_est1_S2W(dt_S2W<consecutive_dt);
    id_tmp_cons = id_tmp(dt_S2W<consecutive_dt,:);
    unique([string(id_tmp_cons)],'rows');
    
    figure('units','normalized','outerposition',[.2 .3 .45 .3]);
    tiledlayout(1,4,'TileSpacing','compact');
    nexttile;   histogram(wbt_freq,freq_Hbins,'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','k'); xlim([4 16]);
    nexttile;   histogram(swp_freq,freq_Hbins,'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','k'); xlim([4 16]);
    nexttile;   scatter(wbt_freq,swp_freq,20,'MarkerFaceColor','k','MarkerEdgeColor','none','MarkerFaceAlpha',.1);
    xlabel('Wingbeat Frequency (Hz)');  ylabel('Estimated Sweep Frequency (Hz)');  axis('square');
    [corr_sw,p_val_sw] = corr(wbt_freq,swp_freq,'type','Spearman');
    title([corr_sw,p_val_sw]);  xlim([4 16]);   ylim([4 16]);   %xlim(prctile(wbt_freq,[1 99]));
    nexttile;   scatter(wbt_freq,swp_freq,20,'MarkerFaceColor','k','MarkerEdgeColor','none','MarkerFaceAlpha',.1);  xlim(prctile(wbt_freq,[2 98])); axis('square');
    sgtitle('Consecutive Sweeps Analysis');
    
    %=== ANALYSIS ON ALL FLIGHTS
    %=== Plot power spectral density of the sweeps and fit
    unique([string(SWPs_u_id)],'rows')
    freq_PSD_sweeps = linspace(0,100,200)';
    expEqn = 'a*exp(-b*x+c)+d';
    x = freq_PSD_sweeps;
    %y1 = mean(PSD_sweeps)';
    y1 = mean(PSD_sweeps./sum(PSD_sweeps,2))';
    f1 = fit(x,y1,expEqn,'Start',[max(y1) 1 0 0],'Lower',[0 -inf 0 0],'Exclude',x<2);
    figure('units','normalized','outerposition',[.1 .1 .35 .3]);
    tiledlayout(1,3,'TileSpacing','tight');
    nexttile;   plot(freq_PSD_sweeps,mean(PSD_sweeps./sum(PSD_sweeps,2)));  set(gca, 'YScale', 'log');   xlim([2 20]);  xlabel('Frequency (Hz)');   ylabel('PSD (norm)');
    nexttile;   plot(x,feval(f1,x),'k--');   hold on; plot(x,y1);           set(gca, 'YScale', 'log');   xlim([2 20]);  xlabel('Frequency (Hz)');   ylabel('PSD (norm)');
    nexttile;   plot(x,y1-feval(f1,x)); hold on; xlim([2 20]);  xlabel('Frequency (Hz)');   ylabel('PSD (norm)');   title('Exp. Decay Subtracted');
    
    %=== Find peak frequency of internal representations for each flight
    peak_swp_psd = [];
    for i=1:size(PSD_sweeps,1)
        y1 = PSD_sweeps(i,:)'./sum(PSD_sweeps(i,:));
        f1 = fit(x,y1,expEqn,'Start',[max(y1) 1 0 0],'Lower',[0 -inf 0 0],'Exclude',x<2);
        y2 = y1-feval(f1,x);
        y2 = y2(x>5 & x<16);
        x2 = x(x>5 & x<16);
        [~,max_tmp] = max(y2);
        peak_swp_psd = [peak_swp_psd; x2(max_tmp)];
    end
    
    figure('units','normalized','outerposition',[.1 .1 .4 .3]);
    tiledlayout(1,4,'TileSpacing','tight');
    nexttile;   histogram(WFR_sweeps,freq_Hbins,'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','k');   xlim([4 16]);
    nexttile;   histogram(peak_swp_psd,freq_Hbins,'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','k'); xlim([4 16]);
    nexttile;   scatter(WFR_sweeps,peak_swp_psd,15,'MarkerFaceColor','k','MarkerEdgeColor','none','MarkerFaceAlpha',.1);
    xlabel('Wingbeat Frequency (Hz)');  ylabel('Sweep Peak Frequency (Hz)');
    [corr_sw,p_val_sw] = corr(WFR_sweeps,peak_swp_psd,'type','Spearman');   axis('square');
    title([corr_sw,p_val_sw]);  xlim([4 16]);   ylim([4 16])   %xlim(prctile(WFR_sweeps,[0 100]));
    nexttile;   scatter(WFR_sweeps,peak_swp_psd,20,'MarkerFaceColor','k','MarkerEdgeColor','none','MarkerFaceAlpha',.1);  xlim(prctile(WFR_sweeps,[2 98])); axis('square');
    sgtitle('Whole Flight Analysis');
    
    mean_and_sem_AF_v0(swp_freq)
    mean_and_sem_AF_v0(peak_swp_psd)
    
    %% === Simulate level of correlation and evaluate results
    
    %=== Evaluate changes in sample size
    r_sim = 0.1;
    sampleSize = 3:1000;
    meanSimCorr = zeros(numel(sampleSize),1);   stdevSimCorr = zeros(numel(sampleSize),1);
    meanSimPval = zeros(numel(sampleSize),1);   stdevSimPval = zeros(numel(sampleSize),1);
    for nn=sampleSize
        n_sst = nn;
        sim_corr = [];sim_pval = [];
        for i=1:50
            wbt_freq_sst = datasample(wbt_freq,n_sst);
            z_tmp = randn(size(wbt_freq_sst));
            swp_freq_sim = (r_sim*zscore(wbt_freq_sst)+sqrt(1-r_sim^2)*z_tmp)+mean(swp_freq);
            [corr_tmp,p_val_tmp] = corr(wbt_freq_sst,swp_freq_sim);
            sim_corr = [sim_corr; corr_tmp];
            sim_pval = [sim_pval; p_val_tmp];
        end
        
        meanSimCorr(nn-2) = mean(sim_corr);
        meanSimPval(nn-2) = mean(sim_pval);
        stdevSimCorr(nn-2) = std(sim_corr);
        stdevSimPval(nn-2) = std(sim_pval);
        
    end
    
    figure('units','normalized','outerposition',[.1 .1 .3 .3]);
    tiledlayout(1,3,'TileSpacing','tight');
    nexttile;   scatter(wbt_freq_sst(1:numel(wbt_freq)),swp_freq_sim(1:numel(wbt_freq)),20,'MarkerFaceColor','k','MarkerEdgeColor','none','MarkerFaceAlpha',.5);  xlim(prctile(wbt_freq_sst,[2 98]));
    nexttile;   plotWinterval_AF_v0(sampleSize,meanSimCorr,meanSimCorr-stdevSimCorr,meanSimCorr+stdevSimCorr,'k');   hold on; plot(xlim,r_sim.*[1 1],'k--'); plot(numel(wbt_freq)*[1 1],ylim);
    nexttile;   plotWinterval_AF_v0(sampleSize,meanSimPval,meanSimPval-stdevSimPval,meanSimPval+stdevSimPval,'k');   hold on; plot(xlim,[.05 .05],'k--');    plot(numel(wbt_freq)*[1 1],ylim);
    
    %=== Evaluate changes in r
    r_sim_values = linspace(0, 0.5, 100);
    n_trials = 50;
    meanSimCorr = zeros(numel(r_sim_values),1);
    stdevSimCorr = zeros(numel(r_sim_values),1);
    meanSimPval = zeros(numel(r_sim_values),1);
    stdevSimPval = zeros(numel(r_sim_values),1);
    meanPowerEst = zeros(numel(r_sim_values),1);
    stdevPowerEst = zeros(numel(r_sim_values),1);
    
    % Use z-scored version of x (fixed variance)
    x = zscore(wbt_freq);
    n_samples = numel(x);
    swp_mean = mean(swp_freq);
    
    for nn = 1:numel(r_sim_values)
        r_sim = r_sim_values(nn);
        sim_corr = zeros(n_trials,1);
        sim_pval = zeros(n_trials,1);
        for i = 1:n_trials
            z = randn(n_samples, 1);
            y = r_sim * x + sqrt(1 - r_sim^2) * z + swp_mean;
            [r, p] = corr(x, y,'type','Spearman');
            sim_corr(i) = r;
            sim_pval(i) = p;
        end
        
        % Store stats
        meanSimCorr(nn) = mean(sim_corr);
        stdevSimCorr(nn) = std(sim_corr);
        meanSimPval(nn) = mean(sim_pval);
        stdevSimPval(nn) = std(sim_pval);
        meanPowerEst(nn) = mean(sim_pval < 0.05);  % Proportion of trials detecting significance
        stdevPowerEst(nn) = std(sim_pval < 0.05);  % Proportion of trials detecting significance
    end
    
    
    figure('units','normalized','outerposition',[.1 .1 .4 .3]);
    tiledlayout(1,4,'TileSpacing','tight');
    nexttile;   scatter(wbt_freq,y,20,'MarkerFaceColor','k','MarkerEdgeColor','none','MarkerFaceAlpha',.3); axis('square');     xlim([7.3 9]);   ylim([6 13]);
    nexttile;   plotWinterval_AF_v0(r_sim_values,meanSimCorr,meanSimCorr-stdevSimCorr,meanSimCorr+stdevSimCorr,'k');   hold on;
    nexttile;   plotWinterval_AF_v0(r_sim_values,meanSimPval,meanSimPval-stdevSimPval,meanSimPval+stdevSimPval,'k');   hold on; plot(xlim,[.05 .05],'k--');
    nexttile;   plotWinterval_AF_v0(r_sim_values,meanPowerEst,meanPowerEst-stdevPowerEst,meanPowerEst+stdevPowerEst,'k');   hold on; plot(xlim,[.05 .05],'k--');
    
    
    
    
end
