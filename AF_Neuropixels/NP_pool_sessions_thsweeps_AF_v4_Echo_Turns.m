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

%% CURVATURE CALCULATION AND ANALYSIS AROUND LOOPS
for hide=1
    
    Fs_Hi = 1/mean(diff(SWPs.bin_time{1,1}));
    %=== Calculate curvature for each flight, starting from the 3D coordinates
    %figure('units','normalized','outerposition',[.3 .3 .4 .3]);
    for i = 1:size(SWPs, 1)
        
        r = SWPs.r_real{i, 1};
        r = smoothdata(r,1,'gaussian',10);
        
        %=== Calculate curvature
        norm_grad3D = diff(r,1,1)./vecnorm(diff(r,1,1),2,2);        % Intuitive definition
        curv3D_1 = [0; vecnorm(diff(norm_grad3D,1,1),2,2); 0];
        [~, R_tmp, ~] = curvature(r);                               % from Are Mjaavatten
        curv3D_2 = 1 ./ R_tmp;
        
        %=== Flag if coming from a session without echolocation
        SWPs.mic_flag(i) = 1;
        if strcmp(SWPs.unique_ID{i,1},'Dataset_2')
            SWPs.mic_flag(i) = 0;
        end
        
        %=== Forse start/stop to low curvature and smooth
        curv3D_1(1:round(0.25*numel(curv3D_1)))= .005;   curv3D_1(round(0.75*numel(curv3D_1)):end)= .005;
        curv3D_2(1:round(0.25*numel(curv3D_2)))= .2;   curv3D_2(round(0.75*numel(curv3D_2)):end)= .2;
        SWPs.curv3D_1{i,1} = smoothdata(curv3D_1,'gaussian',10);
        SWPs.curv3D_2{i,1} = smoothdata(curv3D_2,'gaussian',10);
        
        %=== Assign turn flag
        SWPs.curv_flag{i,1} = SWPs.curv3D_2{i,1}.*0;
        [max_curv,max_loc]= max(SWPs.curv3D_2{i,1});
        if max_curv>1
            SWPs.curv_flag{i,1}(max_loc) = max_curv;
        end
        
        
        %         %=== Uncomment to show for debugging
        %         tiledlayout(1, 3);        % Set up 1x2 layout
        %         nexttile(1);     plot(r(:, 1), r(:, 2));
        %         nexttile(2);    plot(SWPs.curv3D_1{i,1});
        %         nexttile(3);    plot(SWPs.curv3D_2{i,1});
        %         choice = questdlg('Continue', 'Class Verification','Yes', 'No', 'Stop','Yes');
        %         switch choice
        %             case 'Stop'
        %                 break;
        %         end
        
    end
    
    %=== Add a couple more features
    SWPs.has_turn = cellfun(@max,SWPs.curv_flag)>0;                                                                                 % Flag the fligtht if a turn is detected
    SWPs.wbt_freq = cellfun(@(x) instfreq(smoothdata(x,'gaussian',30),Fs_Hi,'Method','hilbert'),SWPs.wbt,'UniformOutput',false);    % Calculate instantaneous wingbeat frequency
    SWPs.speed = cellfun(@(x) [0; diff(smoothdata(x,'gaussian',30))].*Fs_Hi,SWPs.pos_real,'UniformOutput',false);                   % Calculate speed
    SWPs.trn = cellfun(@(x) double(movmax(x,[0 round(Fs_Hi*0.5)])>0),SWPs.curv_flag,'UniformOutput',false);                         % Extend the turn flag for 500 ms before the turn (value 1)
    unique([string(SWPs.unique_ID)],'rows')
    
    %=== Show distribution of max curvature
    figure('units','normalized','outerposition',[.3 .3 .15 .3]);
    histogram(cellfun(@max,SWPs.curv3D_2),10.^([-1:0.1:1]),'Normalization','probability');   set(gca, 'XScale', 'log');
    set(gca, 'YScale', 'log');  xlabel('Max Curvature');    hold on; plot([1 1],ylim,'r--');
    title(sum(SWPs.has_turn)./size(SWPs, 1));
    
    %=== Show some example flights
    figure('units','normalized','outerposition',[.1 .3 .3 .3]);
    tiledlayout(1,3,'TileSpacing','compact');
    for i=[82,94,118]
        
        x_example = SWPs.r_real{i, 1}(:,1)';
        y_example = SWPs.r_real{i, 1}(:,2)';
        c_example = normalize(smoothdata(SWPs.curv3D_2{i,1},'gaussian',50),'range')';         % color values
        
        
        % Prepare colormap
        cmap = colormap(jet(256));
        nColors = size(cmap, 1);
        cIdx = round(rescale(c_example, 1, nColors));
        
        nexttile;   hold on;
        for j = 1:(length(x_example)-1)
            col = cmap(cIdx(j), :);
            patch([x_example(j), x_example(j+1)],[y_example(j), y_example(j+1)],[0, 0], col,'EdgeColor', col, 'LineWidth', 4, 'FaceColor', 'none');
        end
        axis equal; title(i);    xlim([-2.5 2.8]);  ylim([-3 3]);   colorbar;
        
    end
    
    
    %%
    %=== Look at wingbeat and echolocation around turns
    SWPs_loops = SWPs(SWPs.has_turn & SWPs.mic_flag,:);
    interval = [-round(Fs_Hi*1):round(Fs_Hi*1)];
    clk_array = [];  wfr_array = [];    spd_array = [];    crv_array = [];  lfp_array = [];  spg_array = [];  wbt_array = [];
    for i = 1:size(SWPs_loops,1)
        
        [~,sample_tmp] = max(SWPs_loops.curv_flag{i,1});
        
        if sample_tmp>interval(end) && sample_tmp+interval(end)<numel(SWPs_loops.clk{i,1})
            clk_array = [clk_array;   SWPs_loops.clk{i,1}(sample_tmp+interval)'];
            wfr_array = [wfr_array;   SWPs_loops.wbt_freq{i,1}(sample_tmp+interval)'];
            spd_array = [spd_array;   SWPs_loops.speed{i,1}(sample_tmp+interval)'];
            crv_array = [crv_array;   SWPs_loops.curv3D_2{i,1}(sample_tmp+interval)'];
            lfp_array = [lfp_array;   SWPs_loops.tmr_power{i,1}(sample_tmp+interval)'];
            
            %=== Calculate power spectrum
            [PS_LFP,freq_SG]= cwt(SWPs_loops.LFP{i,1}(sample_tmp+interval)',Fs_Hi);
            spg_array = cat(3,spg_array,abs(PS_LFP));
            [PS_WBT,~]= cwt(SWPs_loops.wbt{i,1}(sample_tmp+interval)',Fs_Hi);
            wbt_array = cat(3,wbt_array,abs(PS_WBT));
            
        end
    end
    
    figure('units','normalized','outerposition',[.3 .05 .15 .9]);
    tiledlayout(7,1,'TileSpacing','compact');
    
    nexttile;   sem = std(clk_array) ./ sqrt(size(clk_array,1));
    plotWinterval_AF_v0(interval./Fs_Hi, mean(clk_array), mean(clk_array)+sem, mean(clk_array)-sem, 'r');
    hold on; plot([0 0], ylim, 'k'); xlim('tight'); ylabel('Echolocation Rate (Hz)'); xticks([]);
    nexttile;   ci = bootci(100, @median, wfr_array);
    plotWinterval_AF_v0(interval./Fs_Hi, median(wfr_array), ci(2,:), ci(1,:), 'k');
    hold on; plot([0 0], ylim, 'k'); xlim('tight'); ylabel('Wingbeat Frequency (Hz)'); xticks([]);
    nexttile;   sem = std(spd_array) ./ sqrt(size(spd_array,1));
    plotWinterval_AF_v0(interval./Fs_Hi, mean(spd_array), mean(spd_array)+sem, mean(spd_array)-sem, 'b');
    hold on; plot([0 0], ylim, 'k'); xlim('tight'); ylabel('Velocity (m/s)'); xticks([]);
    nexttile;   ci = bootci(100, @median, crv_array);   sem = iqr(crv_array) ./ sqrt(size(crv_array,1));
    %     plotWinterval_AF_v0(interval./Fs_Hi, median(crv_array), ci(2,:), ci(1,:), 'g');
    plotWinterval_AF_v0(interval./Fs_Hi, median(crv_array), median(crv_array)+sem, median(crv_array)-sem, 'g');
    hold on; plot([0 0], ylim, 'k'); xlim('tight'); ylabel('Curvature (1/m)'); xticks([]);
    nexttile;   ci = bootci(100, @median, lfp_array);
    plotWinterval_AF_v0(interval./Fs_Hi, median(lfp_array), ci(2,:), ci(1,:), 'm');
    hold on; plot([0 0], ylim, 'k'); xlim('tight'); ylabel('Non-oscillatory Power (\muV^2)'); xticks([]);
    nexttile;   imagesc([-1 1], [freq_SG(1), freq_SG(end)], imgaussfilt(squeeze(mean(spg_array,3)), [0.01 0.1]));
    shading interp; colormap(redblue);  hold on; plot([0 0], ylim, 'w--');  set(gca, 'YScale', 'log', 'YDir', 'normal', 'TickLength', [0 0]);   ylim([2 30]); yticks([2 5 10 20 30]);   ylabel('Freq (Hz)'); title('LFP'); xticks([]);
    nexttile;   imagesc([-1 1], [freq_SG(1), freq_SG(end)], imgaussfilt(squeeze(mean(wbt_array,3)), [0.01 0.1]));
    shading interp; colormap(redblue);  hold on; plot([0 0], ylim, 'w--');  set(gca, 'YScale', 'log', 'YDir', 'normal', 'TickLength', [0 0]);   ylim([2 30]); yticks([2 5 10 20 30]);   ylabel('Freq (Hz)'); xlabel('Time from turn (s)');  title('Accel.');
    sgtitle('Loops');
    
    %=== Look at wingbeat and echolocation around turns
    SWPs_strgth = SWPs(~SWPs.has_turn & SWPs.mic_flag,:);
    clk_array = [];  wfr_array = [];    spd_array = [];    crv_array = [];  lfp_array = [];  spg_array = [];
    for i = 1:size(SWPs_strgth,1)
        
        sample_tmp = round(numel(SWPs_strgth.curv_flag{i,1})/2);
        
        if sample_tmp>interval(end) && sample_tmp+interval(end)<numel(SWPs_strgth.clk{i,1})
            clk_array = [clk_array;   SWPs_strgth.clk{i,1}(sample_tmp+interval)'];
            wfr_array = [wfr_array;   SWPs_strgth.wbt_freq{i,1}(sample_tmp+interval)'];
            spd_array = [spd_array;   SWPs_strgth.speed{i,1}(sample_tmp+interval)'];
            crv_array = [crv_array;   SWPs_strgth.curv3D_2{i,1}(sample_tmp+interval)'];
            lfp_array = [lfp_array;   SWPs_strgth.tmr_power{i,1}(sample_tmp+interval)'];
            
            %=== Calculate power spectrum
            [PS_LFP,freq_SG]= cwt(SWPs_strgth.LFP{i,1}(sample_tmp+interval)',Fs_Hi);
            spg_array = cat(3,spg_array,abs(PS_LFP));
            [PS_WBT,~]= cwt(SWPs_strgth.wbt{i,1}(sample_tmp+interval)',Fs_Hi);
            wbt_array = cat(3,wbt_array,abs(PS_WBT));
            
        end
    end
    
    figure('units','normalized','outerposition',[.5 .05 .15 .9]);
    tiledlayout(7,1,'TileSpacing','compact');
    
    nexttile;   sem = std(clk_array) ./ sqrt(size(clk_array,1));
    plotWinterval_AF_v0(interval./Fs_Hi, mean(clk_array), mean(clk_array)+sem, mean(clk_array)-sem, 'r');
    hold on; plot([0 0], ylim, 'k'); xlim('tight'); ylabel('Echolocation Rate (Hz)'); xticks([]);
    nexttile;   ci = bootci(100, @median, wfr_array);
    plotWinterval_AF_v0(interval./Fs_Hi, median(wfr_array), ci(2,:), ci(1,:), 'k');
    hold on; plot([0 0], ylim, 'k'); xlim('tight'); ylabel('Wingbeat Frequency (Hz)'); xticks([]);
    nexttile;   sem = std(spd_array) ./ sqrt(size(spd_array,1));
    plotWinterval_AF_v0(interval./Fs_Hi, mean(spd_array), mean(spd_array)+sem, mean(spd_array)-sem, 'b');
    hold on; plot([0 0], ylim, 'k'); xlim('tight'); ylabel('Velocity (m/s)'); xticks([]);
    nexttile;   ci = bootci(100, @median, crv_array);   sem = iqr(crv_array) ./ sqrt(size(crv_array,1));
    %     plotWinterval_AF_v0(interval./Fs_Hi, median(crv_array), ci(2,:), ci(1,:), 'g');
    plotWinterval_AF_v0(interval./Fs_Hi, median(crv_array), median(crv_array)+sem, median(crv_array)-sem, 'g');
    hold on; plot([0 0], ylim, 'k'); xlim('tight'); ylabel('Curvature (1/m)'); xticks([]);
    nexttile;   ci = bootci(100, @median, lfp_array);
    plotWinterval_AF_v0(interval./Fs_Hi, median(lfp_array), ci(2,:), ci(1,:), 'm');
    hold on; plot([0 0], ylim, 'k'); xlim('tight'); ylabel('Non-oscillatory Power (\muV^2)'); xticks([]);
    nexttile;   imagesc([-1 1], [freq_SG(1), freq_SG(end)], imgaussfilt(squeeze(mean(spg_array,3)), [0.01 0.1]));
    shading interp; colormap(redblue);  hold on; plot([0 0], ylim, 'w--');  set(gca, 'YScale', 'log', 'YDir', 'normal', 'TickLength', [0 0]);   ylim([2 30]); yticks([2 5 10 20 30]);   ylabel('Freq (Hz)'); title('LFP'); xticks([]);
    nexttile;   imagesc([-1 1], [freq_SG(1), freq_SG(end)], imgaussfilt(squeeze(mean(wbt_array,3)), [0.01 0.1]));
    shading interp; colormap(redblue);  hold on; plot([0 0], ylim, 'w--');  set(gca, 'YScale', 'log', 'YDir', 'normal', 'TickLength', [0 0]);   ylim([2 30]); yticks([2 5 10 20 30]);   ylabel('Freq (Hz)'); xlabel('Time from turn (s)');  title('Accel.');
    sgtitle('Straigth Flights');
    
    %=== Look at echolocation distribution along the flight
    figure('units','normalized','outerposition',[.3 .1 .1 .3]);
    SWPs_loops.echo_phase = cellfun(@(x) interp1(1:numel(x),x,linspace(1,numel(x),100))./max(x),SWPs_loops.clk,'UniformOutput',false);
    SWPs_strgth.echo_phase = cellfun(@(x) interp1(1:numel(x),x,linspace(1,numel(x),100))./max(x),SWPs_strgth.clk,'UniformOutput',false);
    hold on;
    data_tmp = vertcat(SWPs_loops.echo_phase{:});
    plotWinterval_AF_v0(1:100,mean(data_tmp,'omitnan'),mean(data_tmp,'omitnan')+std(data_tmp,'omitnan')./sqrt(size(data_tmp,1)),mean(data_tmp,'omitnan')-std(data_tmp,'omitnan')./sqrt(size(data_tmp,1)),'r');
    data_tmp = vertcat(SWPs_strgth.echo_phase{:});
    plotWinterval_AF_v0(1:100,mean(data_tmp,'omitnan'),mean(data_tmp,'omitnan')+std(data_tmp,'omitnan')./sqrt(size(data_tmp,1)),mean(data_tmp,'omitnan')-std(data_tmp,'omitnan')./sqrt(size(data_tmp,1)),'b');
    xlabel('Flight Phase'); ylabel('Average Echolocation (Norm)');
    
    unique([string(SWPs_loops.unique_ID)],'rows')
    unique([string(SWPs_strgth.unique_ID)],'rows')
    
end

%% THETA-SWEEPS AND WINGBEAT
for hide=1
    
    %=== Select subtable of good flights
    warning('off');
    SWPs_sst = SWPs(SWPs.rmsDec_error<2 & SWPs.prc_decoded>0.5,:);
    
    %=== Params and initialize relevant variables
    bin_size_1D = 0.15;             % Spatial bin size
    n_swp_shuffles = 10;             % Default 20
    n_reshape = 50;                 % Default 50
    smooth_f = [1 .3];              % [space bin,time bin]
    N_flights = size(SWPs_sst,1);   % Number of flights
    max_int_shift = round(0.120/mean(diff(SWPs_sst.bin_time{1,1})));    % Random shift for shuffling (within plus/minus 60 ms)
    single_SWP_sh_cell = {};
    
    %=== Useful to define this structure to speed code up
    template = struct( ...
        'flight_ID', [], 'bat_ID', [], 'session_ID', [], ...
        'rms_dec', [], 'prc_dec', [], 'FromFlightWithClicks', [], 'FromFlightWithTurns',[],...
        'raw_posterior', [], 'sft_posterior', [], 'rsp_posterior', [], ...
        'rsz_spk_dsty', [], 'wbt_power', [], 'rsz_LFP', [], 'rsz_wbt', [], ...
        'fract_pos', [], 'rsz_clk', [],'rsz_trn', [],'rsz_crv', [], 'rsz_wbt_phase', [], 'rsz_tmr_phase', [], ...
        'mean_spk_dsty', [], 'mean_fract_pos', [], ...
        'med_jmp_distance', [], 'mean_jmp_distance', [], 'med_max_post', [], ...
        'est_dist1', [], 'est_dist2', []);
    
    %=== Real data, Cut at wingbeat maxima
    single_SWP_rl = table();
    for zz=1:N_flights
        
        %=== Find wingbeat minima
        zero_phs_idx = find(SWPs_sst.wbt_phase{zz,1}(1:end-1).* SWPs_sst.wbt_phase{zz,1}(2:end)<0 & diff(SWPs_sst.wbt_phase{zz,1})<0);  % Segment based on wingbeat
        sweep_strt = zero_phs_idx(1:end-1); sweep_stop = zero_phs_idx(2:end);
        N_spatial_bins = size(SWPs_sst.p_dec_flight{zz,1},1);
        spt_bin_ids = [1:N_spatial_bins]';
        
        %=== Preallocate a structure
        single_SWP_struct = repmat(template, numel(sweep_strt), 1);
        %=== Real data
        for ss=1:numel(sweep_strt)
            swp_interval = sweep_strt(ss):sweep_stop(ss);
            
            %=== Basic features of the sweep
            single_SWP_struct(ss).flight_ID = zz;
            single_SWP_struct(ss).bat_ID = SWPs_sst.unique_ID(zz,2);
            single_SWP_struct(ss).session_ID = SWPs_sst.unique_ID(zz,3);
            single_SWP_struct(ss).rms_dec = SWPs_sst.rmsDec_error(zz);
            single_SWP_struct(ss).prc_dec = SWPs_sst.prc_decoded(zz);
            single_SWP_struct(ss).FromFlightWithClicks = any(SWPs_sst.clk{zz,1});
            single_SWP_struct(ss).FromFlightWithTurns = any(SWPs_sst.trn{zz,1});
            
            %=== Raw, filtered, shifted and rehaped posterior
            single_SWP_struct(ss).raw_posterior = {SWPs_sst.p_dec_flight{zz,1}(:,swp_interval)};
            single_SWP_struct(ss).sft_posterior = {SWPs_sst.p_dec_shifted{zz,1}(:,swp_interval)};
            single_SWP_struct(ss).rsp_posterior = {imresize(single_SWP_struct(ss).sft_posterior{:},[size(single_SWP_struct(ss).sft_posterior{:},1),n_reshape])};
            
            %=== Spike density, wingbeat, LFP, Tamir's phase echolocation and phase
            single_SWP_struct(ss).rsz_spk_dsty = {interp1(SWPs_sst.spk_dsty{zz,1}(swp_interval),linspace(1,numel(SWPs_sst.spk_dsty{zz,1}(swp_interval)),n_reshape)')};
            single_SWP_struct(ss).wbt_power = {SWPs_sst.wbt_power{zz,1}(swp_interval)};
            single_SWP_struct(ss).rsz_LFP = {interp1(SWPs_sst.LFP{zz,1}(swp_interval),linspace(1,numel(SWPs_sst.LFP{zz,1}(swp_interval)),n_reshape)')};
            single_SWP_struct(ss).rsz_wbt = {interp1(SWPs_sst.wbt{zz,1}(swp_interval),linspace(1,numel(SWPs_sst.wbt{zz,1}(swp_interval)),n_reshape)')};
            single_SWP_struct(ss).fract_pos = {SWPs_sst.pos_real{zz,1}(swp_interval)/SWPs_sst.pos_real{zz,1}(end)};
            single_SWP_struct(ss).rsz_clk = {interp1(SWPs_sst.clk{zz,1}(swp_interval),linspace(1,numel(SWPs_sst.clk{zz,1}(swp_interval)),n_reshape)')};
            single_SWP_struct(ss).rsz_trn = {interp1(SWPs_sst.trn{zz,1}(swp_interval),linspace(1,numel(SWPs_sst.trn{zz,1}(swp_interval)),n_reshape)')};
            single_SWP_struct(ss).rsz_crv = {interp1(SWPs_sst.curv3D_2{zz,1}(swp_interval),linspace(1,numel(SWPs_sst.curv3D_2{zz,1}(swp_interval)),n_reshape)')};
            single_SWP_struct(ss).rsz_wbt_phase = {interp1(SWPs_sst.wbt_phase{zz,1}(swp_interval),linspace(1,numel(SWPs_sst.wbt_phase{zz,1}(swp_interval)),n_reshape)')};
            single_SWP_struct(ss).rsz_tmr_phase = {interp1(SWPs_sst.tmr_phase{zz,1}(swp_interval),linspace(1,numel(SWPs_sst.tmr_phase{zz,1}(swp_interval)),n_reshape)')};
            
            %=== Add some features of the single cycle
            single_SWP_struct(ss).mean_spk_dsty = mean(single_SWP_struct(ss).rsz_spk_dsty{:});   % Average spike density
            single_SWP_struct(ss).mean_fract_pos = mean(single_SWP_struct(ss).fract_pos{:});     % Average phase of the flight
            
            %=== Add some features of the posterior
            [max_p,max_loc] = max(single_SWP_struct(ss).rsp_posterior{:},[],1);        % Location of max posterior
            cnt_mass = spt_bin_ids'*single_SWP_struct(ss).rsp_posterior{:};            % Center of mass
            single_SWP_struct(ss).med_jmp_distance = median(abs(diff(cnt_mass)));     % Median jump distance with center of mass
            single_SWP_struct(ss).mean_jmp_distance = mean(abs(diff(max_loc)));       % Mean jump distance with loc of max posterior
            single_SWP_struct(ss).med_max_post = mean(max_p);                         % Average max posterior
            single_SWP_struct(ss).est_dist1 = {(cnt_mass-N_spatial_bins/2)*bin_size_1D};  % Decoding error, reshaped and in m (center of mass)
            single_SWP_struct(ss).est_dist2 = {(max_loc-N_spatial_bins/2)*bin_size_1D};   % Decoding error, reshaped and in m (max posterior)
            
        end
        single_SWP_rl = [single_SWP_rl; struct2table(single_SWP_struct)];
        
    end
    
    %=== Shuffled data, Cut at random points
    for jjj=1:n_swp_shuffles
        single_SWP_sh = table();   counter_sh=1;
        for zz=1:N_flights
            
            %=== Find wingbeat minima
            zero_phs_idx = find(SWPs_sst.wbt_phase{zz,1}(1:end-1).* SWPs_sst.wbt_phase{zz,1}(2:end)<0 & diff(SWPs_sst.wbt_phase{zz,1})<0);  % Segment based on wingbeat
            sweep_strt = zero_phs_idx(1:end-1); sweep_stop = zero_phs_idx(2:end);
            N_spatial_bins = size(SWPs_sst.p_dec_flight{zz,1},1);
            spt_bin_ids = [1:N_spatial_bins]';
            
            %=== Shuffling procedure
            rand_shift =  -round(max_int_shift/2)+randi(max_int_shift,numel(sweep_strt),1);
            sweep_strt_sh = sweep_strt+rand_shift; sweep_stop_sh = sweep_stop+rand_shift;
            sweep_strt = sweep_strt_sh(sweep_strt_sh>1 & sweep_stop_sh<size(SWPs_sst.wbt_phase{zz,1},1));
            sweep_stop = sweep_stop_sh(sweep_strt_sh>1 & sweep_stop_sh<size(SWPs_sst.wbt_phase{zz,1},1));
            
            %=== Preallocate a structure
            single_SWP_struct = repmat(template, numel(sweep_strt), 1);
            %=== Shuffle data
            for ss=1:numel(sweep_strt)
                
                swp_interval = sweep_strt(ss):sweep_stop(ss);
                
                %=== Basic features of the sweep
                single_SWP_struct(ss).flight_ID = zz;
                single_SWP_struct(ss).bat_ID = SWPs_sst.unique_ID(zz,2);
                single_SWP_struct(ss).session_ID = SWPs_sst.unique_ID(zz,3);
                single_SWP_struct(ss).rms_dec = SWPs_sst.rmsDec_error(zz);
                single_SWP_struct(ss).prc_dec = SWPs_sst.prc_decoded(zz);
                single_SWP_struct(ss).FromFlightWithClicks = any(SWPs_sst.clk{zz,1});
                
                %=== Raw, filtered, shifted and rehaped posterior
                single_SWP_struct(ss).raw_posterior = {SWPs_sst.p_dec_flight{zz,1}(:,swp_interval)};
                %single_SWP_struct(ss).sft_posterior(counter_sh) = {imgaussfilt(SWPs_sst.p_dec_shifted{zz,1}(:,swp_interval),smooth_f)};
                single_SWP_struct(ss).sft_posterior = {SWPs_sst.p_dec_shifted{zz,1}(:,swp_interval)};
                single_SWP_struct(ss).rsp_posterior = {imresize(single_SWP_struct(ss).sft_posterior{:},[size(single_SWP_struct(ss).sft_posterior{:},1),n_reshape])};
                
                %=== Spike density, wingbeat, LFP, Tamir's phase echolocation and phase
                single_SWP_struct(ss).rsz_spk_dsty = {interp1(SWPs_sst.spk_dsty{zz,1}(swp_interval),linspace(1,numel(SWPs_sst.spk_dsty{zz,1}(swp_interval)),n_reshape)')};
                single_SWP_struct(ss).wbt_power = {SWPs_sst.wbt_power{zz,1}(swp_interval)};
                single_SWP_struct(ss).rsz_LFP = {interp1(SWPs_sst.LFP{zz,1}(swp_interval),linspace(1,numel(SWPs_sst.LFP{zz,1}(swp_interval)),n_reshape)')};
                single_SWP_struct(ss).rsz_wbt = {interp1(SWPs_sst.wbt{zz,1}(swp_interval),linspace(1,numel(SWPs_sst.wbt{zz,1}(swp_interval)),n_reshape)')};
                single_SWP_struct(ss).fract_pos = {SWPs_sst.pos_real{zz,1}(swp_interval)/SWPs_sst.pos_real{zz,1}(end)};
                single_SWP_struct(ss).rsz_clk = {interp1(SWPs_sst.clk{zz,1}(swp_interval),linspace(1,numel(SWPs_sst.clk{zz,1}(swp_interval)),n_reshape)')};
                single_SWP_struct(ss).rsz_wbt_phase = {interp1(SWPs_sst.wbt_phase{zz,1}(swp_interval),linspace(1,numel(SWPs_sst.wbt_phase{zz,1}(swp_interval)),n_reshape)')};
                single_SWP_struct(ss).rsz_tmr_phase = {interp1(SWPs_sst.tmr_phase{zz,1}(swp_interval),linspace(1,numel(SWPs_sst.tmr_phase{zz,1}(swp_interval)),n_reshape)')};
                
                %=== Add some features of the single cycle
                single_SWP_struct(ss).mean_spk_dsty = mean(single_SWP_struct(ss).rsz_spk_dsty{:});   % Average spike density
                single_SWP_struct(ss).mean_fract_pos = mean(single_SWP_struct(ss).fract_pos{:}); % Average phase of the flight
                
                %=== Add some features of the posterior
                [max_p,max_loc] = max(single_SWP_struct(ss).rsp_posterior{:},[],1);        % Location of max posterior
                cnt_mass = spt_bin_ids'*single_SWP_struct(ss).rsp_posterior{:};            % Center of mass
                single_SWP_struct(ss).med_jmp_distance = median(abs(diff(cnt_mass)));     % Median jump distance with center of mass
                single_SWP_struct(ss).mean_jmp_distance = mean(abs(diff(max_loc)));       % Mean jump distance with loc of max posterior
                single_SWP_struct(ss).med_max_post = mean(max_p);                         % Average max posterior
                single_SWP_struct(ss).est_dist1 = {(cnt_mass-N_spatial_bins/2)*bin_size_1D};  % Decoding error, reshaped and in m (center of mass)
                single_SWP_struct(ss).est_dist2 = {(max_loc-N_spatial_bins/2)*bin_size_1D};   % Decoding error, reshaped and in m (max posterior)
                
            end
            
            single_SWP_sh = [single_SWP_sh; struct2table(single_SWP_struct)];
            
        end
        single_SWP_sh_cell = [single_SWP_sh_cell;{single_SWP_sh}];
    end
    warning('on');
    
    %%
    %=== Extract subset using defined criteria (exclude flight tails, epochs of low firing and flat sweeps)
    min_pos1 = 0.15;
    min_pos2 = 0.85;
    min_rms = 1.3;
    min_prc = 0.7;
    min_mjp = 0.0;
    min_mp = 0.0;
    
    SWP_sst_rl = single_SWP_rl(single_SWP_rl.mean_fract_pos>min_pos1 & single_SWP_rl.mean_fract_pos<min_pos2 & single_SWP_rl.rms_dec<min_rms &...
        single_SWP_rl.prc_dec>min_prc & single_SWP_rl.mean_jmp_distance>min_mjp & single_SWP_rl.med_max_post>min_mp,:);
    N_sweeps_rl = size(SWP_sst_rl,1);
    
    %=== Calculate averages and STDs (REAL DATA)
    est_dist_rl = zeros(n_reshape,N_sweeps_rl);
    rsz_spk_dsty_rl = zeros(n_reshape,N_sweeps_rl);
    rsz_wbt_rl = zeros(n_reshape,N_sweeps_rl);
    rsz_LFP_rl = zeros(n_reshape,N_sweeps_rl);
    rsz_clk_rl = zeros(n_reshape,N_sweeps_rl);
    rsz_wbt_phase_rl = zeros(n_reshape,N_sweeps_rl);
    for i=1:N_sweeps_rl
        est_dist_rl(:,i) = SWP_sst_rl.est_dist1{i,1};
        rsz_spk_dsty_rl(:,i) = SWP_sst_rl.rsz_spk_dsty{i,1};
        rsz_wbt_rl(:,i) = SWP_sst_rl.rsz_wbt{i,1};
        rsz_LFP_rl(:,i) = SWP_sst_rl.rsz_LFP{i,1};
        rsz_clk_rl(:,i) = SWP_sst_rl.rsz_clk{i,1};
        rsz_wbt_phase_rl(:,i) = SWP_sst_rl.rsz_wbt_phase{i,1};
    end
    avg_est_dist_rl = mean(est_dist_rl,2);            sem_est_dist_rl = std(est_dist_rl,[],2)/sqrt(N_sweeps_rl);
    avg_rsz_spk_dsty_rl = mean(rsz_spk_dsty_rl,2);    sem_rsz_spk_dsty_rl = std(rsz_spk_dsty_rl,[],2)/sqrt(N_sweeps_rl);
    avg_rsz_wbt_rl = mean(rsz_wbt_rl,2);              sem_rsz_wbt_rl = std(rsz_wbt_rl,[],2)/sqrt(N_sweeps_rl);
    avg_rsz_LFP_rl = mean(rsz_LFP_rl,2);              sem_rsz_LFP_rl = std(rsz_LFP_rl,[],2)/sqrt(N_sweeps_rl);
    avg_rsz_wbt_phase_rl = mean(rsz_wbt_phase_rl,2);  sem_rsz_wbt_phase_rl = std(rsz_wbt_phase_rl,[],2)/sqrt(N_sweeps_rl);
    avg_rsz_clk_rl = mean(rsz_clk_rl,2);              sem_rsz_clk_rl = std(rsz_clk_rl,[],2)/sqrt(N_sweeps_rl);
    
    %=== Calculate averages and STDs (SHUFFLED DATA)
    avg_est_dist_sh_all = zeros(n_reshape,n_swp_shuffles);
    avg_rsz_spk_dsty_sh_all = zeros(n_reshape,n_swp_shuffles);
    avg_rsz_wbt_sh_all = zeros(n_reshape,n_swp_shuffles);
    avg_rsz_LFP_sh_all = zeros(n_reshape,n_swp_shuffles);
    avg_rsz_wbt_phase_sh_all = zeros(n_reshape,n_swp_shuffles);
    avg_rsz_clk_sh_all = zeros(n_reshape,n_swp_shuffles);
    for jjj=1:n_swp_shuffles
        
        single_SWP_sh = single_SWP_sh_cell{jjj};
        SWP_sst_sh = single_SWP_sh(single_SWP_sh.mean_fract_pos>min_pos1 & single_SWP_sh.mean_fract_pos<min_pos2 & single_SWP_sh.rms_dec<min_rms &...
            single_SWP_sh.prc_dec>min_prc & single_SWP_sh.mean_jmp_distance>min_mjp & single_SWP_sh.med_max_post>min_mp & single_SWP_sh.FromFlightWithClicks,:);
        N_sweeps_sh = size(SWP_sst_sh,1);
        
        est_dist_sh = zeros(n_reshape,N_sweeps_sh);
        rsz_spk_dsty_sh = zeros(n_reshape,N_sweeps_sh);
        rsz_wbt_sh = zeros(n_reshape,N_sweeps_sh);
        rsz_LFP_sh = zeros(n_reshape,N_sweeps_sh);
        rsz_clk_sh = zeros(n_reshape,N_sweeps_sh);
        rsz_wbt_phase_sh = zeros(n_reshape,N_sweeps_sh);
        for i=1:N_sweeps_sh
            est_dist_sh(:,i) = SWP_sst_sh.est_dist1{i,1};
            rsz_spk_dsty_sh(:,i) = SWP_sst_sh.rsz_spk_dsty{i,1};
            rsz_wbt_sh(:,i) = SWP_sst_sh.rsz_wbt{i,1};
            rsz_LFP_sh(:,i) = SWP_sst_sh.rsz_LFP{i,1};
            rsz_clk_sh(:,i) = SWP_sst_sh.rsz_clk{i,1};
            rsz_wbt_phase_sh(:,i) = SWP_sst_sh.rsz_wbt_phase{i,1};
        end
        avg_est_dist_sh_all(:,jjj) = mean(est_dist_sh,2);
        avg_rsz_spk_dsty_sh_all(:,jjj) = mean(rsz_spk_dsty_sh,2);
        avg_rsz_wbt_sh_all(:,jjj) = mean(rsz_wbt_sh,2);
        avg_rsz_LFP_sh_all(:,jjj) = mean(rsz_LFP_sh,2);
        avg_rsz_wbt_phase_sh_all(:,jjj) = mean(rsz_wbt_phase_sh,2);
        avg_rsz_clk_sh_all(:,jjj) = mean(rsz_clk_sh,2);
    end
    
    %=== Average shuffles
    avg_est_dist_sh = mean(avg_est_dist_sh_all,2);            sem_est_dist_sh = std(avg_est_dist_sh_all,[],2)/sqrt(n_swp_shuffles);
    avg_rsz_spk_dsty_sh = mean(avg_rsz_spk_dsty_sh_all,2);    sem_rsz_spk_dsty_sh = std(avg_rsz_spk_dsty_sh_all,[],2)/sqrt(n_swp_shuffles);
    avg_rsz_wbt_sh = mean(avg_rsz_wbt_sh_all,2);              sem_rsz_wbt_sh = std(avg_rsz_wbt_sh_all,[],2)/sqrt(n_swp_shuffles);
    avg_rsz_LFP_sh = mean(avg_rsz_LFP_sh_all,2);              sem_rsz_LFP_sh = std(avg_rsz_LFP_sh_all,[],2)/sqrt(n_swp_shuffles);
    avg_rsz_wbt_phase_sh = mean(avg_rsz_wbt_phase_sh_all,2);  sem_rsz_wbt_phase_sh = std(avg_rsz_wbt_phase_sh_all,[],2)/sqrt(n_swp_shuffles);
    avg_rsz_clk_sh = mean(avg_rsz_clk_sh_all,2);              sem_rsz_clk_sh = std(avg_rsz_clk_sh_all,[],2)/sqrt(n_swp_shuffles);
    
    %=== Define bins and find peak decoding error phase
    rsz_bin_ctrs = linspace(-180,180,n_reshape);
    [max_val,max_loc] = max(avg_est_dist_rl);                rsz_bin_ctrs(max_loc);
    sum(abs(diff(SWP_sst_rl.flight_ID)))+1;
    
    %=== Calculate significance trace
    pVal_est_dist = sum(avg_est_dist_rl<avg_est_dist_sh_all,2)./n_swp_shuffles;
    
    %=== Show averages
    figure('units','normalized','outerposition',[.3 .3 .45 .3]);
    tiledlayout(1,6,'TileSpacing','compact');
    nexttile;   hold on;
    plotWinterval_AF_v0(rsz_bin_ctrs,avg_est_dist_sh,avg_est_dist_sh-sem_est_dist_sh,avg_est_dist_sh+sem_est_dist_sh,'k');
    plotWinterval_AF_v0(rsz_bin_ctrs,avg_est_dist_rl,avg_est_dist_rl-sem_est_dist_rl,avg_est_dist_rl+sem_est_dist_rl,[.9 .5 .2]);
    plot(rsz_bin_ctrs(pVal_est_dist<0.05),max(avg_est_dist_rl)*ones(size(rsz_bin_ctrs(pVal_est_dist<0.05))),'*');
    xticks([-170 0 170]);    title('Average Decoding Error');    plot(0*[1 1],ylim,'k--'); xlim('tight');   xlim([-170; 170]);
    nexttile;   hold on;
    plotWinterval_AF_v0(rsz_bin_ctrs,avg_rsz_spk_dsty_sh,avg_rsz_spk_dsty_sh-sem_rsz_spk_dsty_sh,avg_rsz_spk_dsty_sh+sem_rsz_spk_dsty_sh,'k');
    plotWinterval_AF_v0(rsz_bin_ctrs,avg_rsz_spk_dsty_rl,avg_rsz_spk_dsty_rl-sem_rsz_spk_dsty_rl,avg_rsz_spk_dsty_rl+sem_rsz_spk_dsty_rl,[.9 .5 .2]);
    xticks([-170 0 170]);    title('Average Spike Density');    plot(0*[1 1],ylim,'k--'); xlim('tight');    xlim([-170; 170]);
    nexttile;   hold on;
    plotWinterval_AF_v0(rsz_bin_ctrs,avg_rsz_LFP_sh,avg_rsz_LFP_sh-sem_rsz_LFP_sh,avg_rsz_LFP_sh+sem_rsz_LFP_sh,'k');
    plotWinterval_AF_v0(rsz_bin_ctrs,avg_rsz_LFP_rl,avg_rsz_LFP_rl-sem_rsz_LFP_rl,avg_rsz_LFP_rl+sem_rsz_LFP_rl,[.9 .5 .2]);
    xticks([-170 0 170]);    title('Average LFP');    plot(0*[1 1],ylim,'k--'); xlim('tight');              xlim([-170; 170]);
    nexttile;   hold on;
    plotWinterval_AF_v0(rsz_bin_ctrs,avg_rsz_clk_sh,avg_rsz_clk_sh-sem_rsz_clk_sh,avg_rsz_clk_sh+sem_rsz_clk_sh,'k');
    plotWinterval_AF_v0(rsz_bin_ctrs,avg_rsz_clk_rl,avg_rsz_clk_rl-sem_rsz_clk_rl,avg_rsz_clk_rl+sem_rsz_clk_rl,[.9 .5 .2]);
    xticks([-170 0 170]);    title('Average Call Rate');    plot(0*[1 1],ylim,'k--'); xlim('tight');        xlim([-170; 170]);
    nexttile;   hold on;
    plotWinterval_AF_v0(rsz_bin_ctrs,avg_rsz_wbt_phase_sh,avg_rsz_wbt_phase_sh-sem_rsz_wbt_phase_sh,avg_rsz_wbt_phase_sh+sem_rsz_wbt_phase_sh,'k');
    plotWinterval_AF_v0(rsz_bin_ctrs,avg_rsz_wbt_phase_rl,avg_rsz_wbt_phase_rl-sem_rsz_wbt_phase_rl,avg_rsz_wbt_phase_rl+sem_rsz_wbt_phase_rl,[.9 .5 .2]);
    xticks([-170 0 170]);    title('Wingbeat Phase');    plot(0*[1 1],ylim,'k--'); xlim('tight');           xlim([-170; 170]);
    nexttile;   hold on;
    plotWinterval_AF_v0(rsz_bin_ctrs,avg_rsz_wbt_sh,avg_rsz_wbt_sh-sem_rsz_wbt_sh,avg_rsz_wbt_sh+sem_rsz_wbt_sh,'k');
    plotWinterval_AF_v0(rsz_bin_ctrs,avg_rsz_wbt_rl,avg_rsz_wbt_rl-sem_rsz_wbt_rl,avg_rsz_wbt_rl+sem_rsz_wbt_rl,[.9 .5 .2]);
    xticks([-170 0 170]);    title('Accelerometer');    plot(0*[1 1],ylim,'k--'); xlim('tight');           xlim([-170; 170]);
    sgtitle(['Average of ',num2str(N_sweeps_rl),' cycles, from ',num2str(numel(unique(SWP_sst_rl.flight_ID))),' flights']);
    
    %% ECHOLOCATION ANALYSIS ON FLIGHTS WITHOUT TURNS
    
    %=== Extract subset using defined criteria (exclude flight tails, epochs of low firing and flat sweeps)
    %min_pos1 = 0.00;        % Extended to include beginning of a flight, where echolocation is typically happening
    %min_pos2 = 1.00;        % Extended to include end of a flight, where echolocation is typically happening
    min_pos1 = 0.00;        % Extended to include beginning of a flight, where echolocation is typically happening
    min_pos2 = 0.5;        % Extended to include end of a flight, where echolocation is typically happening
    
    min_rms = 1.3;
    min_prc = 0.7;
    min_mjp = 0.0;
    min_mp = 0.0;
    match_samples = 0;
    
    SWP_sst_rlWOEcho = single_SWP_rl(single_SWP_rl.mean_fract_pos>min_pos1 & single_SWP_rl.mean_fract_pos<min_pos2 & single_SWP_rl.rms_dec<min_rms & ~cellfun(@any,single_SWP_rl.rsz_clk) &...
        single_SWP_rl.prc_dec>min_prc & single_SWP_rl.mean_jmp_distance>min_mjp & single_SWP_rl.med_max_post>min_mp & single_SWP_rl.FromFlightWithClicks & ~single_SWP_rl.FromFlightWithTurns,:);
    N_sweeps_rlWOEcho = size(SWP_sst_rlWOEcho, 1);
    
    SWP_sst_rlWIEcho = single_SWP_rl(single_SWP_rl.mean_fract_pos>min_pos1 & single_SWP_rl.mean_fract_pos<min_pos2 & single_SWP_rl.rms_dec<min_rms & cellfun(@any,single_SWP_rl.rsz_clk) &...
        single_SWP_rl.prc_dec>min_prc & single_SWP_rl.mean_jmp_distance>min_mjp & single_SWP_rl.med_max_post>min_mp & single_SWP_rl.FromFlightWithClicks & ~single_SWP_rl.FromFlightWithTurns,:);
    N_sweeps_rlWIEcho = size(SWP_sst_rlWIEcho,1);
    
    if match_samples
        SWP_sst_rlWOEcho = SWP_sst_rlWOEcho(datasample(1:N_sweeps_rlWOEcho,N_sweeps_rlWIEcho,'replace',false),:);
        N_sweeps_rlWOEcho = N_sweeps_rlWIEcho;
    end
    
    unique([string(SWP_sst_rlWOEcho.flight_ID)],'rows')
    unique([string(SWP_sst_rlWOEcho.session_ID),string(SWP_sst_rlWOEcho.bat_ID)],'rows')

    
    %=== Calculate averages and STDs (REAL DATA, WITH ECHO)
    est_dist_rlWIEcho = zeros(n_reshape, N_sweeps_rlWIEcho);
    rsz_spk_dsty_rlWIEcho = zeros(n_reshape, N_sweeps_rlWIEcho);
    rsz_wbt_rlWIEcho = zeros(n_reshape, N_sweeps_rlWIEcho);
    rsz_LFP_rlWIEcho = zeros(n_reshape, N_sweeps_rlWIEcho);
    rsz_clk_rlWIEcho = zeros(n_reshape, N_sweeps_rlWIEcho);
    rsz_wbt_phase_rlWIEcho = zeros(n_reshape, N_sweeps_rlWIEcho);
    
    for i = 1:N_sweeps_rlWIEcho
        est_dist_rlWIEcho(:,i)         = SWP_sst_rlWIEcho.est_dist1{i,1};
        rsz_spk_dsty_rlWIEcho(:,i)     = SWP_sst_rlWIEcho.rsz_spk_dsty{i,1};
        rsz_wbt_rlWIEcho(:,i)          = SWP_sst_rlWIEcho.rsz_wbt{i,1};
        rsz_LFP_rlWIEcho(:,i)          = SWP_sst_rlWIEcho.rsz_LFP{i,1};
        rsz_clk_rlWIEcho(:,i)          = SWP_sst_rlWIEcho.rsz_clk{i,1};
        rsz_wbt_phase_rlWIEcho(:,i)    = SWP_sst_rlWIEcho.rsz_wbt_phase{i,1};
    end
    
    avg_est_dist_rlWIEcho        = mean(est_dist_rlWIEcho, 2);         sem_est_dist_rlWIEcho        = std(est_dist_rlWIEcho, [], 2) / sqrt(N_sweeps_rlWIEcho);
    avg_rsz_spk_dsty_rlWIEcho    = mean(rsz_spk_dsty_rlWIEcho, 2);     sem_rsz_spk_dsty_rlWIEcho    = std(rsz_spk_dsty_rlWIEcho, [], 2) / sqrt(N_sweeps_rlWIEcho);
    avg_rsz_wbt_rlWIEcho         = mean(rsz_wbt_rlWIEcho, 2);          sem_rsz_wbt_rlWIEcho         = std(rsz_wbt_rlWIEcho, [], 2) / sqrt(N_sweeps_rlWIEcho);
    avg_rsz_LFP_rlWIEcho         = mean(rsz_LFP_rlWIEcho, 2);          sem_rsz_LFP_rlWIEcho         = std(rsz_LFP_rlWIEcho, [], 2) / sqrt(N_sweeps_rlWIEcho);
    avg_rsz_wbt_phase_rlWIEcho   = mean(rsz_wbt_phase_rlWIEcho, 2);    sem_rsz_wbt_phase_rlWIEcho   = std(rsz_wbt_phase_rlWIEcho, [], 2) / sqrt(N_sweeps_rlWIEcho);
    avg_rsz_clk_rlWIEcho         = mean(rsz_clk_rlWIEcho, 2);          sem_rsz_clk_rlWIEcho         = std(rsz_clk_rlWIEcho, [], 2) / sqrt(N_sweeps_rlWIEcho);
    
    %=== Calculate averages and STDs (REAL DATA, WITHOUT ECHO)
    est_dist_rlWOEcho         = zeros(n_reshape, N_sweeps_rlWOEcho);
    rsz_spk_dsty_rlWOEcho     = zeros(n_reshape, N_sweeps_rlWOEcho);
    rsz_wbt_rlWOEcho          = zeros(n_reshape, N_sweeps_rlWOEcho);
    rsz_LFP_rlWOEcho          = zeros(n_reshape, N_sweeps_rlWOEcho);
    rsz_clk_rlWOEcho          = zeros(n_reshape, N_sweeps_rlWOEcho);
    rsz_wbt_phase_rlWOEcho    = zeros(n_reshape, N_sweeps_rlWOEcho);
    
    for i = 1:N_sweeps_rlWOEcho
        est_dist_rlWOEcho(:,i)         = SWP_sst_rlWOEcho.est_dist1{i,1};
        rsz_spk_dsty_rlWOEcho(:,i)     = SWP_sst_rlWOEcho.rsz_spk_dsty{i,1};
        rsz_wbt_rlWOEcho(:,i)          = SWP_sst_rlWOEcho.rsz_wbt{i,1};
        rsz_LFP_rlWOEcho(:,i)          = SWP_sst_rlWOEcho.rsz_LFP{i,1};
        rsz_clk_rlWOEcho(:,i)          = SWP_sst_rlWOEcho.rsz_clk{i,1};
        rsz_wbt_phase_rlWOEcho(:,i)    = SWP_sst_rlWOEcho.rsz_wbt_phase{i,1};
    end
    
    avg_est_dist_rlWOEcho        = mean(est_dist_rlWOEcho, 2);         sem_est_dist_rlWOEcho        = std(est_dist_rlWOEcho, [], 2) / sqrt(N_sweeps_rlWOEcho);
    avg_rsz_spk_dsty_rlWOEcho    = mean(rsz_spk_dsty_rlWOEcho, 2);     sem_rsz_spk_dsty_rlWOEcho    = std(rsz_spk_dsty_rlWOEcho, [], 2) / sqrt(N_sweeps_rlWOEcho);
    avg_rsz_wbt_rlWOEcho         = mean(rsz_wbt_rlWOEcho, 2);          sem_rsz_wbt_rlWOEcho         = std(rsz_wbt_rlWOEcho, [], 2) / sqrt(N_sweeps_rlWOEcho);
    avg_rsz_LFP_rlWOEcho         = mean(rsz_LFP_rlWOEcho, 2);          sem_rsz_LFP_rlWOEcho         = std(rsz_LFP_rlWOEcho, [], 2) / sqrt(N_sweeps_rlWOEcho);
    avg_rsz_wbt_phase_rlWOEcho   = mean(rsz_wbt_phase_rlWOEcho, 2);    sem_rsz_wbt_phase_rlWOEcho   = std(rsz_wbt_phase_rlWOEcho, [], 2) / sqrt(N_sweeps_rlWOEcho);
    avg_rsz_clk_rlWOEcho         = mean(rsz_clk_rlWOEcho, 2);          sem_rsz_clk_rlWOEcho         = std(rsz_clk_rlWOEcho, [], 2) / sqrt(N_sweeps_rlWOEcho);
    
    %=== Plot comparison: WITH ECHO vs WITHOUT ECHO
    figure('units','normalized','outerposition',[.3 .3 .45 .3]);
    tiledlayout(1,6,'TileSpacing','compact');
    nexttile; hold on;
    plotWinterval_AF_v0(rsz_bin_ctrs, avg_est_dist_rlWOEcho-min(avg_est_dist_rlWOEcho), avg_est_dist_rlWOEcho -min(avg_est_dist_rlWOEcho) - sem_est_dist_rlWOEcho, avg_est_dist_rlWOEcho -min(avg_est_dist_rlWOEcho) + sem_est_dist_rlWOEcho, 'r');
    plotWinterval_AF_v0(rsz_bin_ctrs, avg_est_dist_rlWIEcho-min(avg_est_dist_rlWIEcho), avg_est_dist_rlWIEcho -min(avg_est_dist_rlWIEcho) - sem_est_dist_rlWIEcho, avg_est_dist_rlWIEcho -min(avg_est_dist_rlWIEcho) + sem_est_dist_rlWIEcho, 'b');
    xticks([-170 0 170]); title('Decoding Error'); plot([0 0], ylim, 'k--'); xlim([-170 170]);  legend('Without','','With');    ylim([-0.05 0.15]);
    nexttile; hold on;
    plotWinterval_AF_v0(rsz_bin_ctrs, avg_rsz_spk_dsty_rlWOEcho, avg_rsz_spk_dsty_rlWOEcho - sem_rsz_spk_dsty_rlWOEcho, avg_rsz_spk_dsty_rlWOEcho + sem_rsz_spk_dsty_rlWOEcho, 'r');
    plotWinterval_AF_v0(rsz_bin_ctrs, avg_rsz_spk_dsty_rlWIEcho, avg_rsz_spk_dsty_rlWIEcho - sem_rsz_spk_dsty_rlWIEcho, avg_rsz_spk_dsty_rlWIEcho + sem_rsz_spk_dsty_rlWIEcho, 'b');
    xticks([-170 0 170]); title('Spike Density'); plot([0 0], ylim, 'k--'); xlim([-170 170]);
    nexttile; hold on;
    plotWinterval_AF_v0(rsz_bin_ctrs, avg_rsz_LFP_rlWOEcho, avg_rsz_LFP_rlWOEcho - sem_rsz_LFP_rlWOEcho, avg_rsz_LFP_rlWOEcho + sem_rsz_LFP_rlWOEcho, 'r');
    plotWinterval_AF_v0(rsz_bin_ctrs, avg_rsz_LFP_rlWIEcho, avg_rsz_LFP_rlWIEcho - sem_rsz_LFP_rlWIEcho, avg_rsz_LFP_rlWIEcho + sem_rsz_LFP_rlWIEcho, 'b');
    xticks([-170 0 170]); title('LFP'); plot([0 0], ylim, 'k--'); xlim([-170 170]);
    nexttile; hold on;
    plotWinterval_AF_v0(rsz_bin_ctrs, avg_rsz_clk_rlWOEcho, avg_rsz_clk_rlWOEcho - sem_rsz_clk_rlWOEcho, avg_rsz_clk_rlWOEcho + sem_rsz_clk_rlWOEcho, 'r');
    plotWinterval_AF_v0(rsz_bin_ctrs, avg_rsz_clk_rlWIEcho, avg_rsz_clk_rlWIEcho - sem_rsz_clk_rlWIEcho, avg_rsz_clk_rlWIEcho + sem_rsz_clk_rlWIEcho, 'b');
    xticks([-170 0 170]); title('Call Rate'); plot([0 0], ylim, 'k--'); xlim([-170 170]);   ylim([0 0.6]);
    nexttile; hold on;
    plotWinterval_AF_v0(rsz_bin_ctrs, avg_rsz_wbt_phase_rlWOEcho, avg_rsz_wbt_phase_rlWOEcho - sem_rsz_wbt_phase_rlWOEcho, avg_rsz_wbt_phase_rlWOEcho + sem_rsz_wbt_phase_rlWOEcho, 'r');
    plotWinterval_AF_v0(rsz_bin_ctrs, avg_rsz_wbt_phase_rlWIEcho, avg_rsz_wbt_phase_rlWIEcho - sem_rsz_wbt_phase_rlWIEcho, avg_rsz_wbt_phase_rlWIEcho + sem_rsz_wbt_phase_rlWIEcho, 'b');
    xticks([-170 0 170]); title('Wingbeat Phase'); plot([0 0], ylim, 'k--'); xlim([-170 170]);
    nexttile; hold on;
    plotWinterval_AF_v0(rsz_bin_ctrs, avg_rsz_wbt_rlWOEcho, avg_rsz_wbt_rlWOEcho - sem_rsz_wbt_rlWOEcho, avg_rsz_wbt_rlWOEcho + sem_rsz_wbt_rlWOEcho, 'r');
    plotWinterval_AF_v0(rsz_bin_ctrs, avg_rsz_wbt_rlWIEcho, avg_rsz_wbt_rlWIEcho - sem_rsz_wbt_rlWIEcho, avg_rsz_wbt_rlWIEcho + sem_rsz_wbt_rlWIEcho, 'b');
    xticks([-170 0 170]); title('Accelerometer'); plot([0 0], ylim, 'k--'); xlim([-170 170]);
    sgtitle(['Avg of ', num2str(N_sweeps_rlWIEcho), ' (W/ clicks) vs ', num2str(N_sweeps_rlWOEcho), ' (W/O clicks) sweeps']);
    
    %=== Plot decoding error on top of call rate
    figure('units','normalized','outerposition',[.3 .1 .15 .3]);
    hold on;
    plot(rsz_bin_ctrs, normalize(avg_est_dist_rlWIEcho,'range'));
    plot(rsz_bin_ctrs, normalize(avg_rsz_clk_rlWIEcho,'range'));
    xticks([-170 0 170]); title('Accelerometer'); plot([0 0], ylim, 'k--'); xlim([-170 170]);
    
    %% TURN ANALYSIS
    %=== Extract subset using defined criteria (exclude flight tails, epochs of low firing and flat sweeps)
    min_pos1 = 0.05;
    min_pos2 = 0.95;
    min_rms = 1.3;
    min_prc = 0.7;
    min_mjp = 0.0;
    min_mp = 0.0;
    match_samples = 0;
    
    SWP_sst_rlWOTurn = single_SWP_rl(single_SWP_rl.mean_fract_pos>min_pos1 & single_SWP_rl.mean_fract_pos<min_pos2 & single_SWP_rl.rms_dec<min_rms & ~cellfun(@any,single_SWP_rl.rsz_trn) &...
        single_SWP_rl.prc_dec>min_prc & single_SWP_rl.mean_jmp_distance>min_mjp & single_SWP_rl.med_max_post>min_mp & single_SWP_rl.FromFlightWithTurns,:);
    N_sweeps_rlWOTurn = size(SWP_sst_rlWOTurn, 1);
    
    SWP_sst_rlWITurn = single_SWP_rl(single_SWP_rl.mean_fract_pos>min_pos1 & single_SWP_rl.mean_fract_pos<min_pos2 & single_SWP_rl.rms_dec<min_rms & cellfun(@any,single_SWP_rl.rsz_trn) &...
        single_SWP_rl.prc_dec>min_prc & single_SWP_rl.mean_jmp_distance>min_mjp & single_SWP_rl.med_max_post>min_mp & single_SWP_rl.FromFlightWithTurns,:);
    N_sweeps_rlWITurn = size(SWP_sst_rlWITurn,1);
    
    if match_samples
        SWP_sst_rlWOTurn = SWP_sst_rlWOTurn(datasample(1:N_sweeps_rlWOTurn,N_sweeps_rlWITurn,'replace',false),:);
        N_sweeps_rlWOTurn = N_sweeps_rlWITurn;
    end
    
    %=== Calculate averages and STDs (REAL DATA, WITH TURN)
    est_dist_rlWITurn = zeros(n_reshape, N_sweeps_rlWITurn);
    rsz_spk_dsty_rlWITurn = zeros(n_reshape, N_sweeps_rlWITurn);
    rsz_wbt_rlWITurn = zeros(n_reshape, N_sweeps_rlWITurn);
    rsz_LFP_rlWITurn = zeros(n_reshape, N_sweeps_rlWITurn);
    rsz_clk_rlWITurn = zeros(n_reshape, N_sweeps_rlWITurn);
    rsz_crv_rlWITurn = zeros(n_reshape, N_sweeps_rlWITurn);
    rsz_wbt_phase_rlWITurn = zeros(n_reshape, N_sweeps_rlWITurn);
    
    for i = 1:N_sweeps_rlWITurn
        est_dist_rlWITurn(:,i)         = SWP_sst_rlWITurn.est_dist1{i,1};
        rsz_spk_dsty_rlWITurn(:,i)     = SWP_sst_rlWITurn.rsz_spk_dsty{i,1};
        rsz_wbt_rlWITurn(:,i)          = SWP_sst_rlWITurn.rsz_wbt{i,1};
        rsz_LFP_rlWITurn(:,i)          = SWP_sst_rlWITurn.rsz_LFP{i,1};
        rsz_clk_rlWITurn(:,i)          = SWP_sst_rlWITurn.rsz_clk{i,1};
        rsz_crv_rlWITurn(:,i)          = SWP_sst_rlWITurn.rsz_crv{i,1};
        rsz_wbt_phase_rlWITurn(:,i)    = SWP_sst_rlWITurn.rsz_wbt_phase{i,1};
    end
    
    avg_est_dist_rlWITurn        = mean(est_dist_rlWITurn, 2);         sem_est_dist_rlWITurn        = std(est_dist_rlWITurn, [], 2) / sqrt(N_sweeps_rlWITurn);
    avg_rsz_spk_dsty_rlWITurn    = mean(rsz_spk_dsty_rlWITurn, 2);     sem_rsz_spk_dsty_rlWITurn    = std(rsz_spk_dsty_rlWITurn, [], 2) / sqrt(N_sweeps_rlWITurn);
    avg_rsz_wbt_rlWITurn         = mean(rsz_wbt_rlWITurn, 2);          sem_rsz_wbt_rlWITurn         = std(rsz_wbt_rlWITurn, [], 2) / sqrt(N_sweeps_rlWITurn);
    avg_rsz_LFP_rlWITurn         = mean(rsz_LFP_rlWITurn, 2);          sem_rsz_LFP_rlWITurn         = std(rsz_LFP_rlWITurn, [], 2) / sqrt(N_sweeps_rlWITurn);
    avg_rsz_wbt_phase_rlWITurn   = mean(rsz_wbt_phase_rlWITurn, 2);    sem_rsz_wbt_phase_rlWITurn   = std(rsz_wbt_phase_rlWITurn, [], 2) / sqrt(N_sweeps_rlWITurn);
    avg_rsz_clk_rlWITurn         = mean(rsz_clk_rlWITurn, 2);          sem_rsz_clk_rlWITurn         = std(rsz_clk_rlWITurn, [], 2) / sqrt(N_sweeps_rlWITurn);
    avg_rsz_crv_rlWITurn         = mean(rsz_crv_rlWITurn, 2);          sem_rsz_crv_rlWITurn         = std(rsz_crv_rlWITurn, [], 2) / sqrt(N_sweeps_rlWITurn);
    
    %=== Calculate averages and STDs (REAL DATA, WITHOUT TURN)
    est_dist_rlWOTurn         = zeros(n_reshape, N_sweeps_rlWOTurn);
    rsz_spk_dsty_rlWOTurn     = zeros(n_reshape, N_sweeps_rlWOTurn);
    rsz_wbt_rlWOTurn          = zeros(n_reshape, N_sweeps_rlWOTurn);
    rsz_LFP_rlWOTurn          = zeros(n_reshape, N_sweeps_rlWOTurn);
    rsz_clk_rlWOTurn          = zeros(n_reshape, N_sweeps_rlWOTurn);
    rsz_crv_rlWOTurn          = zeros(n_reshape, N_sweeps_rlWOTurn);
    rsz_wbt_phase_rlWOTurn    = zeros(n_reshape, N_sweeps_rlWOTurn);
    
    for i = 1:N_sweeps_rlWOTurn
        est_dist_rlWOTurn(:,i)         = SWP_sst_rlWOTurn.est_dist1{i,1};
        rsz_spk_dsty_rlWOTurn(:,i)     = SWP_sst_rlWOTurn.rsz_spk_dsty{i,1};
        rsz_wbt_rlWOTurn(:,i)          = SWP_sst_rlWOTurn.rsz_wbt{i,1};
        rsz_LFP_rlWOTurn(:,i)          = SWP_sst_rlWOTurn.rsz_LFP{i,1};
        rsz_clk_rlWOTurn(:,i)          = SWP_sst_rlWOTurn.rsz_clk{i,1};
        rsz_crv_rlWOTurn(:,i)          = SWP_sst_rlWOTurn.rsz_crv{i,1};
        rsz_wbt_phase_rlWOTurn(:,i)    = SWP_sst_rlWOTurn.rsz_wbt_phase{i,1};
    end
    
    avg_est_dist_rlWOTurn        = mean(est_dist_rlWOTurn, 2);         sem_est_dist_rlWOTurn        = std(est_dist_rlWOTurn, [], 2) / sqrt(N_sweeps_rlWOTurn);
    avg_rsz_spk_dsty_rlWOTurn    = mean(rsz_spk_dsty_rlWOTurn, 2);     sem_rsz_spk_dsty_rlWOTurn    = std(rsz_spk_dsty_rlWOTurn, [], 2) / sqrt(N_sweeps_rlWOTurn);
    avg_rsz_wbt_rlWOTurn         = mean(rsz_wbt_rlWOTurn, 2);          sem_rsz_wbt_rlWOTurn         = std(rsz_wbt_rlWOTurn, [], 2) / sqrt(N_sweeps_rlWOTurn);
    avg_rsz_LFP_rlWOTurn         = mean(rsz_LFP_rlWOTurn, 2);          sem_rsz_LFP_rlWOTurn         = std(rsz_LFP_rlWOTurn, [], 2) / sqrt(N_sweeps_rlWOTurn);
    avg_rsz_wbt_phase_rlWOTurn   = mean(rsz_wbt_phase_rlWOTurn, 2);    sem_rsz_wbt_phase_rlWOTurn   = std(rsz_wbt_phase_rlWOTurn, [], 2) / sqrt(N_sweeps_rlWOTurn);
    avg_rsz_clk_rlWOTurn         = mean(rsz_clk_rlWOTurn, 2);          sem_rsz_clk_rlWOTurn         = std(rsz_clk_rlWOTurn, [], 2) / sqrt(N_sweeps_rlWOTurn);
    avg_rsz_crv_rlWOTurn         = mean(rsz_crv_rlWOTurn, 2);          sem_rsz_crv_rlWOTurn         = std(rsz_crv_rlWOTurn, [], 2) / sqrt(N_sweeps_rlWOTurn);
    
    %=== Plot comparison: WITH TURN vs WITHOUT TURN
    figure('units','normalized','outerposition',[.3 .3 .5 .3]);
    tiledlayout(1,7,'TileSpacing','compact');
    nexttile; hold on;
    plotWinterval_AF_v0(rsz_bin_ctrs, avg_est_dist_rlWOTurn, avg_est_dist_rlWOTurn - sem_est_dist_rlWOTurn, avg_est_dist_rlWOTurn + sem_est_dist_rlWOTurn, 'r');
    plotWinterval_AF_v0(rsz_bin_ctrs, avg_est_dist_rlWITurn, avg_est_dist_rlWITurn - sem_est_dist_rlWITurn, avg_est_dist_rlWITurn + sem_est_dist_rlWITurn, 'b');
    xticks([-170 0 170]); title('Decoding Error'); plot([0 0], ylim, 'k--'); xlim([-170 170]); legend('Without','','With')
    
    nexttile; hold on;
    plotWinterval_AF_v0(rsz_bin_ctrs, avg_rsz_spk_dsty_rlWOTurn, avg_rsz_spk_dsty_rlWOTurn - sem_rsz_spk_dsty_rlWOTurn, avg_rsz_spk_dsty_rlWOTurn + sem_rsz_spk_dsty_rlWOTurn, 'r');
    plotWinterval_AF_v0(rsz_bin_ctrs, avg_rsz_spk_dsty_rlWITurn, avg_rsz_spk_dsty_rlWITurn - sem_rsz_spk_dsty_rlWITurn, avg_rsz_spk_dsty_rlWITurn + sem_rsz_spk_dsty_rlWITurn, 'b');
    xticks([-170 0 170]); title('Spike Density'); plot([0 0], ylim, 'k--'); xlim([-170 170]);
    
    nexttile; hold on;
    plotWinterval_AF_v0(rsz_bin_ctrs, avg_rsz_LFP_rlWOTurn, avg_rsz_LFP_rlWOTurn - sem_rsz_LFP_rlWOTurn, avg_rsz_LFP_rlWOTurn + sem_rsz_LFP_rlWOTurn, 'r');
    plotWinterval_AF_v0(rsz_bin_ctrs, avg_rsz_LFP_rlWITurn, avg_rsz_LFP_rlWITurn - sem_rsz_LFP_rlWITurn, avg_rsz_LFP_rlWITurn + sem_rsz_LFP_rlWITurn, 'b');
    xticks([-170 0 170]); title('LFP'); plot([0 0], ylim, 'k--'); xlim([-170 170]);
    
    nexttile; hold on;
    plotWinterval_AF_v0(rsz_bin_ctrs, avg_rsz_clk_rlWOTurn, avg_rsz_clk_rlWOTurn - sem_rsz_clk_rlWOTurn, avg_rsz_clk_rlWOTurn + sem_rsz_clk_rlWOTurn, 'r');
    plotWinterval_AF_v0(rsz_bin_ctrs, avg_rsz_clk_rlWITurn, avg_rsz_clk_rlWITurn - sem_rsz_clk_rlWITurn, avg_rsz_clk_rlWITurn + sem_rsz_clk_rlWITurn, 'b');
    xticks([-170 0 170]); title('Call Rate'); plot([0 0], ylim, 'k--'); xlim([-170 170]);
    
    nexttile; hold on;
    plotWinterval_AF_v0(rsz_bin_ctrs, avg_rsz_crv_rlWOTurn, avg_rsz_crv_rlWOTurn - sem_rsz_crv_rlWOTurn, avg_rsz_crv_rlWOTurn + sem_rsz_crv_rlWOTurn, 'r');
    plotWinterval_AF_v0(rsz_bin_ctrs, avg_rsz_crv_rlWITurn, avg_rsz_crv_rlWITurn - sem_rsz_crv_rlWITurn, avg_rsz_crv_rlWITurn + sem_rsz_crv_rlWITurn, 'b');
    xticks([-170 0 170]); title('Curvature'); plot([0 0], ylim, 'k--'); xlim([-170 170]);
    
    nexttile; hold on;
    plotWinterval_AF_v0(rsz_bin_ctrs, avg_rsz_wbt_phase_rlWOTurn, avg_rsz_wbt_phase_rlWOTurn - sem_rsz_wbt_phase_rlWOTurn, avg_rsz_wbt_phase_rlWOTurn + sem_rsz_wbt_phase_rlWOTurn, 'r');
    plotWinterval_AF_v0(rsz_bin_ctrs, avg_rsz_wbt_phase_rlWITurn, avg_rsz_wbt_phase_rlWITurn - sem_rsz_wbt_phase_rlWITurn, avg_rsz_wbt_phase_rlWITurn + sem_rsz_wbt_phase_rlWITurn, 'b');
    xticks([-170 0 170]); title('Wingbeat Phase'); plot([0 0], ylim, 'k--'); xlim([-170 170]);
    
    nexttile; hold on;
    plotWinterval_AF_v0(rsz_bin_ctrs, avg_rsz_wbt_rlWOTurn, avg_rsz_wbt_rlWOTurn - sem_rsz_wbt_rlWOTurn, avg_rsz_wbt_rlWOTurn + sem_rsz_wbt_rlWOTurn, 'r');
    plotWinterval_AF_v0(rsz_bin_ctrs, avg_rsz_wbt_rlWITurn, avg_rsz_wbt_rlWITurn - sem_rsz_wbt_rlWITurn, avg_rsz_wbt_rlWITurn + sem_rsz_wbt_rlWITurn, 'b');
    xticks([-170 0 170]); title('Accelerometer'); plot([0 0], ylim, 'k--'); xlim([-170 170]);
    
    sgtitle(['Avg of ', num2str(N_sweeps_rlWITurn), ' (W/ turns) vs ', num2str(N_sweeps_rlWOTurn), ' (W/O turns) sweeps']);
    
    %=== Plot decoding error on top of call rate
    figure('units','normalized','outerposition',[.3 .1 .15 .3]);
    hold on;
    plot(rsz_bin_ctrs, normalize(avg_est_dist_rlWITurn,'range'));
    plot(rsz_bin_ctrs, normalize(avg_rsz_clk_rlWITurn,'range'));
    xticks([-170 0 170]); title('Decoding error vs echolocation'); plot([0 0], ylim, 'k--'); xlim([-170 170]);
    
    
end

%% LOOK AT AUTOMATICALLY CLASSIFIED SWEEPS FOR ECHOLOCATION AND TURN ANALYSIS
for hide = 1
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
    SWPs_acc = {}; SWPs_sample = []; wbt_phase_all = []; tmr_phase_all = []; tmr_power_all = []; fract_pos_all = [];    SWPs_phs = [];  PSD_sweeps = [];    WFR_sweeps = [];    SWPs_clk = [];  SWPs_trn = [];
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
                SWPs_trn = [SWPs_trn; {SWPs_sst.curv_flag{zz,1}}];
                
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
    
    %%
    %=== Extract good sweeps
    S2W = table();                                          % Sweep to wingbeat table
    S2W.flightID = SWPs_flightID(good_ids);                 % Id of the flights
    S2W.acc = SWPs_acc(good_ids);                           % Accelerometer
    S2W.phs = SWPs_phs(good_ids);                           % Wingbeat Phase
    S2W.trn = SWPs_trn(good_ids);                           % Turn flag
    S2W.trn = cellfun(@(x) find(x>0),S2W.trn,'UniformOutput',false);                           % Turn flag
    S2W.sample = SWPs_sample(good_ids);                     % Sample of the sweep
    S2W.u_id = SWPs_u_id(good_ids,:);                       % Unique identifier
    [~,~,S2W.sessionID] = unique([string(S2W.u_id)],'rows');% Identifier
    
    %=== Extract subset happening in flights with turns
    S2W_sst = S2W(~cellfun(@isempty,S2W.trn),:);
    S2W_sst.dur = cellfun(@(x) numel(x)./Fs_Hi,S2W_sst.acc);
    S2W_sst.trn = cell2mat(S2W_sst.trn);
    
    S2W_sst.turn2end = S2W_sst.dur-S2W_sst.trn./Fs_Hi;
    S2W_sst.turn2str = S2W_sst.trn./Fs_Hi;
    S2W_sst.sweep2trn = (S2W_sst.sample-S2W_sst.trn)./Fs_Hi;
    
    %=== Show distribution
    data = S2W_sst.sweep2trn(S2W_sst.sweep2trn>-1.5 & S2W_sst.sweep2trn<1.5);
    [x_vals, density] = ksdensity(data,'Bandwidth', 0.3);
    
    % Create a histogram (optional, just for comparison)
    figure('units','normalized','outerposition',[.3 .1 .1 .3]);
    histogram(data,[-1.5:0.1:1.5], 'Normalization', 'pdf');    hold on;
    plot(density, x_vals, 'LineWidth', 2);
    xlabel('Time from turn (s)');   ylabel('Sweep Probability');    xlim([-1.5 1.5]);
    plot([0 0],ylim,'k--');
    
    unique(S2W_sst.flightID)
    unique([string(S2W_sst.u_id)],'rows')
    
    %% ==================================================== ANALYSIS OF SWEEP VS ECHOLOCATION =========================================================================
    
    id2keep = good_ids;
    S2W_st = table();                                               % Sweep to wingbeat table
    S2W_st.flightID = SWPs_flightID(id2keep);                       % Id of the flights
    S2W_st.acc = SWPs_acc(id2keep);                                 % Accelerometer
    S2W_st.phs = SWPs_phs(id2keep);                                 % Wingbeat Phase
    S2W_st.clk = SWPs_clk(id2keep);                                 % Echolocation
    S2W_st.tmr = SWPs_tmr(id2keep);                                 % Tamir Phase
    S2W_st.sample = SWPs_sample(id2keep);                           % Sample of the sweep
    S2W_st.u_id = SWPs_u_id(id2keep,:);                             % Unique identifier
    [~,~,S2W_st.sessionID] = unique([string(S2W_st.u_id)],'rows');  % Identifier
    
    %=== Keep only sweeps that start away from flight beginning/end and that come from sessions with mic
    S2W_st = S2W_st(S2W_st.sample>30 & S2W_st.sample<370 & ~any(contains(string(S2W_st.u_id), "14445") | contains(string(S2W_st.u_id), "14611"),2),:);
    unique([string(S2W_st.u_id)],'rows');
    n_rep1 = 50;
    dt_interval = 0.04; % 0.04 default
    interval = [-round(Fs_Hi*dt_interval):round(Fs_Hi*dt_interval)];
    
    acc_accum_gd = zeros(size(S2W_st,1),numel(interval));
    acc_accum_sh_gd = zeros(size(S2W_st,1),numel(interval),n_rep1);
    clk_accum_gd = zeros(size(S2W_st,1),numel(interval));
    clk_accum_sh_gd = zeros(size(S2W_st,1),numel(interval),n_rep1);
    tmr_accum_gd = zeros(size(S2W_st,1),numel(interval));
    tmr_accum_sh_gd = zeros(size(S2W_st,1),numel(interval),n_rep1);
    for i=1:size(S2W_st,1)
        
        acc_accum_gd(i,:) = S2W_st.acc{i,1}(S2W_st.sample(i)+interval,1)';
        clk_accum_gd(i,:) = S2W_st.clk{i,1}(S2W_st.sample(i)+interval,1)';
        tmr_accum_gd(i,:) = S2W_st.tmr{i,1}(S2W_st.sample(i)+interval,1)';
        for jj =1:n_rep1
            sample_shuffle = 50+randi(300);
            acc_accum_sh_gd(i,:,jj) = S2W_st.acc{i,1}(sample_shuffle+interval,1)';
            clk_accum_sh_gd(i,:,jj) = S2W_st.clk{i,1}(sample_shuffle+interval,1)';
            tmr_accum_sh_gd(i,:,jj) = S2W_st.tmr{i,1}(sample_shuffle+interval,1)';
        end
        
    end
    acc_accum_sh_mean_gd = squeeze(mean(acc_accum_sh_gd,3));
    clk_accum_sh_mean_gd = squeeze(mean(clk_accum_sh_gd,3));
    tmr_accum_sh_mean_gd = squeeze(mean(tmr_accum_sh_gd,3));
    
    figure('units','normalized','outerposition',[.1 .1 .3 .3]);
    tiledlayout(1,3,'TileSpacing','tight');
    nexttile;
    plotWinterval_AF_v0(interval./Fs_Hi,mean(acc_accum_gd),mean(acc_accum_gd)-std(acc_accum_gd)./sqrt(size(S2W_st,1)),mean(acc_accum_gd)+std(acc_accum_gd)./sqrt(size(S2W_st,1)),'r');    hold on;
    plotWinterval_AF_v0(interval./Fs_Hi,mean(acc_accum_sh_mean_gd),mean(acc_accum_sh_mean_gd)-std(acc_accum_sh_mean_gd)./sqrt(size(S2W_st,1)),mean(acc_accum_sh_mean_gd)+std(acc_accum_sh_mean_gd)./sqrt(size(S2W_st,1)),'k');
    xlabel('Time from Sweep Peak (s)'); ylabel('Average Acceleration (g)');         ylim([-0.05 0.05]);
    nexttile;
    plotWinterval_AF_v0(interval./Fs_Hi,mean(clk_accum_gd),mean(clk_accum_gd)-std(clk_accum_gd)./sqrt(size(S2W_st,1)),mean(clk_accum_gd)+std(clk_accum_gd)./sqrt(size(S2W_st,1)),'r');    hold on;
    plotWinterval_AF_v0(interval./Fs_Hi,mean(clk_accum_sh_mean_gd),mean(clk_accum_sh_mean_gd)-std(clk_accum_sh_mean_gd)./sqrt(size(S2W_st,1)),mean(clk_accum_sh_mean_gd)+std(clk_accum_sh_mean_gd)./sqrt(size(S2W_st,1)),'k');
    xlabel('Time from Sweep Peak (s)'); ylabel('Average Echolocation Rate (Hz)');   ylim([0.05 0.13]);
    nexttile;
    plotWinterval_AF_v0(interval./Fs_Hi,mean(tmr_accum_gd),mean(tmr_accum_gd)-std(tmr_accum_gd)./sqrt(size(S2W_st,1)),mean(tmr_accum_gd)+std(tmr_accum_gd)./sqrt(size(S2W_st,1)),'r');    hold on;
    plotWinterval_AF_v0(interval./Fs_Hi,mean(tmr_accum_sh_mean_gd),mean(tmr_accum_sh_mean_gd)-std(tmr_accum_sh_mean_gd)./sqrt(size(S2W_st,1)),mean(tmr_accum_sh_mean_gd)+std(tmr_accum_sh_mean_gd)./sqrt(size(S2W_st,1)),'k');
    xlabel('Time from Sweep Peak (s)'); ylabel('Average LFP (uV)');                 ylim([-10 3]);
    sgtitle('Good Sweeps');
    
    id2keep = baad_ids;
    S2W_st = table();                                         % Sweep to wingbeat table
    S2W_st.flightID = SWPs_flightID(id2keep);                 % Id of the flights
    S2W_st.acc = SWPs_acc(id2keep);                           % Accelerometer
    S2W_st.phs = SWPs_phs(id2keep);                           % Wingbeat Phase
    S2W_st.clk = SWPs_clk(id2keep);                           % Echolocation
    S2W_st.tmr = SWPs_tmr(id2keep);                           % Tamir Phase
    S2W_st.sample = SWPs_sample(id2keep);                     % Sample of the sweep
    S2W_st.u_id = SWPs_u_id(id2keep,:);                       % Unique identifier
    [~,~,S2W_st.sessionID] = unique([string(S2W_st.u_id)],'rows');% Identifier
    
    %=== Keep only sweeps that start away from flight beginning/end and that come from sessions with mic
    S2W_st = S2W_st(S2W_st.sample>30 & S2W_st.sample<370 & ~any(contains(string(S2W_st.u_id), "14445") | contains(string(S2W_st.u_id), "14611"),2),:);
    unique([string(S2W_st.u_id)],'rows');
    
    acc_accum_bd = zeros(size(S2W_st,1),numel(interval));
    acc_accum_sh_bd = zeros(size(S2W_st,1),numel(interval),n_rep1);
    clk_accum_bd = zeros(size(S2W_st,1),numel(interval));
    clk_accum_sh_bd = zeros(size(S2W_st,1),numel(interval),n_rep1);
    tmr_accum_bd = zeros(size(S2W_st,1),numel(interval));
    tmr_accum_sh_bd = zeros(size(S2W_st,1),numel(interval),n_rep1);
    for i=1:size(S2W_st,1)
        
        acc_accum_bd(i,:) = S2W_st.acc{i,1}(S2W_st.sample(i)+interval,1)';
        clk_accum_bd(i,:) = S2W_st.clk{i,1}(S2W_st.sample(i)+interval,1)';
        tmr_accum_bd(i,:) = S2W_st.tmr{i,1}(S2W_st.sample(i)+interval,1)';
        for jj =1:n_rep1
            sample_shuffle = 50+randi(300);
            acc_accum_sh_bd(i,:,jj) = S2W_st.acc{i,1}(sample_shuffle+interval,1)';
            clk_accum_sh_bd(i,:,jj) = S2W_st.clk{i,1}(sample_shuffle+interval,1)';
            tmr_accum_sh_bd(i,:,jj) = S2W_st.tmr{i,1}(sample_shuffle+interval,1)';
        end
        
    end
    acc_accum_sh_mean_bd = squeeze(mean(acc_accum_sh_bd,3));
    clk_accum_sh_mean_bd = squeeze(mean(clk_accum_sh_bd,3));
    tmr_accum_sh_mean_bd = squeeze(mean(tmr_accum_sh_bd,3));
    
    figure('units','normalized','outerposition',[.4 .1 .3 .3]);
    tiledlayout(1,3,'TileSpacing','tight');
    nexttile;
    plotWinterval_AF_v0(interval./Fs_Hi,mean(acc_accum_bd),mean(acc_accum_bd)-std(acc_accum_bd)./sqrt(size(S2W_st,1)),mean(acc_accum_bd)+std(acc_accum_bd)./sqrt(size(S2W_st,1)),'b');    hold on;
    plotWinterval_AF_v0(interval./Fs_Hi,mean(acc_accum_sh_mean_bd),mean(acc_accum_sh_mean_bd)-std(acc_accum_sh_mean_bd)./sqrt(size(S2W_st,1)),mean(acc_accum_sh_mean_bd)+std(acc_accum_sh_mean_bd)./sqrt(size(S2W_st,1)),'k');
    xlabel('Time from Sweep Peak (s)'); ylabel('Average Acceleration (g)');         ylim([-0.05 0.05]);
    nexttile;
    plotWinterval_AF_v0(interval./Fs_Hi,mean(clk_accum_bd),mean(clk_accum_bd)-std(clk_accum_bd)./sqrt(size(S2W_st,1)),mean(clk_accum_bd)+std(clk_accum_bd)./sqrt(size(S2W_st,1)),'b');    hold on;
    plotWinterval_AF_v0(interval./Fs_Hi,mean(clk_accum_sh_mean_bd),mean(clk_accum_sh_mean_bd)-std(clk_accum_sh_mean_bd)./sqrt(size(S2W_st,1)),mean(clk_accum_sh_mean_bd)+std(clk_accum_sh_mean_bd)./sqrt(size(S2W_st,1)),'k');
    xlabel('Time from Sweep Peak (s)'); ylabel('Average Echolocation Rate (Hz)');   ylim([0.05 0.13]);
    nexttile;
    plotWinterval_AF_v0(interval./Fs_Hi,mean(tmr_accum_bd),mean(tmr_accum_bd)-std(tmr_accum_bd)./sqrt(size(S2W_st,1)),mean(tmr_accum_bd)+std(tmr_accum_bd)./sqrt(size(S2W_st,1)),'b');    hold on;
    plotWinterval_AF_v0(interval./Fs_Hi,mean(tmr_accum_sh_mean_bd),mean(tmr_accum_sh_mean_bd)-std(tmr_accum_sh_mean_bd)./sqrt(size(S2W_st,1)),mean(tmr_accum_sh_mean_bd)+std(tmr_accum_sh_mean_bd)./sqrt(size(S2W_st,1)),'k');
    xlabel('Time from Sweep Peak (s)'); ylabel('Average LFP (uV)');                 ylim([-10 3]);
    sgtitle('Bad Sweeps');
    
    %=== Good vs bad
    figure('units','normalized','outerposition',[.4 .1 .3 .3]);
    tiledlayout(1,3,'TileSpacing','tight');
    nexttile;
    plotWinterval_AF_v0(interval./Fs_Hi,mean(acc_accum_bd),mean(acc_accum_bd)-std(acc_accum_bd)./sqrt(size(S2W_st,1)),mean(acc_accum_bd)+std(acc_accum_bd)./sqrt(size(S2W_st,1)),'b');    hold on;
    plotWinterval_AF_v0(interval./Fs_Hi,mean(acc_accum_gd),mean(acc_accum_gd)-std(acc_accum_gd)./sqrt(size(S2W_st,1)),mean(acc_accum_gd)+std(acc_accum_gd)./sqrt(size(S2W_st,1)),'r');
    xlabel('Time from Sweep Peak (s)'); ylabel('Average Acceleration (g)');
    nexttile;
    plotWinterval_AF_v0(interval./Fs_Hi,mean(clk_accum_bd),mean(clk_accum_bd)-std(clk_accum_bd)./sqrt(size(S2W_st,1)),mean(clk_accum_bd)+std(clk_accum_bd)./sqrt(size(S2W_st,1)),'b');    hold on;
    plotWinterval_AF_v0(interval./Fs_Hi,mean(clk_accum_gd),mean(clk_accum_gd)-std(clk_accum_gd)./sqrt(size(S2W_st,1)),mean(clk_accum_gd)+std(clk_accum_gd)./sqrt(size(S2W_st,1)),'r');
    xlabel('Time from Sweep Peak (s)'); ylabel('Average Echolocation Rate (Hz)');
    nexttile;
    plotWinterval_AF_v0(interval./Fs_Hi,mean(tmr_accum_bd),mean(tmr_accum_bd)-std(tmr_accum_bd)./sqrt(size(S2W_st,1)),mean(tmr_accum_bd)+std(tmr_accum_bd)./sqrt(size(S2W_st,1)),'b');    hold on;
    plotWinterval_AF_v0(interval./Fs_Hi,mean(tmr_accum_gd),mean(tmr_accum_gd)-std(tmr_accum_gd)./sqrt(size(S2W_st,1)),mean(tmr_accum_gd)+std(tmr_accum_gd)./sqrt(size(S2W_st,1)),'r');
    xlabel('Time from Sweep Peak (s)'); ylabel('Average LFP (uV)');
    sgtitle('Good vs Bad Sweeps');
    
    %=== Quantify and compare echolocation rate around sweeps
    click_rate_gd = mean(clk_accum_gd,2);
    click_rate_bd = mean(clk_accum_bd,2);
    click_rate_sh = mean(clk_accum_sh_gd,[2,3]);
    
    % Data
    means = [mean(click_rate_gd), mean(click_rate_sh)];                        % mean values
    sems  = [std(click_rate_gd)/sqrt(numel(click_rate_sh)), std(click_rate_sh)/sqrt(numel(click_rate_sh))];                        % standard error of the mean
    
    mean_and_sem_AF_v0(click_rate_gd)
    mean_and_sem_AF_v0(click_rate_sh)
    
    % Bar plot
    figure; b = bar(means, 'FaceColor', 'flat');
    hold on;    x = 1:length(means);    errorbar(x, means, sems, 'k', 'LineStyle', 'none', 'LineWidth', 1.5);
    title(signrank(click_rate_gd,click_rate_sh));    xticklabels({'Observed','Random'});
    
    %=== Plot using confidence intervals
    figure('units','normalized','outerposition',[.4 .1 .2 .3]);
    tiledlayout(1,2,'TileSpacing','tight');
    nexttile;
    ci = bootci(100, @mean, acc_accum_bd); plotWinterval_AF_v0(interval./Fs_Hi,mean(acc_accum_bd),ci(2,:), ci(1,:),'b');    hold on;
    ci = bootci(100, @mean, acc_accum_gd); plotWinterval_AF_v0(interval./Fs_Hi,mean(acc_accum_gd),ci(2,:), ci(1,:),'r');
    xlabel('Time from Sweep Peak (s)'); ylabel('Average Acceleration (g)');
    nexttile;
    ci = bootci(100, @mean, tmr_accum_bd); plotWinterval_AF_v0(interval./Fs_Hi,mean(tmr_accum_bd),ci(2,:), ci(1,:),'b');    hold on;
    ci = bootci(100, @mean, tmr_accum_gd); plotWinterval_AF_v0(interval./Fs_Hi,mean(tmr_accum_gd),ci(2,:), ci(1,:),'r');
    xlabel('Time from Sweep Peak (s)'); ylabel('Average LFP (uV)');
    sgtitle('Good vs Bad Sweeps (bootstrapping) CIs');
    
    
end

