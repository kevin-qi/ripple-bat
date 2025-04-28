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

%% LOOK AT LFP
for hide=1
    
    %=== Define frequency range
    freq_smp = 0.005;
    freq = [0:freq_smp:250];
    
    %=== Initialize relevant vaiables
    PSD_r = zeros(size(LFP,1),numel(freq));
    PSD_f = zeros(size(LFP,1),numel(freq));
    PSD_b = zeros(size(LFP,1),numel(freq));
    PSD_a = zeros(size(LFP,1),numel(freq));
    PSD_pre = zeros(size(LFP,1),numel(freq));
    PSD_dur = zeros(size(LFP,1),numel(freq));
    PSD_pst = zeros(size(LFP,1),numel(freq));
    
    %=== Reinterpolate the LFP on the same frequency range
    N_s = size(LFP,1);
    for i=1:size(LFP,1)
        PSD_r(i,:)  = interp1(LFP(i).f_PSD_r,LFP(i).PSD_r,freq,'linear','extrap');
        PSD_f(i,:)  = interp1(LFP(i).f_PSD_f,LFP(i).PSD_f,freq,'linear','extrap');
        PSD_b(i,:)  = interp1(LFP(i).f_PSD_b,LFP(i).PSD_b,freq,'linear','extrap');
        PSD_a(i,:)  = interp1(LFP(i).f_PSD_a,LFP(i).PSD_a,freq,'linear','extrap');
        PSD_pre(i,:)  = interp1(LFP(i).f_PSD_pre,LFP(i).PSD_pre,freq,'linear','extrap');
        PSD_dur(i,:)  = interp1(LFP(i).f_PSD_dur,LFP(i).PSD_dur,freq,'linear','extrap');
        PSD_pst(i,:)  = interp1(LFP(i).f_PSD_pst,LFP(i).PSD_pst,freq,'linear','extrap');
    end
    
    %=== Normalize on the whole spectrum and scale to have unitary area
    PSD_r = PSD_r./sum(PSD_r,2);
    PSD_f = PSD_f./sum(PSD_f,2);
    PSD_b = PSD_b./sum(PSD_b,2);
    PSD_a = PSD_a./sum(PSD_a,2);
    PSD_pre = PSD_pre./sum(PSD_pre,2);
    PSD_dur = PSD_dur./sum(PSD_dur,2);
    PSD_pst = PSD_pst./sum(PSD_pst,2);
    
    %=== Average and smooth
    avg_PSD_r =   smoothdata(mean(PSD_r,1),  'movmedian',1/freq_smp);
    avg_PSD_f =   smoothdata(mean(PSD_f,1),  'movmedian',1/freq_smp);
    avg_PSD_b =   smoothdata(mean(PSD_b,1),  'movmedian',1/freq_smp);
    avg_PSD_a =   smoothdata(mean(PSD_a,1),  'movmedian',1/freq_smp);
    avg_PSD_pre = smoothdata(mean(PSD_pre,1),'movmedian',1/freq_smp);
    avg_PSD_dur = smoothdata(mean(PSD_dur,1),'movmedian',1/freq_smp);
    avg_PSD_pst = smoothdata(mean(PSD_pst,1),'movmedian',1/freq_smp);
    
    std_PSD_r =   smoothdata(std(PSD_r,[],1)./sqrt(N_s),  'movmedian',1/freq_smp);
    std_PSD_f =   smoothdata(std(PSD_f,[],1)./sqrt(N_s),  'movmedian',1/freq_smp);
    std_PSD_b =   smoothdata(std(PSD_b,[],1)./sqrt(N_s),  'movmedian',1/freq_smp);
    std_PSD_a =   smoothdata(std(PSD_a,[],1)./sqrt(N_s),  'movmedian',1/freq_smp);
    std_PSD_pre = smoothdata(std(PSD_pre,[],1)./sqrt(N_s),'movmedian',1/freq_smp);
    std_PSD_dur = smoothdata(std(PSD_dur,[],1)./sqrt(N_s),'movmedian',1/freq_smp);
    std_PSD_pst = smoothdata(std(PSD_pst,[],1)./sqrt(N_s),'movmedian',1/freq_smp);
    
    %=== Plot the PSDs
    figure('units','normalized','outerposition',[0.1 0.3 0.3 0.3]);
    tiledlayout(1,3,'TileSpacing','tight');
    nexttile;   plot(freq,avg_PSD_r,'k','LineWidth',3);    hold on;     plot(freq,avg_PSD_f,'r','LineWidth',3);
    xlim([0 20]);   xlabel('Frequency (Hz)');   ylabel('PSD (Norm.)');  legend('Rest','Flight');        title('LFP');
    nexttile;   plot(freq,avg_PSD_r,'k','LineWidth',3);    hold on;     plot(freq,avg_PSD_b,'g','LineWidth',3);
    xlim([0 20]);   xlabel('Frequency (Hz)');   ylabel('PSD (Norm.)');  legend('Rest','0-Bouts');        title('LFP');
    nexttile;   plot(freq,avg_PSD_a,'k','LineWidth',3);
    xlim([0 20]);   xlabel('Frequency (Hz)');   ylabel('PSD (Norm.)');  title('Accelerometer');
    
    %============================================================================================================================================================================================================
    %=== Some processing on the plotting of the accelerometer spectra
    [max_m,max_loc] = max(PSD_a(:,freq>1 & freq<11),[],2);
    %PSD_a_norm = smoothdata(PSD_a./max_m,2,'gaussian',100);
    PSD_a_norm = PSD_a./max_m;
    [~,sorted_idx] = sort(max_loc);
    col = spring(numel(sorted_idx));

    figure('units','normalized','outerposition',[0.1 0.3 0.3 0.4]);
    tiledlayout(1,2,'TileSpacing','tight');
    nexttile;   
    for i=1:numel(sorted_idx)
    plot(freq,PSD_a_norm(sorted_idx(i),:),'LineWidth',1,'Color',col(i,:));  hold on;
    end
    plot(freq,mean(PSD_a_norm,1),'LineWidth',4,'Color','k');
    xlim([6 10]);
    nexttile;   data_tmp = PSD_a_norm;
    plotWinterval_AF_v0(freq,mean(data_tmp,'omitnan'),mean(data_tmp,'omitnan')-std(data_tmp,[],'omitnan')./sqrt(size(data_tmp,1)),mean(data_tmp,'omitnan')+std(data_tmp,[],'omitnan')./sqrt(size(data_tmp,1)),'r');
    xlim([6 10]);
    %============================================================================================================================================================================================================
    
    %=== Plot the PSDs
    figure('units','normalized','outerposition',[0.1 0.6 0.3 0.3]);
    tiledlayout(1,3,'TileSpacing','tight');
    nexttile;
    plotWinterval_AF_v0(freq,avg_PSD_r,avg_PSD_r-std_PSD_r,avg_PSD_r+std_PSD_r,'k');    hold on;
    plotWinterval_AF_v0(freq,avg_PSD_f,avg_PSD_f-std_PSD_f,avg_PSD_f+std_PSD_f,'r');    hold on;
    xlim([0 20]);   xlabel('Frequency (Hz)');   ylabel('PSD (Norm.)');  title('LFP');
    nexttile;
    plotWinterval_AF_v0(freq,avg_PSD_r,avg_PSD_r-std_PSD_r,avg_PSD_r+std_PSD_r,'k');    hold on;
    plotWinterval_AF_v0(freq,avg_PSD_b,avg_PSD_b-std_PSD_b,avg_PSD_b+std_PSD_b,'g');    hold on;
    xlim([0 20]);   xlabel('Frequency (Hz)');   ylabel('PSD (Norm.)');  title('LFP');
    nexttile;
    plotWinterval_AF_v0(freq,avg_PSD_a,avg_PSD_a-std_PSD_a,avg_PSD_a+std_PSD_a,'k');    hold on;
    xlim([0 20]);   xlabel('Frequency (Hz)');   ylabel('PSD (Norm.)');  title('Accelerometer');
    
    %=== Plot the PSDs
    figure('units','normalized','outerposition',[0.4 0.3 0.1 0.3]);
    plot(freq,avg_PSD_pre,'k','LineWidth',3);    hold on;     plot(freq,avg_PSD_dur,'r','LineWidth',3);   plot(freq,avg_PSD_pst,'b','LineWidth',3);
    xlim([0 20]);   xlabel('Frequency (Hz)');   ylabel('PSD (Norm.)');  legend('LFP (Pre)','LFP (Dur)','LFP (Pst)');
    
    %=== Find peak of wingbeat signal and plot the PSDs
    [~,max_loc] = max(avg_PSD_a);
    wb_freq = freq(max_loc);
    figure('units','normalized','outerposition',[0.4 0.3 0.15 0.3]);
    plot(freq,avg_PSD_a,'k','LineWidth',3);
    xlim([0 20]);   xlabel('Frequency (Hz)');   ylabel('PSD (Norm.)');  title('Accelerometer');
    
    %=== Quantify a few features of the LFP
    LFP_table = struct2table(LFP);
    tt2dt = vertcat(LFP.theta2delta_flgt);
    dur_pre_diff = vertcat(LFP.delta_tht_DurPre);
    fract_flights_wTht = mean(cellfun(@(x,y,z) sum(x>0 & y<0.05 & z>2),LFP_table.delta_tht_DurPre,LFP_table.p_val_tht_DurPre,LFP_table.theta2delta_flgt)./LFP_table.f_num);
    pwr_b = LFP_table.perc_thtpower_b;
    pwr_f = LFP_table.perc_thtpower_f;
    disp([num2str(fract_flights_wTht*sum(LFP_table.f_num),3),' out of ', num2str(sum(LFP_table.f_num),4), ' flights with theta (', num2str(fract_flights_wTht,3),')' ]);
    
    %=== Plot a few features of the single flights
    figure('units','normalized','outerposition',[0.1 0.2 0.3 0.3]);
    tiledlayout(1,3,'TileSpacing','tight');
    nexttile;  histogram(tt2dt,[0:0.1:5],'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','k');
    hold on;    plot([2 2],ylim,'k');   title(['Fraction flights: ', num2str(sum(tt2dt>2)/numel(tt2dt),2)]);  ylabel('Fraction'); xlabel('Theta to Delta Ratio');
    nexttile;  histogram(dur_pre_diff,[-2:0.1:2],'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','k');
    hold on;    plot([0 0],ylim,'k');   title(['Fraction flight: ', num2str(sum(dur_pre_diff>0)/numel(dur_pre_diff),2)]);  ylabel('Fraction'); xlabel('Delta Theta (Flight-Pre)');
    nexttile;
    plot_distr_AF_v0(pwr_b, pwr_f, {'Bouts', 'Flight'}, 'SEM', 'Relative Theta Power'); ylim([0 max([pwr_b;pwr_f])]);
    
    histogram(vertcat(LFP_table.tht_bouts_dur{:}),'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','k');
    
    signrank(pwr_b, pwr_f)
    disp(['Fraction bouts during rest: ', num2str(1-mean([LFP.fract_bouts_flight]),3)]);
    std(1-[LFP.fract_bouts_flight])./sqrt(numel([LFP.fract_bouts_flight]))
    
    %=== Look at Tamir's power
    figure('units','normalized','outerposition',[0.1 0.2 0.2 0.4]);
    tiledlayout(1,2,'TileSpacing','tight');
    nexttile; plot_distr_AF_v0(LFP_table.tmr_power1(:,1), LFP_table.tmr_power1(:,2), {'Flight', 'Rest'}, 'SEM', 'Mean Non-Oscillatory Power');  hold on; 
    title(['p = ',num2str(signrank(LFP_table.tmr_power(:,1), LFP_table.tmr_power(:,2)),3)]);
    nexttile; plot_distr_AF_v0(LFP_table.tmr_power1(:,1), LFP_table.tmr_power1(:,3), {'Flight', 'Rest (Match)'}, 'SEM', 'Mean Non-Oscillatory Power');  hold on; 
    title(['p = ',num2str(signrank(LFP_table.tmr_power(:,1), LFP_table.tmr_power(:,3)),3)]);
    
    figure('units','normalized','outerposition',[0.1 0.2 0.1 0.4]);
    boxchart([LFP_table.tmr_power1(:,1), LFP_table.tmr_power1(:,3)]);
    title(['p = ',num2str(signrank(LFP_table.tmr_power(:,1), LFP_table.tmr_power(:,3)),3)]);
    ylabel('Power');    xticklabels({'Flight','Rest'});    ylim_tmp = ylim;    ylim([0 ylim_tmp(2)]);
    
    %=== Look at spike-tgd LFP
    figure('units','normalized','outerposition',[0.3 0.2 0.15 0.3]);
    hold on;
    data = vertcat(LFP_table.spkTgdLfP_rest{:});
    plotWinterval_AF_v0(LFP_table.spkTgdLfP_time{1,1},mean(data),mean(data)-std(data)./sqrt(size(data,1)),mean(data)+std(data)./sqrt(size(data,1)),'k');
    data = vertcat(LFP_table.spkTgdLfP_flight{:});
    plotWinterval_AF_v0(LFP_table.spkTgdLfP_time{1,1},mean(data),mean(data)-std(data)./sqrt(size(data,1)),mean(data)+std(data)./sqrt(size(data,1)),'r');
    xlim('tight');  xlabel('Time from spike (s)');  ylabel('LFP (uV)'); title('Spike Triggered LFP (All sessions)');    legend('Rest','','Flight');
    
    %=== Look at spike-tgd LFP and Wingbeat from all sessions
    figure('units','normalized','outerposition',[0.3 0.2 0.45 0.3]);
    tiledlayout(1,2,'TileSpacing','tight');
    nexttile;   data = vertcat(LFP_table.spkTgdLfP_rest{:});    data = highpass(data',2,500)';
    plotWinterval_AF_v0(LFP_table.spkTgdLfP_time{1,1},mean(data),mean(data)-std(data)./sqrt(size(data,1)),mean(data)+std(data)./sqrt(size(data,1)),'k');    
    hold on;
    data = vertcat(LFP_table.spkTgdLfP_flight{:});  data = highpass(data',2,500)';
    plotWinterval_AF_v0(LFP_table.spkTgdLfP_time{1,1},mean(data),mean(data)-std(data)./sqrt(size(data,1)),mean(data)+std(data)./sqrt(size(data,1)),'r');    
    plot([0 0],ylim,'k--'); xlim('tight');  xlabel('Time from spike (s)');  ylabel('LFP (uV)'); xlim([-2 2]); 
    
    nexttile;   data = vertcat(LFP_table.spkTgdWbt_flight{:});  data = highpass(data',2,500)';
    plotWinterval_AF_v0(LFP_table.spkTgdLfP_time{1,1},mean(data),mean(data)-std(data)./sqrt(size(data,1)),mean(data)+std(data)./sqrt(size(data,1)),'r');    
    plot([0 0],ylim,'k--'); xlim('tight');  xlabel('Time from spike (s)');  ylabel('Accelerometer (g)'); xlim([-2 2]);
    
end

%% ADD SOME FEATURES TO EACH UNIT
for hide=1
    
    N_cells = size(NP_units,1); % Number of cells
    bin_size_1D = 0.15;
    cosEqn = 'cos(b*x-c)';
    phase_ctrs = (NP_units.phase_bins(1,1:end-1)+mean(diff(NP_units.phase_bins(1,:)))/2);
    
    for nc=1:N_cells
        
        %=== Rearrange Place Cell Data
        n_clus = size(NP_units.f_clus{nc,1},2)-1;
        NP_unitOnClus = cell(1,n_clus);
        for j = 1:n_clus
            NP_unitOnClus{1, j} = struct(); % Initialize as struct array
            
            %=== Store features of the unit within a cluster of interest
            NP_unitOnClus{1, j}.f_lenght = NP_units.f_clus{nc,1}(j+1).f_lenght;
            NP_unitOnClus{1, j}.f_duration = NP_units.f_clus{nc,1}(j+1).f_duration;
            NP_unitOnClus{1, j}.fr = NP_units.fr;
            NP_unitOnClus{1, j}.SI = NP_units.f_clus{nc,1}(j+1).SI_value;
            NP_unitOnClus{1, j}.p_val = NP_units.f_clus{nc,1}(j+1).SI_p_val;
            NP_unitOnClus{1, j}.spkPerflight = NP_units.f_clus{nc,1}(j+1).sum_spk/NP_units.f_clus{nc,1}(j+1).n;
            NP_unitOnClus{1, j}.plc_map = NP_units.f_clus{nc,1}(j+1).map(1,:)';
            NP_unitOnClus{1, j}.plc_ctr = NP_units.f_clus{nc,1}(j+1).binC';
            NP_unitOnClus{1, j}.prob_x = NP_units.f_clus{nc,1}(j+1).prob_x';
            NP_unitOnClus{1, j}.corr_m = NP_units.f_clus{nc,1}(j+1).corr_map_OdEv;
            NP_unitOnClus{1, j}.dist_m = NP_units.f_clus{nc,1}(j+1).dist_map_OdEv;
            NP_unitOnClus{1, j}.stab_m = NP_units.f_clus{nc,1}(j+1).corr_map_1h2h;
            NP_unitOnClus{1, j}.stbd_m = NP_units.f_clus{nc,1}(j+1).dist_map_1h2h;
            NP_unitOnClus{1, j}.peakHz = NP_units.f_clus{nc,1}(j+1).peakHz;
            NP_unitOnClus{1, j}.field_loc = NP_units.f_clus{nc,1}(j+1).field_loc;
            NP_unitOnClus{1, j}.field_loc_m = NP_units.f_clus{nc,1}(j+1).field_loc*bin_size_1D;
            NP_unitOnClus{1, j}.n_fields = NP_units.f_clus{nc,1}(j+1).n_fields;
            NP_unitOnClus{1, j}.f_width = NP_units.f_clus{nc,1}(j+1).f_width;
            NP_unitOnClus{1, j}.phase_max = NP_units.f_clus{nc,1}(j+1).phase_max;
            NP_unitOnClus{1, j}.sff = NP_units.f_clus{nc,1}(j+1).sff;
            NP_unitOnClus{1, j}.map_interp = NP_units.f_clus{nc,1}(j+1).map_interp;
            NP_unitOnClus{1, j}.asymm_frct = sum(NP_units.f_clus{nc,1}(j+1).map(1,1:round(NP_units.f_clus{nc,1}(j+1).field_loc)))./sum(NP_units.f_clus{nc,1}(j+1).map(1,:));
            NP_unitOnClus{1, j}.spk2wbt_rho = NP_units.f_clus{nc,1}(j+1).spk2wbt_rho;
            NP_unitOnClus{1, j}.spk2wbt_pvl = NP_units.f_clus{nc,1}(j+1).spk2wbt_pvl_sh;
            
            NP_unitOnClus{1, j}.place_cond =    NP_unitOnClus{1, j}.spkPerflight>1 &...  % Min Spikes per flight (DEF: 1)
                NP_unitOnClus{1, j}.peakHz>3 &...        % Min Peak Firing Rate (DEF: 3)
                NP_unitOnClus{1, j}.stab_m>0.4 &...      % Min Stability (DEF: .4)
                NP_unitOnClus{1, j}.sff<0.7;             % Min Peakyness (DEF: .7)
            
        end
        
        %=== Classify as place cell or not
        place_cond = zeros(n_clus,1);
        spkPerflight = zeros(n_clus,1);
        pp_corr = zeros(n_clus,1);
        pp_pval = zeros(n_clus,1);
        asymm_frct = zeros(n_clus,1);
        for j = 1:n_clus
            place_cond(j) = NP_unitOnClus{1, j}.place_cond;
            spkPerflight(j) = NP_unitOnClus{1, j}.spkPerflight;
            pp_corr(j) = NP_unitOnClus{1, j}.spk2wbt_rho;
            pp_pval(j) = NP_unitOnClus{1, j}.spk2wbt_pvl;
            asymm_frct(j) = NP_unitOnClus{1, j}.asymm_frct;
        end
        NP_units.place_cond(nc) = any(place_cond);
        NP_units.analz_cond(nc) = any(spkPerflight>1);
        NP_units.pp_analz_cond(nc) = any(~isnan(pp_corr));
        [min_val,min_loc] = min(pp_corr(place_cond & pp_pval<0.05));
        NP_units.phase_prec(nc) = ~isempty(pp_corr(place_cond & pp_pval<0.05));
        if ~isempty(min_val)
            NP_units.phs_prc_corr(nc) = min_val;
            NP_units.phs_prc_asym(nc) = asymm_frct(min_loc);
        else
            NP_units.phs_prc_corr(nc) = NaN;
            NP_units.phs_prc_asym(nc) = NaN;
        end
        
        %=== Add a few more features
        NP_units.num_spikes_f(nc) = sum(~isnan(NP_units.spk_phase2wbt_F{nc,1}));    % Number of spikes in flight
        NP_units.spk_phase2wbt_R2(nc) = NaN;
        NP_units.spk_phase2wbt_fit_par(nc,:) = NaN(1,2);
        NP_units.spk2wbt_p(nc) = NaN;
        NP_units.spk2wbt_pref_phase(nc) = NaN;
        NP_units.spk_phase2wbt_corrV(nc) = NaN;
        NP_units.spk_phase2wbt_corrP(nc) = NaN;
        NP_units.spk_phase2tmr_R2(nc) = NaN;
        NP_units.spk_phase2tmr_fit_par(nc,:) = NaN(1,2);
        NP_units.spk2tmr_p(nc) = NaN;
        NP_units.spk2tmr_pref_phase(nc) = NaN;
        NP_units.spk_phase2tmr_corrV(nc) = NaN;
        NP_units.spk_phase2tmr_corrP(nc) = NaN;
        if numel(NP_units.spk_phase2wbt_F{nc,1})>30
            
            %=== Fit with cosine function
            x = [phase_ctrs,phase_ctrs+2*pi];
            y = repmat(smoothdata(NP_units.spk_phase2wbt_F_counts{nc,:}','movmean',3),1,2);
            x1 = x(~isnan(y));  y1 = y(~isnan(y));
            y1 = normalize(y1-mean(y1),'range',[-1 1]);
            [f1,gof1] = fit(x1',y1',cosEqn,'Start',[2 0],'Lower',[-inf 0]);
            [f2,gof2] = fit(x1',y1',cosEqn,'Start',[1 0],'Lower',[-inf 0]);
%             [f1,gof1] = fit(x1',y1',cosEqn,'Start',[2 0],'Lower',[0 -pi]);
%             [f2,gof2] = fit(x1',y1',cosEqn,'Start',[1 0],'Lower',[0 -pi]);
            
            %=== Keep best fit
            if gof1.rsquare>gof2.rsquare
                NP_units.spk_phase2wbt_R2(nc) = gof1.rsquare;
                NP_units.spk_phase2wbt_fit_par(nc,:) = coeffvalues(f1);
                [NP_units.spk_phase2wbt_corrV(nc),NP_units.spk_phase2wbt_corrP(nc)] = corr(y1',f1(x1));
            else
                NP_units.spk_phase2wbt_R2(nc) = gof2.rsquare;
                NP_units.spk_phase2wbt_fit_par(nc,:) = coeffvalues(f2);
                [NP_units.spk_phase2wbt_corrV(nc),NP_units.spk_phase2wbt_corrP(nc)] = corr(y1',f2(x1));
            end
            NP_units.spk2wbt_p(nc) = circ_rtest(NP_units.spk_phase2wbt_F{nc,1});    % p value Rayleight test
            [~,max_loc] = max(smoothdata(NP_units.spk_phase2wbt_F_counts{nc,:},'movmean',3));
            NP_units.spk2wbt_pref_phase(nc) = phase_ctrs(max_loc);                  % Preferred phase
%                                 %=== Plotting (troubleshoot)
%                                 if gof1.rsquare>gof2.rsquare
%                                     plot(f1);    hold on;  bar(x1,y1,1,'EdgeColor','none'); hold off;   xticks()
%                                     title(gof1.rsquare);
%                                 else
%                                     plot(f2);    hold on;  bar(x1,y1,1,'EdgeColor','none'); hold off;
%                                     title(gof2.rsquare);
%                                 end
%                                   
%                                 choice = questdlg('Continue', 'Class Verification','Yes', 'No', 'Stop','Yes');
%                                 switch choice
%                                     case 'Stop'
%                                         break;
%                                 end
            
            
            % Fit also Tamir's phase distribution
            %=== Fit with cosine function
            x = [phase_ctrs,phase_ctrs+2*pi];
            y = repmat(smoothdata(NP_units.spk_phase2tmr_F_counts{nc,:}','movmean',3),1,2);
            x1 = x(~isnan(y));  y1 = y(~isnan(y));
            y1 = normalize(y1-mean(y1),'range',[-1 1]);
            [f1,gof1] = fit(x1',y1',cosEqn,'Start',[2 0],'Lower',[-inf 0]);
            [f2,gof2] = fit(x1',y1',cosEqn,'Start',[1 0],'Lower',[-inf 0]);
            
            %=== Keep best fit
            if gof1.rsquare>gof2.rsquare
                NP_units.spk_phase2tmr_R2(nc) = gof1.rsquare;
                NP_units.spk_phase2tmr_fit_par(nc,:) = coeffvalues(f1);
                [NP_units.spk_phase2tmr_corrV(nc),NP_units.spk_phase2tmr_corrP(nc)] = corr(y1',f1(x1));
            else
                NP_units.spk_phase2tmr_R2(nc) = gof2.rsquare;
                NP_units.spk_phase2tmr_fit_par(nc,:) = coeffvalues(f2);
                [NP_units.spk_phase2tmr_corrV(nc),NP_units.spk_phase2tmr_corrP(nc)] = corr(y1',f2(x1));
            end
            NP_units.spk2tmr_p(nc) = circ_rtest(NP_units.spk_phase2tmr_F{nc,1});    % p value Rayleight test
            [~,max_loc] = max(smoothdata(NP_units.spk_phase2tmr_F_counts{nc,:},'movmean',3));
            NP_units.spk2tmr_pref_phase(nc) = phase_ctrs(max_loc);
            
        end
        
    end
    
end

%% LOOK AT PHASE LOCKING USING CIRCULAR STATISTICS
for hide=1

    %=== Add group ID to each cell and 
    [~,~,NP_units.groupID] =  unique(string(NP_units.unique_ID),'rows');
    NP_units.p_val_r_spk_phase2wbt_omni = cellfun(@circ_otest,NP_units.spk_phase2wbt_F);
    NP_units.p_val_r_spk_phase2tmr_omni = cellfun(@circ_otest,NP_units.spk_phase2tmr_F);
    NP_units.n_spikes_tmr = cellfun(@numel,NP_units.spk_phase2tmr_F);
    
    %=== Define Subset of Interest
    %NP_sst = NP_units(NP_units.fr>0 & NP_units.nspkpf>0 & NP_units.num_spikes_f>50,:);
    NP_sst = NP_units(NP_units.fr>0 & NP_units.nspkpf>0 & NP_units.place_cond & NP_units.num_spikes_f>30,:);
    N_cells = size(NP_sst,1); % Number of cells
    
    %=== Add condition on Phase Locked Neurons
    NP_sst.wbt_lock1 = NP_sst.p_val_r_spk_phase2wbt_boots<0.05;
    NP_sst.tmr_lock1 = NP_sst.p_val_r_spk_phase2tmr_boots<0.05;
    NP_sst.wbt_lock2 = NP_sst.p_val_r_spk_phase2wbt_omni<0.05;
    NP_sst.tmr_lock2 = NP_sst.p_val_r_spk_phase2tmr_omni<0.05;
    NP_sst.wbt_lock3 = NP_sst.spk_phase2wbt_R2>0;
    NP_sst.tmr_lock3 = NP_sst.spk_phase2tmr_R2>0;
    display([sum(NP_sst.wbt_lock1)./N_cells,sum(NP_sst.tmr_lock1)./N_cells]);
    display([sum(NP_sst.wbt_lock2)./N_cells,sum(NP_sst.tmr_lock2)./N_cells]);
    display([sum(NP_sst.wbt_lock3)./N_cells,sum(NP_sst.tmr_lock3)./N_cells]);
    
    %=== Calculate angular standard deviation
    NP_sst.ang_std_spk_phase2wbt = cellfun(@circ_std,NP_sst.spk_phase2wbt_F);
    NP_sst.ang_std_spk_phase2tmr = cellfun(@circ_std,NP_sst.spk_phase2tmr_F);
    
    %=== Phase locking at the population level
    mrv_wbt1 = [];   mrv_tmr1 = [];
    mrv_wbt2 = [];   mrv_tmr2 = [];
    mrv_wbt3 = [];   mrv_tmr3 = [];
    sps_wbt1 = [];   sps_tmr1 = [];
    sps_wbt2 = [];   sps_tmr2 = [];
    sps_wbt3 = [];   sps_tmr3 = [];
    for i=unique(NP_sst.groupID)'
        all_wbt_phases = vertcat(NP_sst.spk_phase2wbt_F{NP_sst.groupID==i & NP_sst.wbt_lock1});
        all_tmr_phases = vertcat(NP_sst.spk_phase2tmr_F{NP_sst.groupID==i & NP_sst.tmr_lock1});
        if ~isempty(all_wbt_phases(~isnan(all_wbt_phases)))
            mrv_wbt1 = [mrv_wbt1; circ_r(all_wbt_phases(~isnan(all_wbt_phases)))];
            mrv_tmr1 = [mrv_tmr1; circ_r(all_tmr_phases(~isnan(all_tmr_phases)))];
            sps_wbt1 = [sps_wbt1; numel(all_wbt_phases(~isnan(all_wbt_phases)))];
            sps_tmr1 = [sps_tmr1; numel(all_tmr_phases(~isnan(all_tmr_phases)))];
        end
        
        all_wbt_phases = vertcat(NP_sst.spk_phase2wbt_F{NP_sst.groupID==i & NP_sst.wbt_lock2});
        all_tmr_phases = vertcat(NP_sst.spk_phase2tmr_F{NP_sst.groupID==i & NP_sst.tmr_lock2});
        if ~isempty(all_wbt_phases(~isnan(all_wbt_phases)))
            mrv_wbt2 = [mrv_wbt2; circ_r(all_wbt_phases(~isnan(all_wbt_phases)))];
            mrv_tmr2 = [mrv_tmr2; circ_r(all_tmr_phases(~isnan(all_tmr_phases)))];
            sps_wbt2 = [sps_wbt2; numel(all_wbt_phases(~isnan(all_wbt_phases)))];
            sps_tmr2 = [sps_tmr2; numel(all_tmr_phases(~isnan(all_tmr_phases)))];
        end
        
        all_wbt_phases = vertcat(NP_sst.spk_phase2wbt_F{NP_sst.groupID==i & NP_sst.wbt_lock3});
        all_tmr_phases = vertcat(NP_sst.spk_phase2tmr_F{NP_sst.groupID==i & NP_sst.tmr_lock3});
        if ~isempty(all_wbt_phases(~isnan(all_wbt_phases)))
            mrv_wbt3 = [mrv_wbt3; circ_r(all_wbt_phases(~isnan(all_wbt_phases)))];
            mrv_tmr3 = [mrv_tmr3; circ_r(all_tmr_phases(~isnan(all_tmr_phases)))];
            sps_wbt3 = [sps_wbt3; numel(all_wbt_phases(~isnan(all_wbt_phases)))];
            sps_tmr3 = [sps_tmr3; numel(all_tmr_phases(~isnan(all_tmr_phases)))];
        end
    end
    
    %=== Quantify fraction of locked cells and mean vector lenght
    Fraction_summary = groupsummary(NP_sst,'groupID','mean',{'wbt_lock1','tmr_lock1','wbt_lock2','tmr_lock2','wbt_lock3','tmr_lock3'});
    Wbt_summary1 = groupsummary(NP_sst(NP_sst.wbt_lock1,:),'groupID','mean',{'r_spk_phase2wbt','m_spk_phase2wbt','ang_std_spk_phase2wbt','z_spk_phase2wbt','z_spk_phase2tmr','n_spk_phase2wbt'});
    Tmr_summary1 = groupsummary(NP_sst(NP_sst.tmr_lock1,:),'groupID','mean',{'r_spk_phase2tmr','m_spk_phase2tmr','ang_std_spk_phase2tmr','z_spk_phase2wbt','z_spk_phase2tmr','n_spk_phase2wbt'});
    Wbt_summary2 = groupsummary(NP_sst(NP_sst.wbt_lock2,:),'groupID','mean',{'r_spk_phase2wbt','m_spk_phase2wbt','ang_std_spk_phase2wbt','z_spk_phase2wbt','z_spk_phase2tmr','n_spk_phase2wbt'});
    Tmr_summary2 = groupsummary(NP_sst(NP_sst.tmr_lock2,:),'groupID','mean',{'r_spk_phase2tmr','m_spk_phase2tmr','ang_std_spk_phase2tmr','z_spk_phase2wbt','z_spk_phase2tmr','n_spk_phase2wbt'});
    Wbt_summary3 = groupsummary(NP_sst(NP_sst.wbt_lock3,:),'groupID','mean',{'r_spk_phase2wbt','m_spk_phase2wbt','ang_std_spk_phase2wbt','z_spk_phase2wbt','z_spk_phase2tmr','n_spk_phase2wbt'});
    Tmr_summary3 = groupsummary(NP_sst(NP_sst.tmr_lock3,:),'groupID','mean',{'r_spk_phase2tmr','m_spk_phase2tmr','ang_std_spk_phase2tmr','z_spk_phase2wbt','z_spk_phase2tmr','n_spk_phase2wbt'});

    %=== Display relevant numbers
    mean([Fraction_summary.mean_wbt_lock2, Fraction_summary.mean_tmr_lock2])
    mean_and_sem_AF_v0(Fraction_summary.mean_wbt_lock2)
    mean_and_sem_AF_v0(Wbt_summary2.mean_r_spk_phase2wbt)
    mean_and_sem_AF_v0(Wbt_summary2.mean_z_spk_phase2wbt)
    mean_and_sem_AF_v0(Wbt_summary2.mean_n_spk_phase2wbt)
    
    mean([Fraction_summary.mean_wbt_lock1, Fraction_summary.mean_tmr_lock1])
    mean_and_sem_AF_v0(Wbt_summary1.mean_r_spk_phase2wbt)
    mean_and_sem_AF_v0(Wbt_summary1.mean_z_spk_phase2wbt)
    mean_and_sem_AF_v0(Wbt_summary1.mean_n_spk_phase2wbt)
    
    mean(Fraction_summary.mean_wbt_lock3)
    
    %=== Look at the dispersion of the prefferred phase within a session
    wbt_std = [];   tmr_std = [];
    for i=unique(NP_sst.groupID)'
       NP_session_wbt = NP_sst(NP_sst.groupID==i & NP_sst.wbt_lock2,:);
       NP_session_tmr = NP_sst(NP_sst.groupID==i & NP_sst.tmr_lock2,:);
       if size(NP_session_wbt,1)>1 && size(NP_session_tmr,1)>1
           wbt_std = [wbt_std; circ_std(NP_session_wbt.m_spk_phase2wbt)];
           tmr_std = [tmr_std; circ_std(NP_session_tmr.m_spk_phase2tmr)];
       end
    end
    figure('units','normalized','outerposition',[.2 .3 .15 .3]);
    plot_distr_AF_v0(wbt_std , tmr_std, {'Wingbeat', 'Non-oscillatory'}, 'SEM', 'Within Session Angular Standard Deviation'); %ylim_vals = ylim;   ylim([0 ylim_vals(2)]);
    title(signrank(wbt_std, tmr_std))
    
    %=== Accumulate Phase distributions
    spk2wbt_counts = zeros(numel(NP_sst.spk_phase2wbt_F_counts{1,:}),N_cells);
    all2wbt_counts = zeros(numel(NP_sst.spk_phase2wbt_F_counts{1,:}),N_cells);
    spk2wbt_nrmcts = zeros(numel(NP_sst.spk_phase2wbt_F_counts{1,:}),N_cells);
    spk2tmr_counts = zeros(numel(NP_sst.spk_phase2tmr_F_counts{1,:}),N_cells);
    all2tmr_counts = zeros(numel(NP_sst.spk_phase2tmr_F_counts{1,:}),N_cells);
    spk2tmr_nrmcts = zeros(numel(NP_sst.spk_phase2tmr_F_counts{1,:}),N_cells);
    phase_bins = NP_sst.phase_bins(1,:);
    phase_ctrs = NP_sst.phase_bins(1,1:end-1)+mean(diff(NP_sst.phase_bins(1,:)))/2;
    for nc=1:N_cells
        
        spk2wbt_counts(:,nc) = NP_sst.spk_phase2wbt_F_counts{nc,:};
        all2wbt_counts(:,nc) = NP_sst.all_phase2wbt_counts{nc,:};
        spk2wbt_nrmcts(:,nc) = NP_sst.spk_phase2wbt_F_counts{nc,:}-NP_sst.all_phase2wbt_counts{nc,:};
        spk2tmr_counts(:,nc) = NP_sst.spk_phase2tmr_F_counts{nc,:};
        all2tmr_counts(:,nc) = NP_sst.all_phase2tmr_counts{nc,:};
        spk2tmr_nrmcts(:,nc) = NP_sst.spk_phase2tmr_F_counts{nc,:}-NP_sst.all_phase2tmr_counts{nc,:};
        
    end
    
    %=== Define Phase Locked Neurons
    %wbt_locked_cells = find(NP_sst.p_val_r_spk_phase2wbt_omni<0.05);
    %tmr_locked_cells = find(NP_sst.p_val_r_spk_phase2tmr_omni<0.05);
    wbt_locked_cells = find(NP_sst.wbt_lock2);
    tmr_locked_cells = find(NP_sst.tmr_lock2);
    
    %=== Show distribution of mean vector lenghts and phase preferences
%     figure('units','normalized','outerposition',[.2 .3 .15 .4]);
%     plot_distr_AF_v0(NP_sst.r_spk_phase2wbt(wbt_locked_cells), NP_sst.r_spk_phase2tmr(tmr_locked_cells), {'Wingbeat', 'Non-oscillatory'}, 'SEM', 'Mean resultant vector length (significant)'); ylim_vals = ylim;   ylim([0 ylim_vals(2)]);
%     
    
    %=== Plot some examples of phase locked cells
    figure('units','normalized','outerposition',[0 0 .5 1]);
    tiledlayout(10,10,'TileSpacing','tight');
    [~,tmp_idx] = sort(NP_sst.r_spk_phase2wbt(wbt_locked_cells),'descend');
    wbt_locked_cells2plot = wbt_locked_cells(tmp_idx);
    for i=1:min(100,numel(wbt_locked_cells2plot))
        nexttile;   polarhistogram(NP_sst.spk_phase2wbt_F{wbt_locked_cells2plot(i),:},18,'facealpha',.5,'edgecolor','none','FaceColor','k','Normalization','probability'); 
        rlim_tmp = rlim;    hold on;
        polarplot(circ_mean(NP_sst.spk_phase2wbt_F{wbt_locked_cells2plot(i),:}).*[1 1], [0 1], 'r', 'LineWidth', 3);
        rlim(rlim_tmp);         rticks([]);   thetaticklabels([]);  thetaticks([0 90 180 270]);
    end
    %=== Plot some examples of phase locked cells
    figure('units','normalized','outerposition',[.5 0 .5 1]);
    tiledlayout(10,10,'TileSpacing','tight');
    [~,tmp_idx] = sort(NP_sst.r_spk_phase2tmr(tmr_locked_cells),'descend');
    tmr_locked_cells2plot = tmr_locked_cells(tmp_idx);
    for i=1:min(100,numel(tmr_locked_cells2plot))
        nexttile;   polarhistogram(NP_sst.spk_phase2tmr_F{tmr_locked_cells2plot(i),:},18,'facealpha',.5,'edgecolor','none','FaceColor','b','Normalization','probability'); 
        rlim_tmp = rlim;    hold on;
        polarplot(circ_mean(NP_sst.spk_phase2tmr_F{tmr_locked_cells2plot(i),:}).*[1 1], [0 1], 'r', 'LineWidth', 3);
        rlim(rlim_tmp);         rticks([]);   thetaticklabels([]);  thetaticks([0 90 180 270]);
    end
    
    %=====================================================================================================================================================================
    
    %=== Look at spike-tgd LFP and Wingbeat from phase locked cells
    figure('units','normalized','outerposition',[0.3 0.2 0.45 0.3]);
    tiledlayout(1,3,'TileSpacing','tight');
    nexttile;   data = horzcat(NP_sst.spkTgdLfP_rest{find(NP_sst.tmr_lock3),:});
    data = highpass(data,2,500)';
    plotWinterval_AF_v0(LFP_table.spkTgdLfP_time{1,1},mean(data),mean(data)-std(data)./sqrt(size(data,1)),mean(data)+std(data)./sqrt(size(data,1)),'k');
    xlabel('Time from spike (s)');  ylabel('LFP (uV)');         xlim('tight');  plot([0 0],ylim,'k--'); title('Rest');  xlim([-2 2]);
    nexttile;   data = horzcat(NP_sst.spkTgdLfP_flight{find(NP_sst.tmr_lock3),:});
    data = highpass(data,2,500)';
    plotWinterval_AF_v0(LFP_table.spkTgdLfP_time{1,1},mean(data),mean(data)-std(data)./sqrt(size(data,1)),mean(data)+std(data)./sqrt(size(data,1)),'k');
    xlabel('Time from spike (s)');  ylabel('LFP (uV)');         xlim('tight');  plot([0 0],ylim,'k--');  title('Flight');   xlim([-2 2]);
    nexttile;   data = horzcat(NP_sst.spkTgdWbt_flight{find(NP_sst.wbt_lock3),:});
    data = highpass(data,2,500)';
    plotWinterval_AF_v0(LFP_table.spkTgdLfP_time{1,1},mean(data),mean(data)-std(data)./sqrt(size(data,1)),mean(data)+std(data)./sqrt(size(data,1)),'r');
    xlabel('Time from spike (s)');  ylabel('Acceleration (g)'); xlim('tight');  plot([0 0],ylim,'k--');   title('Flight');  xlim([-2 2]);
    
    st_LFR = [];    st_LFP = [];    st_WBT = [];
    for i=unique(NP_sst.groupID)'
       NP_session_wbt = NP_sst(NP_sst.groupID==i & NP_sst.wbt_lock3,:);
       NP_session_tmr = NP_sst(NP_sst.groupID==i & NP_sst.tmr_lock3,:);
       if size(NP_session_wbt,1)>1 && size(NP_session_tmr,1)>1
           st_LFR = [st_LFR ,median(highpass(horzcat(NP_session_tmr.spkTgdLfP_rest{:}),2,500),2)];
           st_LFP = [st_LFP ,median(highpass(horzcat(NP_session_tmr.spkTgdLfP_flight{:}),2,500),2)];
           st_WBT = [st_WBT ,median(highpass(horzcat(NP_session_wbt.spkTgdWbt_flight{:}),2,500),2)];
       end
    end
    figure('units','normalized','outerposition',[0.3 0.2 0.45 0.3]);
    tiledlayout(1,3,'TileSpacing','tight');
    nexttile;   data = st_LFR';
    plotWinterval_AF_v0(LFP_table.spkTgdLfP_time{1,1},mean(data),mean(data)-std(data)./sqrt(size(data,1)),mean(data)+std(data)./sqrt(size(data,1)),'k');
    xlabel('Time from spike (s)');  ylabel('LFP (uV)');         xlim('tight');  plot([0 0],ylim,'k--'); title('Rest');
    nexttile;   data = st_LFP';
    plotWinterval_AF_v0(LFP_table.spkTgdLfP_time{1,1},mean(data),mean(data)-std(data)./sqrt(size(data,1)),mean(data)+std(data)./sqrt(size(data,1)),'k');
    xlabel('Time from spike (s)');  ylabel('LFP (uV)');         xlim('tight');  plot([0 0],ylim,'k--');  title('Flight');
    nexttile;   data = st_WBT';
    plotWinterval_AF_v0(LFP_table.spkTgdLfP_time{1,1},mean(data),mean(data)-std(data)./sqrt(size(data,1)),mean(data)+std(data)./sqrt(size(data,1)),'r');
    xlabel('Time from spike (s)');  ylabel('Acceleration (g)'); xlim('tight');  plot([0 0],ylim,'k--');   title('Flight');
    
    %=== Pool data from all the cells
    figure('units','normalized','outerposition',[0.3 0.2 0.45 0.3]);
    tiledlayout(1,3,'TileSpacing','tight');
    nexttile;   data = horzcat(NP_sst.spkTgdLfP_rest{:});
    data = highpass(data,2,500)';
    plotWinterval_AF_v0(LFP_table.spkTgdLfP_time{1,1},mean(data),mean(data)-std(data)./sqrt(size(data,1)),mean(data)+std(data)./sqrt(size(data,1)),'k');
    xlabel('Time from spike (s)');  ylabel('LFP (uV)');         xlim('tight');  plot([0 0],ylim,'k--'); title('Rest');
    nexttile;   data = horzcat(NP_sst.spkTgdLfP_flight{:});
    data = highpass(data,2,500)';
    plotWinterval_AF_v0(LFP_table.spkTgdLfP_time{1,1},mean(data),mean(data)-std(data)./sqrt(size(data,1)),mean(data)+std(data)./sqrt(size(data,1)),'k');
    xlabel('Time from spike (s)');  ylabel('LFP (uV)');         xlim('tight');  plot([0 0],ylim,'k--');  title('Flight');
    nexttile;   data = horzcat(NP_sst.spkTgdWbt_flight{:});
    data = highpass(data,2,500)';
    plotWinterval_AF_v0(LFP_table.spkTgdLfP_time{1,1},mean(data),mean(data)-std(data)./sqrt(size(data,1)),mean(data)+std(data)./sqrt(size(data,1)),'r');
    xlabel('Time from spike (s)');  ylabel('Acceleration (g)'); xlim('tight');  plot([0 0],ylim,'k--');   title('Flight');
    
end

%% THETA-SWEEPS AND WINGBEAT
for hide=1
    
    %=== Select subtable of good flights
    warning('off');
    SWPs_sst = SWPs(SWPs.rmsDec_error<1.3 & SWPs.prc_decoded>0.7,:);
    
    %=== Params and initialize relevant variables
    bin_size_1D = 0.15;             % Spatial bin size
    n_swp_shuffles = 20;             % Default 20
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
        'rsz_spk_dsty', [], 'wbt_power', [],'tmr_power', [], 'rsz_LFP', [], 'rsz_wbt', [], ...
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
            
            %=== Raw, filtered, shifted and rehaped posterior
            single_SWP_struct(ss).raw_posterior = {SWPs_sst.p_dec_flight{zz,1}(:,swp_interval)};
            single_SWP_struct(ss).sft_posterior = {SWPs_sst.p_dec_shifted{zz,1}(:,swp_interval)};
            single_SWP_struct(ss).rsp_posterior = {imresize(single_SWP_struct(ss).sft_posterior{:},[size(single_SWP_struct(ss).sft_posterior{:},1),n_reshape])};
            
            %=== Spike density, wingbeat, LFP, Tamir's phase echolocation and phase
            single_SWP_struct(ss).rsz_spk_dsty = {interp1(SWPs_sst.spk_dsty{zz,1}(swp_interval),linspace(1,numel(SWPs_sst.spk_dsty{zz,1}(swp_interval)),n_reshape)')};
            single_SWP_struct(ss).wbt_power = {SWPs_sst.wbt_power{zz,1}(swp_interval)};
            single_SWP_struct(ss).tmr_power = mean(SWPs_sst.tmr_power{zz,1}(swp_interval));
            single_SWP_struct(ss).rsz_LFP = {interp1(SWPs_sst.LFP{zz,1}(swp_interval),linspace(1,numel(SWPs_sst.LFP{zz,1}(swp_interval)),n_reshape)')};
            single_SWP_struct(ss).rsz_wbt = {interp1(SWPs_sst.wbt{zz,1}(swp_interval),linspace(1,numel(SWPs_sst.wbt{zz,1}(swp_interval)),n_reshape)')};
            single_SWP_struct(ss).fract_pos = {SWPs_sst.pos_real{zz,1}(swp_interval)/SWPs_sst.pos_real{zz,1}(end)};
            single_SWP_struct(ss).rsz_clk = {interp1(SWPs_sst.clk{zz,1}(swp_interval),linspace(1,numel(SWPs_sst.clk{zz,1}(swp_interval)),n_reshape)')};
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
    
    %=== Real data, Cut at tamir maxima
    single_SWP_tm = table();
    for zz=1:N_flights
        
        %=== Find wingbeat minima
        zero_phs_idx = find(SWPs_sst.tmr_phase{zz,1}(1:end-1).* SWPs_sst.tmr_phase{zz,1}(2:end)<0 & diff(SWPs_sst.tmr_phase{zz,1})<0);  % Segment based on wingbeat
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
            
            %=== Raw, filtered, shifted and rehaped posterior
            single_SWP_struct(ss).raw_posterior = {SWPs_sst.p_dec_flight{zz,1}(:,swp_interval)};
            single_SWP_struct(ss).sft_posterior = {SWPs_sst.p_dec_shifted{zz,1}(:,swp_interval)};
            single_SWP_struct(ss).rsp_posterior = {imresize(single_SWP_struct(ss).sft_posterior{:},[size(single_SWP_struct(ss).sft_posterior{:},1),n_reshape])};
            
            %=== Spike density, wingbeat, LFP, Tamir's phase echolocation and phase
            single_SWP_struct(ss).rsz_spk_dsty = {interp1(SWPs_sst.spk_dsty{zz,1}(swp_interval),linspace(1,numel(SWPs_sst.spk_dsty{zz,1}(swp_interval)),n_reshape)')};
            single_SWP_struct(ss).wbt_power = {SWPs_sst.wbt_power{zz,1}(swp_interval)};
            single_SWP_struct(ss).tmr_power = mean(SWPs_sst.tmr_power{zz,1}(swp_interval));
            single_SWP_struct(ss).rsz_LFP = {interp1(SWPs_sst.LFP{zz,1}(swp_interval),linspace(1,numel(SWPs_sst.LFP{zz,1}(swp_interval)),n_reshape)')};
            single_SWP_struct(ss).rsz_wbt = {interp1(SWPs_sst.wbt{zz,1}(swp_interval),linspace(1,numel(SWPs_sst.wbt{zz,1}(swp_interval)),n_reshape)')};
            single_SWP_struct(ss).fract_pos = {SWPs_sst.pos_real{zz,1}(swp_interval)/SWPs_sst.pos_real{zz,1}(end)};
            single_SWP_struct(ss).rsz_clk = {interp1(SWPs_sst.clk{zz,1}(swp_interval),linspace(1,numel(SWPs_sst.clk{zz,1}(swp_interval)),n_reshape)')};
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
        single_SWP_tm = [single_SWP_tm; struct2table(single_SWP_struct)];
        
    end
    
    %=== Shuffled data, Cut at random points
    rng(1);
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
    %min_tmr_pwr = prctile([single_SWP_rl.tmr_power;single_SWP_tm.tmr_power],25);
    min_tmr_pwr = prctile(vertcat(SWPs_sst.tmr_power{:}),25);
    
    %=== Extract wingbeat and Tamir's set
    SWP_sst_rl = single_SWP_rl(single_SWP_rl.mean_fract_pos>min_pos1 & single_SWP_rl.mean_fract_pos<min_pos2 & single_SWP_rl.rms_dec<min_rms &...
        single_SWP_rl.prc_dec>min_prc & single_SWP_rl.mean_jmp_distance>min_mjp & single_SWP_rl.med_max_post>min_mp & single_SWP_rl.tmr_power>min_tmr_pwr,:);
    SWP_sst_tm = single_SWP_tm(single_SWP_tm.mean_fract_pos>min_pos1 & single_SWP_tm.mean_fract_pos<min_pos2 & single_SWP_tm.rms_dec<min_rms &...
        single_SWP_tm.prc_dec>min_prc & single_SWP_tm.mean_jmp_distance>min_mjp & single_SWP_tm.med_max_post>min_mp & single_SWP_tm.tmr_power>min_tmr_pwr,:);
    
    %=== Match the same flights
    match_flight_ID = intersect(SWP_sst_rl.flight_ID,SWP_sst_tm.flight_ID);
    SWP_sst_rl = SWP_sst_rl(arrayfun(@(x) any(x == match_flight_ID),SWP_sst_rl.flight_ID),:);
    SWP_sst_tm = SWP_sst_tm(arrayfun(@(x) any(x == match_flight_ID),SWP_sst_tm.flight_ID),:);
    N_sweeps_rl = size(SWP_sst_rl,1);
    N_sweeps_tm = size(SWP_sst_tm,1);
    
    %=== Calculate averages and STDs (REAL DATA, wingbeat)
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
    
    %=== Calculate averages and STDs (REAL DATA, Tamir's)
    est_dist_tm = zeros(n_reshape,N_sweeps_tm);
    rsz_spk_dsty_tm = zeros(n_reshape,N_sweeps_tm);
    rsz_wbt_tm = zeros(n_reshape,N_sweeps_tm);
    rsz_LFP_tm = zeros(n_reshape,N_sweeps_tm);
    rsz_clk_tm = zeros(n_reshape,N_sweeps_tm);
    rsz_wbt_phase_tm = zeros(n_reshape,N_sweeps_tm);
    for i=1:N_sweeps_tm
        est_dist_tm(:,i) = SWP_sst_tm.est_dist1{i,1};
        rsz_spk_dsty_tm(:,i) = SWP_sst_tm.rsz_spk_dsty{i,1};
        rsz_wbt_tm(:,i) = SWP_sst_tm.rsz_wbt{i,1};
        rsz_LFP_tm(:,i) = SWP_sst_tm.rsz_LFP{i,1};
        rsz_clk_tm(:,i) = SWP_sst_tm.rsz_clk{i,1};
        rsz_wbt_phase_tm(:,i) = SWP_sst_tm.rsz_wbt_phase{i,1};
    end
    avg_est_dist_tm = mean(est_dist_tm,2);            sem_est_dist_tm = std(est_dist_tm,[],2)/sqrt(N_sweeps_tm);
    avg_rsz_spk_dsty_tm = mean(rsz_spk_dsty_tm,2);    sem_rsz_spk_dsty_tm = std(rsz_spk_dsty_tm,[],2)/sqrt(N_sweeps_tm);
    avg_rsz_wbt_tm = mean(rsz_wbt_tm,2);              sem_rsz_wbt_tm = std(rsz_wbt_tm,[],2)/sqrt(N_sweeps_tm);
    avg_rsz_LFP_tm = mean(rsz_LFP_tm,2);              sem_rsz_LFP_tm = std(rsz_LFP_tm,[],2)/sqrt(N_sweeps_tm);
    avg_rsz_wbt_phase_tm = mean(rsz_wbt_phase_tm,2);  sem_rsz_wbt_phase_tm = std(rsz_wbt_phase_tm,[],2)/sqrt(N_sweeps_tm);
    avg_rsz_clk_tm = mean(rsz_clk_tm,2);              sem_rsz_clk_tm = std(rsz_clk_tm,[],2)/sqrt(N_sweeps_tm);
    
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
            single_SWP_sh.prc_dec>min_prc & single_SWP_sh.mean_jmp_distance>min_mjp & single_SWP_sh.med_max_post>min_mp,:);
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
    pVal_est_dist_rl = sum(avg_est_dist_rl<avg_est_dist_sh_all,2)./n_swp_shuffles;
    pVal_est_dist_tm = sum(avg_est_dist_tm<avg_est_dist_sh_all,2)./n_swp_shuffles;
    
    % === Show averages 
    phlim = [-180 180];
    
    %=== Show averages (Wingbeat)
    figure('units','normalized','outerposition',[.3 .4 .45 .3]);
    tiledlayout(1,6,'TileSpacing','compact');
    nexttile;   hold on;
    plotWinterval_AF_v0(rsz_bin_ctrs,avg_est_dist_sh,avg_est_dist_sh-sem_est_dist_sh,avg_est_dist_sh+sem_est_dist_sh,'k');
    plotWinterval_AF_v0(rsz_bin_ctrs,avg_est_dist_rl,avg_est_dist_rl-sem_est_dist_rl,avg_est_dist_rl+sem_est_dist_rl,[.9 .5 .2]);
    plot(rsz_bin_ctrs(pVal_est_dist_rl<0.05),max(avg_est_dist_rl)*ones(size(rsz_bin_ctrs(pVal_est_dist_rl<0.05))),'*');
    xticks(sort([phlim, 0]));    title('Average Decoding Error');    plot(0*[1 1],ylim,'k--'); xlim('tight');   xlim(phlim);
    nexttile;   hold on;
    plotWinterval_AF_v0(rsz_bin_ctrs,avg_rsz_spk_dsty_sh,avg_rsz_spk_dsty_sh-sem_rsz_spk_dsty_sh,avg_rsz_spk_dsty_sh+sem_rsz_spk_dsty_sh,'k');
    plotWinterval_AF_v0(rsz_bin_ctrs,avg_rsz_spk_dsty_rl,avg_rsz_spk_dsty_rl-sem_rsz_spk_dsty_rl,avg_rsz_spk_dsty_rl+sem_rsz_spk_dsty_rl,[.9 .5 .2]);
    xticks(sort([phlim, 0]));    title('Average Spike Density');    plot(0*[1 1],ylim,'k--'); xlim('tight');    xlim(phlim);
    nexttile;   hold on;
    plotWinterval_AF_v0(rsz_bin_ctrs,avg_rsz_LFP_sh,avg_rsz_LFP_sh-sem_rsz_LFP_sh,avg_rsz_LFP_sh+sem_rsz_LFP_sh,'k');
    plotWinterval_AF_v0(rsz_bin_ctrs,avg_rsz_LFP_rl,avg_rsz_LFP_rl-sem_rsz_LFP_rl,avg_rsz_LFP_rl+sem_rsz_LFP_rl,[.9 .5 .2]);
    xticks(sort([phlim, 0]));    title('Average LFP');    plot(0*[1 1],ylim,'k--'); xlim('tight');              xlim(phlim);
    nexttile;   hold on;
    plotWinterval_AF_v0(rsz_bin_ctrs,avg_rsz_clk_sh,avg_rsz_clk_sh-sem_rsz_clk_sh,avg_rsz_clk_sh+sem_rsz_clk_sh,'k');
    plotWinterval_AF_v0(rsz_bin_ctrs,avg_rsz_clk_rl,avg_rsz_clk_rl-sem_rsz_clk_rl,avg_rsz_clk_rl+sem_rsz_clk_rl,[.9 .5 .2]);
    xticks(sort([phlim, 0]));    title('Average Call Rate');    plot(0*[1 1],ylim,'k--'); xlim('tight');        xlim(phlim);
    nexttile;   hold on;
    plotWinterval_AF_v0(rsz_bin_ctrs,avg_rsz_wbt_phase_sh,avg_rsz_wbt_phase_sh-sem_rsz_wbt_phase_sh,avg_rsz_wbt_phase_sh+sem_rsz_wbt_phase_sh,'k');
    plotWinterval_AF_v0(rsz_bin_ctrs,avg_rsz_wbt_phase_rl,avg_rsz_wbt_phase_rl-sem_rsz_wbt_phase_rl,avg_rsz_wbt_phase_rl+sem_rsz_wbt_phase_rl,[.9 .5 .2]);
    xticks(sort([phlim, 0]));    title('Wingbeat Phase');    plot(0*[1 1],ylim,'k--'); xlim('tight');           xlim(phlim);
    nexttile;   hold on;
    plotWinterval_AF_v0(rsz_bin_ctrs,avg_rsz_wbt_sh,avg_rsz_wbt_sh-sem_rsz_wbt_sh,avg_rsz_wbt_sh+sem_rsz_wbt_sh,'k');
    plotWinterval_AF_v0(rsz_bin_ctrs,avg_rsz_wbt_rl,avg_rsz_wbt_rl-sem_rsz_wbt_rl,avg_rsz_wbt_rl+sem_rsz_wbt_rl,[.9 .5 .2]);
    xticks(sort([phlim, 0]));    title('Accelerometer');    plot(0*[1 1],ylim,'k--'); xlim('tight');           xlim(phlim);
    sgtitle(['Average of ',num2str(N_sweeps_rl),' cycles, from ',num2str(numel(unique(SWP_sst_rl.flight_ID))),' flights']);
    
    %=== Show averages (Tamir)
    figure('units','normalized','outerposition',[.3 .1 .45 .3]);
    tiledlayout(1,6,'TileSpacing','compact');
    nexttile;   hold on;
    plotWinterval_AF_v0(rsz_bin_ctrs,avg_est_dist_sh,avg_est_dist_sh-sem_est_dist_sh,avg_est_dist_sh+sem_est_dist_sh,'k');
    plotWinterval_AF_v0(rsz_bin_ctrs,avg_est_dist_tm,avg_est_dist_tm-sem_est_dist_tm,avg_est_dist_tm+sem_est_dist_tm,[.5 .0 .5]);
    plot(rsz_bin_ctrs(pVal_est_dist_tm<0.05),max(avg_est_dist_tm)*ones(size(rsz_bin_ctrs(pVal_est_dist_tm<0.05))),'*');
    xticks(sort([phlim, 0]));    title('Average Decoding Error');    plot(0*[1 1],ylim,'k--'); xlim('tight');   xlim(phlim);
    nexttile;   hold on;
    plotWinterval_AF_v0(rsz_bin_ctrs,avg_rsz_spk_dsty_sh,avg_rsz_spk_dsty_sh-sem_rsz_spk_dsty_sh,avg_rsz_spk_dsty_sh+sem_rsz_spk_dsty_sh,'k');
    plotWinterval_AF_v0(rsz_bin_ctrs,avg_rsz_spk_dsty_tm,avg_rsz_spk_dsty_tm-sem_rsz_spk_dsty_tm,avg_rsz_spk_dsty_tm+sem_rsz_spk_dsty_tm,[.5 .0 .5]);
    xticks(sort([phlim, 0]));    title('Average Spike Density');    plot(0*[1 1],ylim,'k--'); xlim('tight');    xlim(phlim);
    nexttile;   hold on;
    plotWinterval_AF_v0(rsz_bin_ctrs,avg_rsz_LFP_sh,avg_rsz_LFP_sh-sem_rsz_LFP_sh,avg_rsz_LFP_sh+sem_rsz_LFP_sh,'k');
    plotWinterval_AF_v0(rsz_bin_ctrs,avg_rsz_LFP_tm,avg_rsz_LFP_tm-sem_rsz_LFP_tm,avg_rsz_LFP_tm+sem_rsz_LFP_tm,[.5 .0 .5]);
    xticks(sort([phlim, 0]));    title('Average LFP');    plot(0*[1 1],ylim,'k--'); xlim('tight');              xlim(phlim);
    nexttile;   hold on;
    plotWinterval_AF_v0(rsz_bin_ctrs,avg_rsz_clk_sh,avg_rsz_clk_sh-sem_rsz_clk_sh,avg_rsz_clk_sh+sem_rsz_clk_sh,'k');
    plotWinterval_AF_v0(rsz_bin_ctrs,avg_rsz_clk_tm,avg_rsz_clk_tm-sem_rsz_clk_tm,avg_rsz_clk_tm+sem_rsz_clk_tm,[.5 .0 .5]);
    xticks(sort([phlim, 0]));    title('Average Call Rate');    plot(0*[1 1],ylim,'k--'); xlim('tight');        xlim(phlim);
    nexttile;   hold on;
    plotWinterval_AF_v0(rsz_bin_ctrs,avg_rsz_wbt_phase_sh,avg_rsz_wbt_phase_sh-sem_rsz_wbt_phase_sh,avg_rsz_wbt_phase_sh+sem_rsz_wbt_phase_sh,'k');
    plotWinterval_AF_v0(rsz_bin_ctrs,avg_rsz_wbt_phase_tm,avg_rsz_wbt_phase_tm-sem_rsz_wbt_phase_tm,avg_rsz_wbt_phase_tm+sem_rsz_wbt_phase_tm,[.5 .0 .5]);
    xticks(sort([phlim, 0]));    title('Wingbeat Phase');    plot(0*[1 1],ylim,'k--'); xlim('tight');           xlim(phlim);
    nexttile;   hold on;
    plotWinterval_AF_v0(rsz_bin_ctrs,avg_rsz_wbt_sh,avg_rsz_wbt_sh-sem_rsz_wbt_sh,avg_rsz_wbt_sh+sem_rsz_wbt_sh,'k');
    plotWinterval_AF_v0(rsz_bin_ctrs,avg_rsz_wbt_tm,avg_rsz_wbt_tm-sem_rsz_wbt_tm,avg_rsz_wbt_tm+sem_rsz_wbt_tm,[.5 .0 .5]);
    xticks(sort([phlim, 0]));    title('Accelerometer');    plot(0*[1 1],ylim,'k--'); xlim('tight');           xlim(phlim);
    sgtitle(['Average of ',num2str(N_sweeps_tm),' cycles, from ',num2str(numel(unique(SWP_sst_tm.flight_ID))),' flights']);
    
    %=== Show averages (Both)
    figure('units','normalized','outerposition',[.3 .4 .1 .3]);
    hold on;
    plotWinterval_AF_v0(rsz_bin_ctrs,avg_est_dist_sh,avg_est_dist_sh-sem_est_dist_sh,avg_est_dist_sh+sem_est_dist_sh,'k');
    plotWinterval_AF_v0(rsz_bin_ctrs,avg_est_dist_rl,avg_est_dist_rl-sem_est_dist_rl,avg_est_dist_rl+sem_est_dist_rl,[.9 .5 .2]);
    plotWinterval_AF_v0(rsz_bin_ctrs,avg_est_dist_tm,avg_est_dist_tm-sem_est_dist_tm,avg_est_dist_tm+sem_est_dist_tm,[.5 .0 .5]);
    xticks(sort([phlim, 0]));    title('Average Decoding Error');    plot(0*[1 1],ylim,'k--'); xlim('tight');   xlim(phlim);
    
    %=== Comparisons at increasing intervals around phase 0
    figure('units','normalized','outerposition',[0 .4 1 .3]);
    tiledlayout(1,18,'TileSpacing','compact');
    for i=1:18
        
        anglim = i*10.*[-1 1];
        baseline_sh = mean(mean(avg_est_dist_sh_all(rsz_bin_ctrs>anglim(1) & rsz_bin_ctrs<anglim(2),:),1));
        avg_dist_rl = (mean(est_dist_rl(rsz_bin_ctrs>anglim(1) & rsz_bin_ctrs<anglim(2),:),1)'-baseline_sh)./std(mean(avg_est_dist_sh_all(rsz_bin_ctrs>anglim(1) & rsz_bin_ctrs<anglim(2),:),1));
        avg_dist_tm = (mean(est_dist_tm(rsz_bin_ctrs>anglim(1) & rsz_bin_ctrs<anglim(2),:),1)'-baseline_sh)./std(mean(avg_est_dist_sh_all(rsz_bin_ctrs>anglim(1) & rsz_bin_ctrs<anglim(2),:),1));
        nexttile;
        if i~=9
            bar([mean(avg_dist_rl),mean(avg_dist_tm)]);
        else
            bar([mean(avg_dist_rl),mean(avg_dist_tm)],'FaceColor','r');
        end
        hold on;    errorbar(1:2, [mean(avg_dist_rl),mean(avg_dist_tm)], [std(avg_dist_rl)/sqrt(numel(avg_dist_rl)),std(avg_dist_tm)/sqrt(numel(avg_dist_tm))], 'k.', 'LineWidth', 1.5)
        xticks(1:2);    xticklabels({'Wingbeat', 'Non-oscillatory'});  ylabel('Average Decoding error (normalized)'); title(ranksum(avg_dist_rl,avg_dist_tm),anglim);
        signrank(avg_dist_rl)
        signrank(avg_dist_tm)
        
        
    end
    
end

%% LOOK AT AUTOMATICALLY CLASSIFIED SWEEPS FOR WINGBEAT AND LFP ANALYSIS
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
    SWPs_acc = {}; SWPs_sample = []; wbt_phase_all = []; tmr_phase_all = []; tmr_power_all = []; fract_pos_all = [];    SWPs_phs = [];  PSD_sweeps = [];    WFR_sweeps = [];    SWPs_clk = [];  SWPs_lfp = [];  
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
                SWPs_tmr = [SWPs_tmr; {SWPs_sst.tmr_phase{zz,1}}];
                SWPs_lfp = [SWPs_lfp; {SWPs_sst.LFP{zz,1}}];
                %SWPs_tmr = [SWPs_tmr; {bandpass(SWPs_sst.LFP{zz,1},[1 10],Fs_Hi)}];
                
                
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
    
    %=== Plot average of good vs. bad sweeps
    figure('units','normalized','outerposition',[.3 .3 .15 .3]);
    tiledlayout(1,2,'TileSpacing','compact');
    nexttile;   imagesc([0 size(SWPs_posterior,2)/Fs_Hi]-size(SWPs_posterior,2)/(2*Fs_Hi),[0  size(SWPs_posterior,1)*bin_size_1D],mean(SWPs_posterior(:,:, good_ids),3,'omitnan'),prctile(mean(SWPs_posterior(:,:, good_ids),3,'omitnan'),[1 99],'all')');    
    colormap(flipud(gray)); set(gca,'YDir','normal');   title('Valid');
    nexttile;   imagesc([0 size(SWPs_posterior,2)/Fs_Hi]-size(SWPs_posterior,2)/(2*Fs_Hi),[0  size(SWPs_posterior,1)*bin_size_1D],mean(SWPs_posterior(:,:, baad_ids),3,'omitnan'),prctile(mean(SWPs_posterior(:,:, baad_ids),3,'omitnan'),[1 99],'all')');    
    colormap(flipud(gray)); set(gca,'YDir','normal');   title('Invalid');
    colorbar
    
  
    
    
    %% ==================================================== ANALYSIS OF SWEEP VS ACCELEROMETER OR LFP =========================================================================
    
    id2keep = good_ids;
    S2W_st = table();                                               % Sweep to wingbeat table
    S2W_st.flightID = SWPs_flightID(id2keep);                       % Id of the flights
    S2W_st.acc = SWPs_acc(id2keep);                                 % Accelerometer
    S2W_st.phs = SWPs_phs(id2keep);                                 % Wingbeat Phase
    S2W_st.clk = SWPs_clk(id2keep);                                 % Echolocation
    S2W_st.tmr = SWPs_tmr(id2keep);                                 % Tamir Phase
    S2W_st.lfp = SWPs_lfp(id2keep);                                 % LFP
    S2W_st.sample = SWPs_sample(id2keep);                           % Sample of the sweep
    S2W_st.u_id = SWPs_u_id(id2keep,:);                             % Unique identifier
    [~,~,S2W_st.sessionID] = unique([string(S2W_st.u_id)],'rows');  % Identifier
    
    ss1 = 30; ss2= numel(S2W_st.acc{1,1})-ss1;  % 30 default
    
    %=== Keep only sweeps that start away from flight beginning/end
    S2W_st = S2W_st(S2W_st.sample>ss1 & S2W_st.sample<ss2,:);
    unique([string(S2W_st.u_id)],'rows');
    n_rep1 = 50;
    dt_interval = 0.06; % 0.04 default
    interval = [-round(Fs_Hi*dt_interval):round(Fs_Hi*dt_interval)];
    
    acc_accum_gd = zeros(size(S2W_st,1),numel(interval));
    acc_accum_sh_gd = zeros(size(S2W_st,1),numel(interval),n_rep1);
    lfp_accum_gd = zeros(size(S2W_st,1),numel(interval));
    lfp_accum_sh_gd = zeros(size(S2W_st,1),numel(interval),n_rep1);
    wbt_accum_gd = zeros(size(S2W_st,1),numel(interval));
    tmr_accum_gd = zeros(size(S2W_st,1),numel(interval));
    for i=1:size(S2W_st,1)
        
        acc_accum_gd(i,:) = S2W_st.acc{i,1}(S2W_st.sample(i)+interval,1)';
        lfp_accum_gd(i,:) = S2W_st.lfp{i,1}(S2W_st.sample(i)+interval,1)';
        
        wbt_accum_gd(i,:) = S2W_st.phs{i,1}(S2W_st.sample(i)+interval,1)';
        tmr_accum_gd(i,:) = S2W_st.tmr{i,1}(S2W_st.sample(i)+interval,1)';
        
        for jj =1:n_rep1
            sample_shuffle = 100+randi(200);
            acc_accum_sh_gd(i,:,jj) = S2W_st.acc{i,1}(sample_shuffle+interval,1)';
            lfp_accum_sh_gd(i,:,jj) = S2W_st.lfp{i,1}(sample_shuffle+interval,1)';
        end
        
    end
    acc_accum_sh_mean_gd = squeeze(mean(acc_accum_sh_gd,3));
    tmr_accum_sh_mean_gd = squeeze(mean(lfp_accum_sh_gd,3));
    
    figure('units','normalized','outerposition',[.1 .1 .2 .3]);
    tiledlayout(1,2,'TileSpacing','tight');
    nexttile;
    plotWinterval_AF_v0(interval./Fs_Hi,mean(acc_accum_gd),mean(acc_accum_gd)-std(acc_accum_gd)./sqrt(size(S2W_st,1)),mean(acc_accum_gd)+std(acc_accum_gd)./sqrt(size(S2W_st,1)),'r');    hold on;
    plotWinterval_AF_v0(interval./Fs_Hi,mean(acc_accum_sh_mean_gd),mean(acc_accum_sh_mean_gd)-std(acc_accum_sh_mean_gd)./sqrt(size(S2W_st,1)),mean(acc_accum_sh_mean_gd)+std(acc_accum_sh_mean_gd)./sqrt(size(S2W_st,1)),'k');
    xlabel('Time from Sweep Peak (s)'); ylabel('Average Acceleration (g)');         ylim([-0.05 0.05]);
    nexttile;
    plotWinterval_AF_v0(interval./Fs_Hi,mean(lfp_accum_gd),mean(lfp_accum_gd)-std(lfp_accum_gd)./sqrt(size(S2W_st,1)),mean(lfp_accum_gd)+std(lfp_accum_gd)./sqrt(size(S2W_st,1)),'r');    hold on;
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
    S2W_st.lfp = SWPs_lfp(id2keep);                           % LFP
    S2W_st.sample = SWPs_sample(id2keep);                     % Sample of the sweep
    S2W_st.u_id = SWPs_u_id(id2keep,:);                       % Unique identifier
    [~,~,S2W_st.sessionID] = unique([string(S2W_st.u_id)],'rows');% Identifier
    
    %=== Keep only sweeps that start away from flight beginning/end
    S2W_st = S2W_st(S2W_st.sample>ss1 & S2W_st.sample<ss2,:);
    unique([string(S2W_st.u_id)],'rows');
    
    acc_accum_bd = zeros(size(S2W_st,1),numel(interval));
    acc_accum_sh_bd = zeros(size(S2W_st,1),numel(interval),n_rep1);
    lfp_accum_bd = zeros(size(S2W_st,1),numel(interval));
    lfp_accum_sh_bd = zeros(size(S2W_st,1),numel(interval),n_rep1);
    wbt_accum_bd = zeros(size(S2W_st,1),numel(interval));
    tmr_accum_bd = zeros(size(S2W_st,1),numel(interval));
    for i=1:size(S2W_st,1)
        
        acc_accum_bd(i,:) = S2W_st.acc{i,1}(S2W_st.sample(i)+interval,1)';
        lfp_accum_bd(i,:) = S2W_st.lfp{i,1}(S2W_st.sample(i)+interval,1)';
        
        wbt_accum_bd(i,:) = S2W_st.phs{i,1}(S2W_st.sample(i)+interval,1)';
        tmr_accum_bd(i,:) = S2W_st.tmr{i,1}(S2W_st.sample(i)+interval,1)';
        
        for jj =1:n_rep1
            sample_shuffle = 100+randi(200);
            acc_accum_sh_bd(i,:,jj) = S2W_st.acc{i,1}(sample_shuffle+interval,1)';
            lfp_accum_sh_bd(i,:,jj) = S2W_st.lfp{i,1}(sample_shuffle+interval,1)';
        end
        
    end
    acc_accum_sh_mean_bd = squeeze(mean(acc_accum_sh_bd,3));
    tmr_accum_sh_mean_bd = squeeze(mean(lfp_accum_sh_bd,3));
    
    figure('units','normalized','outerposition',[.4 .1 .2 .3]);
    tiledlayout(1,2,'TileSpacing','tight');
    nexttile;
    plotWinterval_AF_v0(interval./Fs_Hi,mean(acc_accum_bd),mean(acc_accum_bd)-std(acc_accum_bd)./sqrt(size(S2W_st,1)),mean(acc_accum_bd)+std(acc_accum_bd)./sqrt(size(S2W_st,1)),'b');    hold on;
    plotWinterval_AF_v0(interval./Fs_Hi,mean(acc_accum_sh_mean_bd),mean(acc_accum_sh_mean_bd)-std(acc_accum_sh_mean_bd)./sqrt(size(S2W_st,1)),mean(acc_accum_sh_mean_bd)+std(acc_accum_sh_mean_bd)./sqrt(size(S2W_st,1)),'k');
    xlabel('Time from Sweep Peak (s)'); ylabel('Average Acceleration (g)');         ylim([-0.05 0.05]);
    nexttile;
    plotWinterval_AF_v0(interval./Fs_Hi,mean(lfp_accum_bd),mean(lfp_accum_bd)-std(lfp_accum_bd)./sqrt(size(S2W_st,1)),mean(lfp_accum_bd)+std(lfp_accum_bd)./sqrt(size(S2W_st,1)),'b');    hold on;
    plotWinterval_AF_v0(interval./Fs_Hi,mean(tmr_accum_sh_mean_bd),mean(tmr_accum_sh_mean_bd)-std(tmr_accum_sh_mean_bd)./sqrt(size(S2W_st,1)),mean(tmr_accum_sh_mean_bd)+std(tmr_accum_sh_mean_bd)./sqrt(size(S2W_st,1)),'k');
    xlabel('Time from Sweep Peak (s)'); ylabel('Average LFP (uV)');                 ylim([-10 3]);
    sgtitle('Bad Sweeps');
    
    %=== Good vs bad
    figure('units','normalized','outerposition',[.4 .1 .2 .3]);
    tiledlayout(1,2,'TileSpacing','tight');
    nexttile;
    plotWinterval_AF_v0(interval./Fs_Hi,mean(acc_accum_bd),mean(acc_accum_bd)-std(acc_accum_bd)./sqrt(size(S2W_st,1)),mean(acc_accum_bd)+std(acc_accum_bd)./sqrt(size(S2W_st,1)),'b');    hold on;
    plotWinterval_AF_v0(interval./Fs_Hi,mean(acc_accum_gd),mean(acc_accum_gd)-std(acc_accum_gd)./sqrt(size(S2W_st,1)),mean(acc_accum_gd)+std(acc_accum_gd)./sqrt(size(S2W_st,1)),'r');
    xlabel('Time from Sweep Peak (s)'); ylabel('Average Acceleration (g)');
    nexttile;
    plotWinterval_AF_v0(interval./Fs_Hi,mean(lfp_accum_bd),mean(lfp_accum_bd)-std(lfp_accum_bd)./sqrt(size(S2W_st,1)),mean(lfp_accum_bd)+std(lfp_accum_bd)./sqrt(size(S2W_st,1)),'b');    hold on;
    plotWinterval_AF_v0(interval./Fs_Hi,mean(lfp_accum_gd),mean(lfp_accum_gd)-std(lfp_accum_gd)./sqrt(size(S2W_st,1)),mean(lfp_accum_gd)+std(lfp_accum_gd)./sqrt(size(S2W_st,1)),'r');
    xlabel('Time from Sweep Peak (s)'); ylabel('Average LFP (uV)');
    sgtitle('Good vs Bad Sweeps');
    
    trace2probe = mean(acc_accum_gd);
    p_val_trace2probe = zeros(size(trace2probe));
    for j=1:numel(trace2probe)
        p_val_trace2probe(j) = sum(trace2probe(j)> acc_accum_bd(:,j))./size(acc_accum_bd,1);
    end
    
    %=== Plot using confidence intervals
    figure('units','normalized','outerposition',[.4 .4 .2 .3]);
    tiledlayout(1,2,'TileSpacing','tight');
    nexttile;
    ci = bootci(100, @mean, acc_accum_bd); plotWinterval_AF_v0(interval./Fs_Hi,mean(acc_accum_bd),ci(2,:), ci(1,:),'b');    hold on;
    ci = bootci(100, @mean, acc_accum_gd); plotWinterval_AF_v0(interval./Fs_Hi,mean(acc_accum_gd),ci(2,:), ci(1,:),'r');
    xlabel('Time from Sweep Peak (s)'); ylabel('Average Acceleration (g)');
    nexttile;
    ci = bootci(100, @mean, lfp_accum_bd); plotWinterval_AF_v0(interval./Fs_Hi,mean(lfp_accum_bd),ci(2,:), ci(1,:),'b');    hold on;
    ci = bootci(100, @mean, lfp_accum_gd); plotWinterval_AF_v0(interval./Fs_Hi,mean(lfp_accum_gd),ci(2,:), ci(1,:),'r');
    xlabel('Time from Sweep Peak (s)'); ylabel('Average LFP (uV)');
    sgtitle('Good vs Bad Sweeps (bootstrapping) CIs');
    
    tmp_idx = round(numel(interval)/2);
    
    circ_wwtest(tmr_accum_gd(:,tmp_idx), wbt_accum_gd(:,tmp_idx))
    
end

