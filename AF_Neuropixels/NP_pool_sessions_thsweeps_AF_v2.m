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
    
    %=== Normalize on the whole spectrum and scale between 0 and 1
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
    plot_distr_AF_v0(pwr_b, pwr_f, {'Bouts', 'Flight'}, 'SEM', 'Relative Theta Power');
    
    signrank(pwr_b, pwr_f)
    
    disp(['Fraction bouts during rest: ', num2str(1-mean([LFP.fract_bouts_flight]),3)]);
    
    std(1-[LFP.fract_bouts_flight])./sqrt(numel([LFP.fract_bouts_flight]))
    
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
        if numel(NP_units.spk_phase2wbt_F{nc,1})>10
            
            %=== Fit with cosine function
            x = [phase_ctrs,phase_ctrs+2*pi];
            y = repmat(smoothdata(NP_units.spk_phase2wbt_F_counts{nc,:}','movmean',3),1,2);
            x1 = x(~isnan(y));  y1 = y(~isnan(y));
            y1 = normalize(y1-mean(y1),'range',[-1 1]);
            [f1,gof1] = fit(x1',y1',cosEqn,'Start',[2 0],'Lower',[-inf 0]);
            [f2,gof2] = fit(x1',y1',cosEqn,'Start',[1 0],'Lower',[-inf 0]);
            
            %=== Keep best fit
            if gof1.rsquare>gof2.rsquare
                NP_units.spk_phase2wbt_R2(nc) = gof1.rsquare;
                NP_units.spk_phase2wbt_fit_par(nc,:) = coeffvalues(f1);
            else
                NP_units.spk_phase2wbt_R2(nc) = gof2.rsquare;
                NP_units.spk_phase2wbt_fit_par(nc,:) = coeffvalues(f2);
            end
            NP_units.spk2wbt_p(nc) = circ_rtest(NP_units.spk_phase2wbt_F{nc,1});    % p value Rayleight test
            [~,max_loc] = max(smoothdata(NP_units.spk_phase2wbt_F_counts{nc,:},'movmean',3));
            NP_units.spk2wbt_pref_phase(nc) = phase_ctrs(max_loc);                  % Preferred phase
            
            %                     %=== Plotting (troubleshoot)
            %                     if gof1.rsquare>gof2.rsquare
            %                         plot(f1);    hold on;  bar(x1,y1,1,'EdgeColor','none'); hold off;   xticks()
            %                     else
            %                         plot(f2);    hold on;  bar(x1,y1,1,'EdgeColor','none'); hold off;
            %                     end
            %                       title(rad2deg(phase_ctrs(max_loc)));
            %                     choice = questdlg('Continue', 'Class Verification','Yes', 'No', 'Stop','Yes');
            %                     switch choice
            %                         case 'Stop'
            %                             break;
            %                     end
        else
            NP_units.spk_phase2wbt_R2(nc) = NaN;
            NP_units.spk_phase2wbt_fit_par(nc,:) = NaN(1,2);
            NP_units.spk2wbt_p(nc) = NaN;
            NP_units.spk2wbt_pref_phase(nc) = NaN;
        end
        
    end
    
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

%% LOOK AT AUTOCORRELOGRAMS AND PHASE PREFERENCES
for hide=1
    
    %=== Define subset of cells of interest
    %NP_sst = NP_units(NP_units.fr>0 & NP_units.nspkpf>0 & NP_units.place_cond & NP_units.spk2wbt_p<0.05,:);
    NP_sst = NP_units(NP_units.fr>0 & NP_units.nspkpf>0 & NP_units.place_cond & NP_units.num_spikes_f>50,:);
    N_cells = size(NP_sst,1); % Number of cells
    
    %=== Initialize relevant variables
    AC = zeros(numel(NP_sst.AC{1,1}),N_cells);
    AC_ctrs = NP_sst.AC_bins{1,1};
    R2 = zeros(N_cells,1);
    cos_par = zeros(N_cells,4);
    spk2wbt_counts = zeros(numel(NP_sst.spk_phase2wbt_F_counts{1,:}),N_cells);
    all2wbt_counts = zeros(numel(NP_sst.spk_phase2wbt_F_counts{1,:}),N_cells);
    spk2wbt_nrmcts = zeros(numel(NP_sst.spk_phase2wbt_F_counts{1,:}),N_cells);
    phase_bins = NP_sst.phase_bins(1,:);
    phase_ctrs = NP_sst.phase_bins(1,1:end-1)+mean(diff(NP_sst.phase_bins(1,:)))/2;
    
    %=== Accumulate data across cells
    for nc=1:N_cells
        
        AC(:,nc) = NP_sst.AC{nc,1};
        spk2wbt_counts(:,nc) = NP_sst.spk_phase2wbt_F_counts{nc,:};
        all2wbt_counts(:,nc) = NP_sst.all_phase2wbt_counts{nc,:};
        spk2wbt_nrmcts(:,nc) = NP_sst.spk_phase2wbt_F_counts{nc,:}-NP_sst.all_phase2wbt_counts{nc,:};
        
    end
    
    %=== Impose specific conditions
    %cond = all([NP_sst.spk_phase2wbt_R2>0.25, NP_sst.spk_phase2wbt_fit_par(:,1)>0],2);
    cond = all([NP_sst.spk_phase2wbt_R2>0, NP_sst.num_spikes_f>50],2);
    
    %=== Fit average autocorrelogram
    expEqn = 'a*exp(-b*x+c)+d';
    x = AC_ctrs(AC_ctrs>0);
    y1 = mean(AC(AC_ctrs>0,:   ),2,'omitnan');
    y2 = mean(AC(AC_ctrs>0,cond),2,'omitnan');
    f1 = fit(x,y1,expEqn,'Start',[max(y1) 1 0 0],'Lower',[0 -inf 0 0],'Exclude',x<0.125);
    f2 = fit(x,y2,expEqn,'Start',[max(y2) 1 0 0],'Lower',[0 -inf 0 0],'Exclude',x<0.125);
    [psd_f1,pxx] = pwelch(y1-f1(x),[],[],[],1/mean(diff(x)));
    [psd_f2,~] = pwelch(y2-f2(x),[],[],[],1/mean(diff(x)));
    [~,max_loc] = max(psd_f2);
    
    %=== Plot average autocorrelogram
    figure('units','normalized','outerposition',[0.4 .5 0.25 .5]);
    tiledlayout(2,3,'TileSpacing','tight');
    nexttile;   area(AC_ctrs,mean(AC,2,'omitnan'),'FaceColor','k','EdgeColor','none','FaceAlpha',0.5);
    yticks([]); xlabel('Time lag (s)'); hold on;
    plot(repmat(1/wb_freq*[-3:3],2,1),repmat(ylim,7,1)','k--');
    xlim([0 max(AC_ctrs)]); title('All Place Cells');
    nexttile;   plot(x,smoothdata(y1-f1(x),'movmean',3),'LineWidth',3,'Color','k');   xlim([0.05 max(AC_ctrs)]);    hold on;
    plot(repmat(1/wb_freq*[-3:3],2,1),repmat(ylim,7,1)','k--'); title('Residual');  yticks([]); xlabel('Time lag (s)');
    sgtitle([num2str(sum(cond)),' neurons']);
    nexttile;   plot(pxx,10*log10(psd_f1),'LineWidth',3,'Color','k'); hold on; plot(wb_freq*[1 1],ylim,'k--');    xlabel('Freq');
    nexttile;   area(AC_ctrs,mean(AC(:,cond),2,'omitnan'),'FaceColor','b','EdgeColor','none','FaceAlpha',0.5);
    yticks([]); xlabel('Time lag (s)'); hold on;    plot(x,f2(x));
    plot(repmat(1/wb_freq*[-3:3],2,1),repmat(ylim,7,1)','k--');
    xlim([0 max(AC_ctrs)]); title('Good Cells');
    nexttile;   plot(x,smoothdata(y2-f2(x),'movmean',3),'LineWidth',3,'Color','b');   xlim([0.05 max(AC_ctrs)]);
    yticks([]); xlabel('Time lag (s)'); hold on;
    plot(repmat(1/wb_freq*[-3:3],2,1),repmat(ylim,7,1)','k--'); title('Residual');
    nexttile;   plot(pxx,10*log10(psd_f2),'LineWidth',3,'Color','b'); hold on; plot(wb_freq*[1 1],ylim,'k--');    xlabel('Freq');
    sgtitle([num2str(sum(cond)),'/',num2str(N_cells),' neurons']);
    
    disp(['Peak frequency autocorrelogram: ', num2str(pxx(max_loc),2), ' Hz']);
    
    %=== Plot average spike to wingbeat phase
    figure('units','normalized','outerposition',[0.05 .1 0.25 .6]);
    tiledlayout(3,3,'TileSpacing','tight');
    phase2plot = repmat(median(spk2wbt_counts,2,'omitnan')',1,2); delta_h = (max(phase2plot)-mean(phase2plot));
    nexttile;   bar([phase_ctrs,phase_ctrs+2*pi],smoothdata(phase2plot,'movmean',2),1,'k','edgecolor','none');  ylim(mean(phase2plot)+delta_h*[-2 2]);
    xticks(-pi:pi:3*pi);    xticklabels({'-180', '0', '180', '360', '540'});   xlabel('Wingbeat');  title('All Place Cells');
    phase2plot = repmat(median(spk2wbt_counts(:,cond),2,'omitnan')',1,2); delta_h = (max(phase2plot)-mean(phase2plot));
    nexttile;   bar([phase_ctrs,phase_ctrs+2*pi],smoothdata(phase2plot,'movmean',2),1,'b','edgecolor','none');  ylim(mean(phase2plot)+delta_h*[-2 2]);
    xticks(-pi:pi:3*pi);    xticklabels({'-180', '0', '180', '360', '540'});   xlabel('Wingbeat');  title('Good Cells');
    phase2plot = repmat(median(all2wbt_counts(:,cond),2,'omitnan')',1,2); delta_h = (max(phase2plot)-mean(phase2plot));
    nexttile;   bar([phase_ctrs,phase_ctrs+2*pi],smoothdata(phase2plot,'movmean',2),1,'g','edgecolor','none');  ylim(mean(phase2plot)+delta_h*[-2 2]);
    xticks(-pi:pi:3*pi);    xticklabels({'-180', '0', '180', '360', '540'});   xlabel('Wingbeat');  title('Random Spikes');
    phase2plot = repmat(median(spk2wbt_nrmcts,2,'omitnan')',1,2); delta_h = (max(phase2plot)-mean(phase2plot));
    nexttile;   bar([phase_ctrs,phase_ctrs+2*pi],smoothdata(phase2plot,'movmean',2),1,'k','edgecolor','none');  ylim(mean(phase2plot)+delta_h*[-2 2]);
    xticks(-pi:pi:3*pi);    xticklabels({'-180', '0', '180', '360', '540'});   xlabel('Wingbeat');  title('Real - Random');
    phase2plot = repmat(median(spk2wbt_nrmcts(:,cond),2,'omitnan')',1,2); delta_h = (max(phase2plot)-mean(phase2plot));
    nexttile;   bar([phase_ctrs,phase_ctrs+2*pi],smoothdata(phase2plot,'movmean',2),1,'b','edgecolor','none');  ylim(mean(phase2plot)+delta_h*[-2 2]);
    xticks(-pi:pi:3*pi);    xticklabels({'-180', '0', '180', '360', '540'});   xlabel('Wingbeat');  title('Real - Random');
    phase2plot = repmat(median(all2wbt_counts(:,cond),2,'omitnan')',1,2); delta_h = (max(phase2plot)-mean(phase2plot));
    nexttile;   polarhistogram(NP_sst.spk2wbt_pref_phase,linspace(-pi,pi,20),'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','k');              title('Preferred Phase (all)');
    nexttile;   polarhistogram(NP_sst.spk_phase2wbt_fit_par(:,2),linspace(-pi,pi,20),'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','k');      title('Fit param (all)');
    nexttile;   polarhistogram(NP_sst.spk_phase2wbt_fit_par(cond,2),linspace(-pi,pi,20),'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','b');   title('Fit param (good)');
    nexttile;   polarhistogram(NP_sst.spk2wbt_pref_phase(cond),linspace(-pi,pi,20),'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','r');        title('Preferred Phase (good)');
    sgtitle([num2str(sum(cond)),'/',num2str(N_cells),' neurons']);
    
    %=== Display some quantifications
    preferred_phase = NP_sst.spk2wbt_pref_phase(cond);
    average_phase_cts = mean(spk2wbt_counts(:,cond),2,'omitnan')';
    [~,max_loc] = max(average_phase_cts);
    [circ_s, circ_s0] = circ_std(preferred_phase);
    disp(['Fraction of Cells: ', num2str(sum(cond)/N_cells)]);
    disp(['Median phase of max firing: ', num2str(rad2deg(circ_median(preferred_phase)),2),' p/m ', num2str(rad2deg(circ_s0)./sqrt(numel(preferred_phase)),2)    , ' degrees']);
    disp(['Peak of mean firing phase: ', num2str(rad2deg(phase_ctrs(max_loc)),2),' degrees']);
    
    %=== Fit with cosine
    x = [phase_ctrs,phase_ctrs+2*pi];
    y = repmat(median(spk2wbt_counts(:,cond),2,'omitnan')',1,2);
    x1 = x(~isnan(y));  y1 = y(~isnan(y));
    y1 = normalize(y1-mean(y1),'range',[-1 1]);
    [f1,gof1] = fit(x1',y1',cosEqn,'Start',[1 0],'Lower',[-inf 0]);
    fit_coeff = coeffvalues(f1);
    fit_coeff_int = confint(f1);
    
    figure;
    plot(f1);    hold on;  bar(x1,y1,1,'EdgeColor','none','FaceColor','k'); legend('off');
    xticks(-pi:pi:3*pi);    xticklabels({'-180', '0', '180', '360', '540'});   
    title(['Max at: ', num2str(rad2deg(fit_coeff(2)),3), '( ', num2str( rad2deg(fit_coeff_int(2,1)),3),'-',num2str( rad2deg(fit_coeff_int(2,2)),3),')'])
    
    %=== Plot phase locking for some example cells
    rng(1);
    subset = find(NP_sst.spk_phase2wbt_R2>0 & NP_sst.spk_phase2wbt_fit_par(:,1)<1.5);
    subset = subset(randperm(numel(subset)));
    %subset = find(NP_sst.spk2wbt_p<0.05);
    figure('units','normalized','outerposition',[.5 0 .5 1]);
    tiledlayout(10,10,'TileSpacing','tight');
    for i=1:min(100,numel(subset))
        nexttile;   bar([phase_ctrs,phase_ctrs+2*pi],repmat(smoothdata(NP_sst.spk_phase2wbt_F_counts{subset(i),:}','movmean',3),1,2),1,'k','edgecolor','none');
        xticks(-pi:pi:3*pi);    xticklabels({'-180', '0', '180', '360', '540'});
    end
    title('Example Neurons');
    
    %=== Phase procession
    NP_phlock = NP_sst(cond,:);
    sum(NP_phlock.phase_prec)/numel(NP_phlock.phase_prec);
    figure('units','normalized','outerposition',[.5 .2 .17 .3])
    histogram(NP_phlock.phs_prc_asym,[0:0.1:1]);   xlabel('Asymmetry level');
    
end

%% THETA-SWEEPS AND WINGBEAT
for hide=1
    
    %=== Select subtable of good flights
    SWPs_sst = SWPs(SWPs.rmsDec_error<2 & SWPs.prc_decoded>0.5,:);
    
    %=== Params and initialize relevant variables
    n_reshape = 50;                 % Default 30
    smooth_f = [1 .3];              % [space bin,time bin]
    N_flights = size(SWPs_sst,1);   % Number of flights
    single_SWP = table();   counter=1;           
    warning('off');
    %=== Cut at wingbeat maxima
    for zz=1:N_flights
        
        zero_phs_idx = find(SWPs_sst.wbt_phase{zz,1}(1:end-1).* SWPs_sst.wbt_phase{zz,1}(2:end)<0 & diff(SWPs_sst.wbt_phase{zz,1})<0);  % Segment based on wingbeat
        %zero_phs_idx = find(SWPs_sst.tmr_phase{zz,1}(1:end-1).* SWPs_sst.tmr_phase{zz,1}(2:end)<0 & diff(SWPs_sst.tmr_phase{zz,1})>0); % Segment based on Tamir's phase
        %plot(SWPs_sst.wbt_phase{zz,1});   hold on; stem(zero_phs_idx, ones(size(zero_phs_idx )));  plot(SWPs_sst.wbt{zz,1});
        sweep_strt = zero_phs_idx(1:end-1); sweep_stop = zero_phs_idx(2:end);
        N_spatial_bins = size(SWPs_sst.p_dec_flight{zz,1},1);
        spt_bin_ids = [1:N_spatial_bins]';
        
        for ss=1:numel(sweep_strt)
            
            %=== Basic features of the sweep
            single_SWP.flight_ID(counter) = zz;
            single_SWP.rms_dec(counter) = SWPs_sst.rmsDec_error(zz);
            single_SWP.prc_dec(counter) = SWPs_sst.prc_decoded(zz);
            
            %=== Raw, filtered, shifted and rehaped posterior
            single_SWP.raw_posterior(counter) = {SWPs_sst.p_dec_flight{zz,1}(:,sweep_strt(ss):sweep_stop(ss))};
            single_SWP.reg_posterior(counter) = {imgaussfilt(SWPs_sst.p_dec_flight{zz,1}(:,sweep_strt(ss):sweep_stop(ss)),smooth_f)};
            single_SWP.sft_posterior(counter) = {imgaussfilt(SWPs_sst.p_dec_shifted{zz,1}(:,sweep_strt(ss):sweep_stop(ss)),smooth_f)};
            single_SWP.rsp_posterior(counter) = {imresize(single_SWP.sft_posterior{counter,1},[size(single_SWP.sft_posterior{counter,1},1),n_reshape])};
            
            valid_bins = sum(SWPs_sst.p_dec_flight{zz,1}(:,sweep_strt(ss):sweep_stop(ss)))>0.8;
            
            %=== Spike density, wingbeat, LFP, Tamir's phase echolocation and phase
            single_SWP.spk_dsty(counter) = {SWPs_sst.spk_dsty{zz,1}(sweep_strt(ss):sweep_stop(ss))};
            single_SWP.rsz_spk_dsty(counter) = {interp1(single_SWP.spk_dsty{counter,1},linspace(1,numel(single_SWP.spk_dsty{counter,1}),n_reshape)')};
            single_SWP.wbt_power(counter) = {SWPs_sst.wbt_power{zz,1}(sweep_strt(ss):sweep_stop(ss))};
            single_SWP.LFP(counter) = {SWPs_sst.LFP{zz,1}(sweep_strt(ss):sweep_stop(ss))};
            single_SWP.rsz_LFP(counter) = {interp1(single_SWP.LFP{counter,1},linspace(1,numel(single_SWP.LFP{counter,1}),n_reshape)')};
            single_SWP.wbt(counter) = {SWPs_sst.wbt{zz,1}(sweep_strt(ss):sweep_stop(ss))};
            single_SWP.rsz_wbt(counter) = {interp1(single_SWP.wbt{counter,1},linspace(1,numel(single_SWP.wbt{counter,1}),n_reshape)')};
            single_SWP.fract_pos(counter) = {SWPs_sst.pos_real{zz,1}(sweep_strt(ss):sweep_stop(ss))/SWPs_sst.pos_real{zz,1}(end)};
            single_SWP.clk(counter) = {SWPs_sst.clk{zz,1}(sweep_strt(ss):sweep_stop(ss))};
            single_SWP.rsz_clk(counter) = {interp1(single_SWP.clk{counter,1},linspace(1,numel(single_SWP.clk{counter,1}),n_reshape)')};
            single_SWP.wbt_phase(counter) = {SWPs_sst.wbt_phase{zz,1}(sweep_strt(ss):sweep_stop(ss))};
            single_SWP.rsz_wbt_phase(counter) = {interp1(single_SWP.wbt_phase{counter,1},linspace(1,numel(single_SWP.wbt_phase{counter,1}),n_reshape)')};
            single_SWP.tmr_phase(counter) = {SWPs_sst.tmr_phase{zz,1}(sweep_strt(ss):sweep_stop(ss))};
            single_SWP.rsz_tmr_phase(counter) = {interp1(single_SWP.tmr_phase{counter,1},linspace(1,numel(single_SWP.tmr_phase{counter,1}),n_reshape)')};
            
            %=== Add some features of the single cycle
            single_SWP.mean_spk_dsty(counter) = mean(single_SWP.spk_dsty{counter,1});   % Average spike density
            single_SWP.mean_wbt_power(counter) = mean(single_SWP.wbt_power{counter,1}); % Average wingbeat power
            single_SWP.mean_fract_pos(counter) = mean(single_SWP.fract_pos{counter,1}); % Average phase of the flight
            
            %=== Add some features of the posterior
            [max_p,max_loc] = max(single_SWP.sft_posterior{counter,1},[],1);        % Location of max posterior
            cnt_mass = spt_bin_ids'*single_SWP.sft_posterior{counter,1};            % Center of mass
            single_SWP.med_jmp_distance(counter) = median(abs(diff(cnt_mass)));     % Median jump distance with center of mass
            single_SWP.mean_jmp_distance(counter) = mean(abs(diff(max_loc)));       % Mean jump distance with loc of max posterior
%             pos_sprd = 0;
%             for bb=1:numel(cnt_mass)
%                 pos_sprd = pos_sprd+sqrt((spt_bin_ids'-cnt_mass(bb)).^2*single_SWP.sft_posterior{counter,1}(:,bb));
%             end
%             single_SWP.avg_pos_spread(counter) = pos_sprd/numel(cnt_mass);          % Average positional spread
            single_SWP.med_max_post(counter) = mean(max_p);                         % Average max posterior
            single_SWP.max_loc(counter) = {max_loc-N_spatial_bins/2};               % Decoding error (max posterior)
            single_SWP.cnt_mas(counter) = {cnt_mass-N_spatial_bins/2};              % Decoding error (center of mass)
            single_SWP.est_dist1(counter) = {interp1((cnt_mass-N_spatial_bins/2)*bin_size_1D,linspace(1,numel(cnt_mass),n_reshape)')};  % Decoding error, reshaped and in m (center of mass)
            single_SWP.est_dist2(counter) = {interp1((max_loc-N_spatial_bins/2)*bin_size_1D,linspace(1,numel(cnt_mass),n_reshape)')};   % Decoding error, reshaped and in m (max posterior)
            
            counter=counter+1;
        end
        
    end
    warning('on');
    
    %%
    %=== Extract subset using defined criteria (exclude flight tails, epochs of low firing and flat sweeps)
    %SWP_sst = single_SWP(single_SWP.mean_fract_pos>0.15 & single_SWP.mean_fract_pos<0.85 & single_SWP.mean_spk_dsty>prctile(single_SWP.mean_spk_dsty,10) & single_SWP.med_jmp_distance>0.0 & single_SWP.prc_dec>0.7 & single_SWP.rms_dec<1.5,:);
    SWP_sst = single_SWP(single_SWP.mean_fract_pos>0.15 & single_SWP.mean_fract_pos<0.85 & single_SWP.mean_spk_dsty>prctile(single_SWP.mean_spk_dsty,10) & single_SWP.mean_jmp_distance>0.1 & single_SWP.prc_dec>0.5 & single_SWP.rms_dec<1.5 & single_SWP.med_max_post>0.1,:);
    N_sweeps = size(SWP_sst,1);
    
    %=== Calculate averages and STDs
    est_dist = zeros(n_reshape,N_sweeps);
    rsz_spk_dsty = zeros(n_reshape,N_sweeps);
    rsz_wbt = zeros(n_reshape,N_sweeps);
    rsz_LFP = zeros(n_reshape,N_sweeps);
    rsz_clk = zeros(n_reshape,N_sweeps);
    rsz_wbt_phase = zeros(n_reshape,N_sweeps);
    for i=1:N_sweeps
        est_dist(:,i) = SWP_sst.est_dist2{i,1};
        rsz_spk_dsty(:,i) = SWP_sst.rsz_spk_dsty{i,1};
        rsz_wbt(:,i) = SWP_sst.rsz_wbt{i,1};
        rsz_LFP(:,i) = SWP_sst.rsz_LFP{i,1};
        rsz_clk(:,i) = SWP_sst.rsz_clk{i,1};
        rsz_wbt_phase(:,i) = SWP_sst.rsz_wbt_phase{i,1};
    end
    avg_est_dist = mean(est_dist,2);            sem_est_dist = std(est_dist,[],2)/sqrt(N_sweeps);
    avg_rsz_spk_dsty = mean(rsz_spk_dsty,2);    sem_rsz_spk_dsty = std(rsz_spk_dsty,[],2)/sqrt(N_sweeps);
    avg_rsz_wbt = mean(rsz_wbt,2);              sem_rsz_wbt = std(rsz_wbt,[],2)/sqrt(N_sweeps);
    avg_rsz_LFP = mean(rsz_LFP,2);              sem_rsz_LFP = std(rsz_LFP,[],2)/sqrt(N_sweeps);
    avg_rsz_wbt_phase = mean(rsz_wbt_phase,2);  sem_rsz_wbt_phase = std(rsz_wbt_phase,[],2)/sqrt(N_sweeps);
    avg_rsz_clk = mean(rsz_clk,2);              sem_rsz_clk = std(rsz_clk,[],2)/sqrt(N_sweeps);
    
    %=== Define bins and find peak decoding error phase 
    rsz_bin_ctrs = linspace(-180,180,n_reshape);
    [max_val,max_loc] = max(avg_est_dist);                rsz_bin_ctrs(max_loc);
    sum(abs(diff(SWP_sst.flight_ID)))+1;
    
    %=== Show averages
    figure('units','normalized','outerposition',[.3 .3 .45 .6]);
    tiledlayout(2,5,'TileSpacing','compact');
    nexttile;
    plotWinterval_AF_v0(rsz_bin_ctrs,avg_est_dist,avg_est_dist-sem_est_dist,avg_est_dist+sem_est_dist,'k');
    hold on;    plot(0*[1 1],ylim,'k--'); xlim('tight');
    xticks([-180 0 180]);    title('Average Decoding Error');    %yticks([]);
    nexttile;
    plotWinterval_AF_v0(rsz_bin_ctrs,avg_rsz_spk_dsty,avg_rsz_spk_dsty-sem_rsz_spk_dsty,avg_rsz_spk_dsty+sem_rsz_spk_dsty,'r');
    hold on;    plot(0*[1 1],ylim,'k--'); xlim('tight');
    xticks([-180 0 180]);    title('Average Spike Density');    %yticks([]);
    nexttile;
    plotWinterval_AF_v0(rsz_bin_ctrs,avg_rsz_LFP,avg_rsz_LFP-sem_rsz_LFP,avg_rsz_LFP+sem_rsz_LFP,'g');
    hold on;    plot(0*[1 1],ylim,'k--'); xlim('tight');
    xticks([-180 0 180]);    title('Average LFP');    %yticks([]);
    nexttile;
    plotWinterval_AF_v0(rsz_bin_ctrs,avg_rsz_clk,avg_rsz_clk-sem_rsz_clk,avg_rsz_clk+sem_rsz_clk,'b');
    hold on;    plot(0*[1 1],ylim,'k--'); xlim('tight');
    xticks([-180 0 180]);    title('Average Call Rate');    %yticks([]);
    nexttile;
    plotWinterval_AF_v0(rsz_bin_ctrs,avg_rsz_wbt_phase,avg_rsz_wbt_phase-sem_rsz_wbt_phase,avg_rsz_wbt_phase+sem_rsz_wbt_phase,'m');
    hold on;    plot(0*[1 1],ylim,'k--'); xlim('tight');
    xticks([-180 0 180]);    title('Wingbeat Phase');    %yticks([]);
    nexttile;
    plotWinterval_AF_v0(rsz_bin_ctrs,avg_rsz_wbt,avg_rsz_wbt-sem_rsz_wbt,avg_rsz_wbt+sem_rsz_wbt,'k');
    hold on;    plot(0*[1 1],ylim,'k--'); xlim('tight');
    xticks([-180 0 180]);    ylabel('Accelerometer');    %yticks([]);
    nexttile;
    plotWinterval_AF_v0(rsz_bin_ctrs,avg_rsz_wbt,avg_rsz_wbt-sem_rsz_wbt,avg_rsz_wbt+sem_rsz_wbt,'k');
    hold on;    plot(0*[1 1],ylim,'k--'); xlim('tight');
    xticks([-180 0 180]);    ylabel('Accelerometer');    %yticks([]);
    nexttile;
    plotWinterval_AF_v0(rsz_bin_ctrs,avg_rsz_wbt,avg_rsz_wbt-sem_rsz_wbt,avg_rsz_wbt+sem_rsz_wbt,'k');
    hold on;    plot(0*[1 1],ylim,'k--'); xlim('tight');
    xticks([-180 0 180]);    ylabel('Accelerometer');    %yticks([]);
    nexttile;
    plotWinterval_AF_v0(rsz_bin_ctrs,avg_rsz_wbt,avg_rsz_wbt-sem_rsz_wbt,avg_rsz_wbt+sem_rsz_wbt,'k');
    hold on;    plot(0*[1 1],ylim,'k--'); xlim('tight');
    xticks([-180 0 180]);    ylabel('Accelerometer');    %yticks([]);
    nexttile;
    plotWinterval_AF_v0(rsz_bin_ctrs,avg_rsz_wbt,avg_rsz_wbt-sem_rsz_wbt,avg_rsz_wbt+sem_rsz_wbt,'k');
    hold on;    plot(0*[1 1],ylim,'k--'); xlim('tight');
    xticks([-180 0 180]);    ylabel('Accelerometer');    %yticks([]);
    sgtitle(['Average of ',num2str(N_sweeps),' cycles']);
    
end

%% THETA SWEEPS ANALAYSIS (TEMPLATE MATCHING)
for hide=1
    
    %=== Select subtable of good flights
    SWPs_sst = SWPs(SWPs.rmsDec_error<1.5,:);
    N_flights = size(SWPs_sst,1);
    smooth_f = [1 .3];                                                  % [space bin,time bin]
    
    %=== Create template for finding candidate theta sweeps via template matching
    template_samples = 12;
    sweep_binspan = 6;
    conv_threshold = 2;
    template = normalize(interp1(normpdf(-3:0.1:3),linspace(1,61,template_samples)),'range')*bin_size_1D*sweep_binspan;
    % plot(template)
    
    %=== Init vectors
    SWPs_posterior = [];  SWPs_wbtphase = [];  SWPs_amp = [];
    
    %=== Find theta sweeps by template matching
    warning('off');
    for zz=1:N_flights
        
        %=== Number of spatial bins
        N_spatial_bins = size(SWPs_sst.p_dec_flight{zz,1},1);
        spt_bin_ids = [1:N_spatial_bins]';
        
        %=== Shifted posterior and wingbeat phase
        sft_posterior = imgaussfilt(SWPs_sst.p_dec_shifted{zz,1},smooth_f);
        wbt_phase = SWPs_sst.wbt_phase{zz,1};
        %wbt_phase = wrapTo2Pi(SWPs_sst.wbt_phase{zz,1})-pi;
        
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
        
        %=== For every sweep, store the posterior probability and the wingbeat phase
        for ss=1:numel(candidate_swp)
            if candidate_swp(ss)-template_samples>0 && candidate_swp(ss)+template_samples< size(sft_posterior,2) && round(N_spatial_bins/2)-1*sweep_binspan>0 && round(N_spatial_bins/2)+3*sweep_binspan<N_spatial_bins
                tmp_pst = sft_posterior(round(N_spatial_bins/2)+[-1*sweep_binspan:2*sweep_binspan],candidate_swp(ss)+[-template_samples:template_samples]);
                tmp_pst(isnan(tmp_pst))=0;
                SWPs_posterior = cat(3,SWPs_posterior,tmp_pst);
                SWPs_wbtphase = [SWPs_wbtphase; wbt_phase(candidate_swp(ss))];
                SWPs_amp = [SWPs_amp; candidate_swp_amp(ss)];
            end
        end
        %SWPs_sst.pos_real{zz,1}/SWPs_sst.pos_real{zz,1}(end);
        
    end
    warning('on');
    
    %=== Calculate a score for each sweep
    N_cand = size(SWPs_posterior,3);
    SWPs_score = zeros(N_cand,1);
    for i=1:N_cand
        tmp_pst = squeeze(SWPs_posterior(:,:,i));
        tmp_pst(isnan( tmp_pst)) = 0;
        [max_prb,max_loc] = max(tmp_pst,[],1);
        SWPs_score(i) = max(abs(diff(max_loc)));
    end
    
    %     %=== Plot a few example sweeps
    %     figure('units','normalized','outerposition',[0 0 1 1]);
    %     tiledlayout(10,30,'TileSpacing','compact');
    %     for i=1:min(N_cand,300)
    %         tmp_pst = squeeze(SWPs_posterior(:,:,randi(N_cand)));
    %         nexttile;   imagesc(tmp_pst,prctile(tmp_pst, [1 99],'all')');
    %         colormap(flipud(gray)); axis off;  set(gca,'YDir','normal');
    %         title(num2str(SWPs_score(i),2));
    %     end
    
    %====================================================================================================================================================
    
    % === TRAIN THE NEURAL NETWORK HERE IF NEEDED === %
    % theta_sweep_manual_class_AF_v0
    
    %====================================================================================================================================================
    %=== Load the pretrained network if needed
    if ~exist('trainedNet.mat','var')
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
    
%     %=== Show some examples
%     rng(3)
%     figure('units','normalized','outerposition',[0.1 0.2 .4 0.5]);
%     tiledlayout(10,20,'TileSpacing','tight');
%     for i=1:min(N_good,200)
%         tmp_pst = squeeze(SWPs_posterior(:,:,good_ids(randi(N_good))));
%         nexttile;   imagesc(tmp_pst,prctile(tmp_pst, [1 99],'all')');
%         colormap(flipud(gray)); axis off;  set(gca,'YDir','normal');
%     end
%     
%     sgtitle('Good Sweeps');
%     figure('units','normalized','outerposition',[.5 0 .5 1]);
%     tiledlayout(5,10,'TileSpacing','compact');
%     for i=1:min(N_baad,50)
%         tmp_pst = squeeze(SWPs_posterior(:,:,baad_ids(randi(N_baad))));
%         nexttile;   imagesc(tmp_pst,prctile(tmp_pst, [1 99],'all')');
%         colormap(flipud(gray)); axis off;  set(gca,'YDir','normal');
%     end
%     sgtitle('Bad Sweeps');
    
    %====================================================================================================================================================
    %=== Calculate average sweep posterior and phase distribution
    cond = good_ids;
    n_bins_phase = 11;  %Default 11
    bin_smoothing =5;
    phase_bins_tmp = linspace(-pi,pi,n_bins_phase);
    phase_ctrs_tmp = phase_bins_tmp(1:end-1)+mean(diff(phase_bins_tmp))/2;
    %phase2wbt_tmp = smoothdata(histcounts(SWPs_wbtphase(cond),phase_bins_tmp,'Normalization','probability'),'gaussian',bin_smoothing);
    phase2wbt_tmp = smoothdata(histcounts(SWPs_wbtphase(cond),phase_bins_tmp,'Normalization','probability'),'movmean',bin_smoothing);

    
    %=== Fit values
    x = [phase_ctrs_tmp,phase_ctrs_tmp+2*pi];
    y = repmat(phase2wbt_tmp,1,2);
    x1 = x(~isnan(y));  y1 = y(~isnan(y));
    y1 = normalize(y1-mean(y1),'range',[-1 1]);
    [f1,gof1] = fit(x1',y1',cosEqn,'Start',[1 0],'Lower',[-inf 0]);
    fit_coeff = coeffvalues(f1);
    fit_coeff_int = confint(f1);
    
    %=== Plot data
    figure('units','normalized','outerposition',[.3 .3 .3 .3]);
    tiledlayout(1,3,'TileSpacing','compact');
    nexttile;   imagesc(mean(SWPs_posterior(:,:, cond),3,'omitnan'),prctile(mean(SWPs_posterior(:,:, cond),3,'omitnan'),[1 99],'all')');    %hold on; plot(SWP_sst.max_loc{i,1});
    colormap(flipud(gray));
    axis off;  set(gca,'YDir','normal');
    nexttile;   polarhistogram(SWPs_wbtphase( cond),linspace(-pi,pi,18),'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','k');
    phase2plot = repmat(phase2wbt_tmp,1,2); delta_h = (max(phase2plot)-mean(phase2plot));
    nexttile;   bar([phase_ctrs_tmp,phase_ctrs_tmp+2*pi],phase2plot,1,'k','edgecolor','none');  ylim(mean(phase2plot)+delta_h*[-2 2]);
    xticks(-pi:pi:3*pi);    xticklabels({'-180', '0', '180', '360', '540'});   xlabel('Wingbeat');
    
    figure;
    plot(f1);    hold on;  bar(x1,y1,1,'EdgeColor','none','FaceColor','k'); legend('off');
    xticks(-pi:pi:3*pi);    xticklabels({'-180', '0', '180', '360', '540'});   
    title(['Max at: ', num2str(rad2deg(fit_coeff(2)),3), '( ', num2str( rad2deg(fit_coeff_int(2,1)),3),'-',num2str( rad2deg(fit_coeff_int(2,2)),3),')'])
    numel(cond)
   
end

%% ECHOLOCATION ANALYSIS
for hide=1
    
    %=== Look at echolocation phase
    N_flights = size(SWPs,1);
    click_phase = [];
    for zz=1:N_flights
        cond = [SWPs.pos_real{zz,1}/SWPs.pos_real{zz,1}(end)>0.15 & SWPs.pos_real{zz,1}/SWPs.pos_real{zz,1}(end)<0.85 & SWPs.clk{zz,1}];
        click_phase  = [click_phase; SWPs.wbt_phase{zz,1}(cond)];
    end
    
    %=== Plot phase
    figure('units','normalized','outerposition',[.3 .3 .2 .25]);
    tiledlayout(1,2,'TileSpacing','compact');
    nexttile;   polarhistogram(click_phase,linspace(-pi,pi,90),'Normalization','probability','facealpha',.9,'edgecolor','none','FaceColor','k');
    nexttile;   histogram([click_phase;click_phase+2*pi],unique([linspace(-pi,pi,40),linspace(pi,3*pi,40)]),'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','k');
    yticks([]); xticks(-pi:pi:3*pi);    xticklabels({'-180', '0', '180', '360', '540'});   xlabel('Wingbeat');
    
    %=== Select subtable of good flights for theta sweeps + echolocation
    SWPs_sst = SWPs(SWPs.rmsDec_error<2 & cellfun(@any,SWPs.clk) & SWPs.prc_decoded>0.5,:);
    
    %=== Params and initialize relevant variables
    n_reshape = 50;                 % Default 30
    smooth_f = [1 .3];              % [space bin,time bin]
    N_flights = size(SWPs_sst,1);
    
    single_SWP = table();
    counter=1;
    warning('off');
    %=== Cut at wingbeat maxima
    for zz=1:N_flights
        
        zero_phs_idx = find(SWPs_sst.wbt_phase{zz,1}(1:end-1).* SWPs_sst.wbt_phase{zz,1}(2:end)<0 & diff(SWPs_sst.wbt_phase{zz,1})<0);  % Segment based on wingbeat
        %zero_phs_idx = find(SWPs_sst.tmr_phase{zz,1}(1:end-1).* SWPs_sst.tmr_phase{zz,1}(2:end)<0 & diff(SWPs_sst.tmr_phase{zz,1})>0); % Segment based on Tamir's phase
        %plot(SWPs_sst.wbt_phase{zz,1});   hold on; stem(zero_phs_idx, ones(size(zero_phs_idx )));  plot(SWPs_sst.wbt{zz,1});
        sweep_strt = zero_phs_idx(1:end-1); sweep_stop = zero_phs_idx(2:end);
        N_spatial_bins = size(SWPs_sst.p_dec_flight{zz,1},1);
        spt_bin_ids = [1:N_spatial_bins]';
        
        for ss=1:numel(sweep_strt)
            
            single_SWP.flight_ID(counter) = zz;
            single_SWP.rms_dec(counter) = SWPs_sst.rmsDec_error(zz);
            single_SWP.prc_dec(counter) = SWPs_sst.prc_decoded(zz);
            single_SWP.reg_posterior(counter) = {imgaussfilt(SWPs_sst.p_dec_flight{zz,1}(:,sweep_strt(ss):sweep_stop(ss)),smooth_f)};
            single_SWP.sft_posterior(counter) = {imgaussfilt(SWPs_sst.p_dec_shifted{zz,1}(:,sweep_strt(ss):sweep_stop(ss)),smooth_f)};
            single_SWP.rsp_posterior(counter) = {imresize(single_SWP.sft_posterior{counter,1},[size(single_SWP.sft_posterior{counter,1},1),n_reshape])};
            
            single_SWP.spk_dsty(counter) = {SWPs_sst.spk_dsty{zz,1}(sweep_strt(ss):sweep_stop(ss))};
            single_SWP.rsz_spk_dsty(counter) = {interp1(single_SWP.spk_dsty{counter,1},linspace(1,numel(single_SWP.spk_dsty{counter,1}),n_reshape)')};
            single_SWP.wbt_power(counter) = {SWPs_sst.wbt_power{zz,1}(sweep_strt(ss):sweep_stop(ss))};
            single_SWP.LFP(counter) = {SWPs_sst.LFP{zz,1}(sweep_strt(ss):sweep_stop(ss))};
            single_SWP.rsz_LFP(counter) = {interp1(single_SWP.LFP{counter,1},linspace(1,numel(single_SWP.LFP{counter,1}),n_reshape)')};
            single_SWP.wbt(counter) = {SWPs_sst.wbt{zz,1}(sweep_strt(ss):sweep_stop(ss))};
            single_SWP.rsz_wbt(counter) = {interp1(single_SWP.wbt{counter,1},linspace(1,numel(single_SWP.wbt{counter,1}),n_reshape)')};
            single_SWP.fract_pos(counter) = {SWPs_sst.pos_real{zz,1}(sweep_strt(ss):sweep_stop(ss))/SWPs_sst.pos_real{zz,1}(end)};
            
            single_SWP.clk(counter) = {SWPs_sst.clk{zz,1}(sweep_strt(ss):sweep_stop(ss))};
            single_SWP.rsz_clk(counter) = {interp1(single_SWP.clk{counter,1},linspace(1,numel(single_SWP.clk{counter,1}),n_reshape)')};
            
            single_SWP.wbt_phase(counter) = {SWPs_sst.wbt_phase{zz,1}(sweep_strt(ss):sweep_stop(ss))};
            single_SWP.rsz_wbt_phase(counter) = {interp1(single_SWP.wbt_phase{counter,1},linspace(1,numel(single_SWP.wbt_phase{counter,1}),n_reshape)')};
            
            single_SWP.mean_spk_dsty(counter) = mean(single_SWP.spk_dsty{counter,1});
            single_SWP.mean_wbt_power(counter) = mean(single_SWP.wbt_power{counter,1});
            single_SWP.mean_fract_pos(counter) = mean(single_SWP.fract_pos{counter,1});
            
            [max_p,max_loc] = max(single_SWP.sft_posterior{counter,1},[],1);
            cnt_mass = spt_bin_ids'*single_SWP.sft_posterior{counter,1};
            single_SWP.med_jmp_distance(counter) = median(abs(diff(cnt_mass)));
            pos_sprd = 0;
            for bb=1:numel(cnt_mass)
                pos_sprd = pos_sprd+sqrt((spt_bin_ids'-cnt_mass(bb)).^2*single_SWP.sft_posterior{counter,1}(:,bb));
            end
            single_SWP.avg_pos_spread(counter) = pos_sprd/numel(cnt_mass);
            single_SWP.med_max_post(counter) = mean(max_p);
            single_SWP.max_loc(counter) = {max_loc-N_spatial_bins/2};
            single_SWP.cnt_mas(counter) = {cnt_mass-N_spatial_bins/2};
            single_SWP.est_dist(counter) = {interp1((cnt_mass-N_spatial_bins/2)*bin_size_1D,linspace(1,numel(cnt_mass),n_reshape)')};
            
            counter=counter+1;
        end
        
    end
    warning('on');
    
    no_match = 0;   % If matching the number of cycles for w/wo echolocation
    
    %=== Extract subset using defined criteria
    SWP_sst = single_SWP(single_SWP.mean_fract_pos>0.15 & single_SWP.mean_fract_pos<0.85 & single_SWP.mean_spk_dsty>prctile(single_SWP.mean_spk_dsty,20) & single_SWP.med_jmp_distance>0.2 & single_SWP.prc_dec>0.7 & single_SWP.rms_dec<1.5,:);
    SWP_sst_yclk = SWP_sst( cellfun(@any,SWP_sst.clk),:);
    SWP_sst_nclk = SWP_sst(~cellfun(@any,SWP_sst.clk),:);
    
    %=== Calculate averages and STDs
    N_sweeps = size(SWP_sst_yclk,1);
    est_dist = zeros(size(SWP_sst_yclk.est_dist{1,1},1),size(SWP_sst_yclk,1));
    rsz_spk_dsty = zeros(size(SWP_sst_yclk.rsz_spk_dsty{1,1},1),size(SWP_sst_yclk,1));
    rsz_wbt = zeros(size(SWP_sst_yclk.rsz_wbt{1,1},1),size(SWP_sst_yclk,1));
    rsz_clk = zeros(size(SWP_sst_yclk.rsz_clk{1,1},1),size(SWP_sst_yclk,1));
    for i=1:N_sweeps
        est_dist(:,i) = SWP_sst_yclk.est_dist{i,1};
        rsz_spk_dsty(:,i) = SWP_sst_yclk.rsz_spk_dsty{i,1};
        rsz_wbt(:,i) = SWP_sst_yclk.rsz_wbt{i,1};
        rsz_clk(:,i) = SWP_sst_yclk.rsz_clk{i,1};
    end
    avg_est_dist = mean(est_dist,2);                  avg_rsz_spk_dsty = mean(rsz_spk_dsty,2);                  avg_rsz_wbt = mean(rsz_wbt,2);                  avg_rsz_clk = mean(rsz_clk,2);
    sem_est_dist = std(est_dist,[],2)/sqrt(N_sweeps); sem_rsz_spk_dsty = std(rsz_spk_dsty,[],2)/sqrt(N_sweeps); sem_rsz_wbt = std(rsz_wbt,[],2)/sqrt(N_sweeps); sem_rsz_clk = std(rsz_clk,[],2)/sqrt(N_sweeps);
    
    %=== Calculate averages and STDs
    N_sweeps_nclk = size(SWP_sst_nclk,1);
    est_dist_nclk = zeros(size(SWP_sst_nclk.est_dist{1,1},1),size(SWP_sst_nclk,1));
    rsz_spk_dsty_nclk = zeros(size(SWP_sst_nclk.rsz_spk_dsty{1,1},1),size(SWP_sst_nclk,1));
    rsz_wbt_nclk = zeros(size(SWP_sst_nclk.rsz_wbt{1,1},1),size(SWP_sst_nclk,1));
    rsz_clk_nclk = zeros(size(SWP_sst_nclk.rsz_clk{1,1},1),size(SWP_sst_nclk,1));
    for i=1:N_sweeps_nclk
        est_dist_nclk(:,i) = SWP_sst_nclk.est_dist{i,1};
        rsz_spk_dsty_nclk(:,i) = SWP_sst_nclk.rsz_spk_dsty{i,1};
        rsz_wbt_nclk(:,i) = SWP_sst_nclk.rsz_wbt{i,1};
        rsz_clk_nclk(:,i) = SWP_sst_nclk.rsz_clk{i,1};
    end
    rng(3); %3
    if no_match
        avg_est_dist_nclk = mean(est_dist_nclk,2);                  avg_rsz_spk_dsty_nclk = mean(rsz_spk_dsty_nclk,2);                  avg_rsz_wbt_nclk = mean(rsz_wbt_nclk,2);                  avg_rsz_clk_nclk = mean(rsz_clk_nclk,2);
        sem_est_dist_nclk = std(est_dist_nclk,[],2)/sqrt(N_sweeps_nclk); sem_rsz_spk_dsty_nclk = std(rsz_spk_dsty_nclk,[],2)/sqrt(N_sweeps_nclk); sem_rsz_wbt_nclk = std(rsz_wbt_nclk,[],2)/sqrt(N_sweeps_nclk); sem_rsz_clk_nclk = std(rsz_clk_nclk,[],2)/sqrt(N_sweeps_nclk);
    else
        subsample = datasample(1:N_sweeps_nclk,N_sweeps,'Replace',false);
        avg_est_dist_nclk = mean(est_dist_nclk(:,subsample),2);                  avg_rsz_spk_dsty_nclk = mean(rsz_spk_dsty_nclk(:,subsample),2);                  avg_rsz_wbt_nclk = mean(rsz_wbt_nclk(:,subsample),2);                  avg_rsz_clk_nclk = mean(rsz_clk_nclk(:,subsample),2);
        sem_est_dist_nclk = std(est_dist_nclk(:,subsample),[],2)/sqrt(N_sweeps); sem_rsz_spk_dsty_nclk = std(rsz_spk_dsty_nclk(:,subsample),[],2)/sqrt(N_sweeps); sem_rsz_wbt_nclk = std(rsz_wbt_nclk(:,subsample),[],2)/sqrt(N_sweeps); sem_rsz_clk_nclk = std(rsz_clk_nclk(:,subsample),[],2)/sqrt(N_sweeps);
    end
    
    %=== Define bins
    rsz_bin_ctrs = linspace(-180,180,n_reshape);

    %=== Show averages
    figure('units','normalized','outerposition',[.3 .3 .3 .6]);
    tiledlayout(2,3,'TileSpacing','compact');
    nexttile;   plotWinterval_AF_v0(rsz_bin_ctrs  ,avg_est_dist,avg_est_dist-sem_est_dist,avg_est_dist+sem_est_dist,'k');
    hold on;    plot(0*[1 1],ylim,'k--'); xlim('tight');    title('Average Decoding Error');    xticks([]);
    plotWinterval_AF_v0(rsz_bin_ctrs  ,avg_est_dist_nclk,avg_est_dist_nclk-sem_est_dist_nclk,avg_est_dist_nclk+sem_est_dist_nclk,'k');
    nexttile;   plotWinterval_AF_v0(rsz_bin_ctrs  ,avg_rsz_spk_dsty,avg_rsz_spk_dsty-sem_rsz_spk_dsty,avg_rsz_spk_dsty+sem_rsz_spk_dsty,'r');
    hold on;    plot(0*[1 1],ylim,'k--'); xlim('tight');    title('Average Spike Density');  xticks([]);
    plotWinterval_AF_v0(rsz_bin_ctrs  ,avg_rsz_spk_dsty_nclk,avg_rsz_spk_dsty_nclk-sem_rsz_spk_dsty_nclk,avg_rsz_spk_dsty_nclk+sem_rsz_spk_dsty_nclk,'r');
    nexttile;   plotWinterval_AF_v0(rsz_bin_ctrs  ,avg_rsz_clk,avg_rsz_clk-sem_rsz_clk,avg_rsz_clk+sem_rsz_clk,'b');
    hold on;    plot(0*[1 1],ylim,'k--'); xlim('tight');    title('Average Call Rate');      xticks([]);
    plotWinterval_AF_v0(rsz_bin_ctrs  ,avg_rsz_clk_nclk,avg_rsz_clk_nclk-sem_rsz_clk_nclk,avg_rsz_clk_nclk+sem_rsz_clk_nclk,'b');
    nexttile;   plotWinterval_AF_v0(rsz_bin_ctrs  ,avg_rsz_wbt,avg_rsz_wbt-sem_rsz_wbt,avg_rsz_wbt+sem_rsz_wbt,'k');
    hold on;    plot(0*[1 1],ylim,'k--'); xlim('tight');
    xticks([-180 0 180]);    ylabel('Accelerometer');
    plotWinterval_AF_v0(rsz_bin_ctrs  ,avg_rsz_wbt_nclk,avg_rsz_wbt_nclk-sem_rsz_wbt_nclk,avg_rsz_wbt_nclk+sem_rsz_wbt_nclk,'k');
    nexttile;   plotWinterval_AF_v0(rsz_bin_ctrs  ,avg_rsz_wbt,avg_rsz_wbt-sem_rsz_wbt,avg_rsz_wbt+sem_rsz_wbt,'k');
    hold on;    plot(0*[1 1],ylim,'k--'); xlim('tight');
    xticks([-180 0 180]);    ylabel('Accelerometer');
    plotWinterval_AF_v0(rsz_bin_ctrs  ,avg_rsz_wbt_nclk,avg_rsz_wbt_nclk-sem_rsz_wbt_nclk,avg_rsz_wbt_nclk+sem_rsz_wbt_nclk,'k');
    nexttile;   plotWinterval_AF_v0(rsz_bin_ctrs  ,avg_rsz_wbt,avg_rsz_wbt-sem_rsz_wbt,avg_rsz_wbt+sem_rsz_wbt,'k');
    hold on;    plot(0*[1 1],ylim,'k--'); xlim('tight');
    xticks([-180 0 180]);    ylabel('Accelerometer');
    plotWinterval_AF_v0(rsz_bin_ctrs  ,avg_rsz_wbt_nclk,avg_rsz_wbt_nclk-sem_rsz_wbt_nclk,avg_rsz_wbt_nclk+sem_rsz_wbt_nclk,'k');
    sgtitle(['Average of ',num2str(N_sweeps),' cycles']);
    
end