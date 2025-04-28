%% POOL DATA RELATED TO LFP AND THETA SWEEPS
%=== Load data
for hide=1
    clr;
    
    %=== Default settings
    set(groot,'defaultAxesFontSize',12);
    set(groot,'defaultHistogramEdgeColor','none','defaultHistogramFaceAlpha',0.5);
    
    %=== Sessions to load (comment sessions you don't want to include)
    sessions2include = {...
        'Dataset_4','14521','250319';...
        'Dataset_4','14521','250320';...
        'Dataset_4','14521','250321';...
        'Dataset_4','14521','250322';...
        'Dataset_4','14521','250323';...
        'Dataset_4','14640','250319';...
        %'Dataset_4','14640','250320';...
        'Dataset_4','14640','250321';...
        'Dataset_4','14640','250322';...
        'Dataset_4','14640','250323';...
        'Dataset_4','14640','250324';...
        
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
        };

    
    %=== Load data and aggregate them in the Multi_Day structure
    Folder = cd;    FileList = dir(fullfile(Folder, '**', 'Analyzed_NPs_*'));   LFP = [];   NP_units = [];  SWPs = []; ids = [];
    for nc = 1:length(FileList)
        cd(FileList(nc).folder);
        load(FileList(nc).name);
        
        if ~isfield(NP_unit,'fclus')
            [NP_unit.f_clus] = deal(NaN);
        end
        
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
    
    %=== Plot the PSDs (single session averages)
    figure('units','normalized','outerposition',[0.1 0.3 0.3 0.3]);
    tiledlayout(1,3,'TileSpacing','tight');
    nexttile;   plot(freq,smoothdata(PSD_f',1,'gaussian',1/freq_smp),'LineWidth',2); set(gca, 'YScale', 'log');   xlim([1 20]);   xlabel('Frequency (Hz)');   ylabel('PSD (Norm.)');    title('LFP (flight)');
    nexttile;   plot(freq,smoothdata(PSD_r',1,'gaussian',1/freq_smp),'LineWidth',2); set(gca, 'YScale', 'log');   xlim([1 20]);   xlabel('Frequency (Hz)');   ylabel('PSD (Norm.)');    title('LFP (rest)');
    nexttile;   plot(freq,smoothdata(PSD_a',1,'gaussian',1/freq_smp),'LineWidth',2); set(gca, 'YScale', 'log');   xlim([1 20]);   xlabel('Frequency (Hz)');   ylabel('PSD (Norm.)');    title('Accelerometer');
    
    figure('units','normalized','outerposition',[0 0 .5 1]);
    tiledlayout(5,6,'TileSpacing','tight');
    for i=1:size(PSD_f',2)
    
        nexttile;   hold on; 
        plot(freq,smoothdata(PSD_r(i,:),'gaussian',100),'k','LineWidth',2); 
        plot(freq,smoothdata(PSD_f(i,:),'gaussian',100),'r','LineWidth',2);
        %plot(freq,smoothdata(PSD_b(i,:),'gaussian',200),'b','LineWidth',2);
        %plot(freq,smoothdata(PSD_pst(i,:),'gaussian',200),'g','LineWidth',2); 
        %plot(freq,smoothdata(PSD_a(i,:),'gaussian',50),'r','LineWidth',3);
        
        xlim([1 20]);   %ylim([1e-5 1e-3]);  
        xlabel('Frequency (Hz)');   ylabel('PSD (Norm.)');    title(LFP(i).unique_ID  );
        set(gca, 'YScale', 'log');   
         
    end
   
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
    signrank(pwr_b, pwr_f)
    disp(['Fraction bouts during rest: ', num2str(1-mean([LFP.fract_bouts_flight]),3)]);
    std(1-[LFP.fract_bouts_flight])./sqrt(numel([LFP.fract_bouts_flight]))
    
    figure
    histogram(vertcat(LFP_table.tht_bouts_dur{:}),'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','k');
    
    %=== Look at Tamir's power
    figure('units','normalized','outerposition',[0.1 0.2 0.1 0.4]);
    plot_distr_AF_v0(LFP_table.tmr_power(:,1), LFP_table.tmr_power(:,2), {'Flight', 'Rest'}, 'SEM', 'Median Non-Oscillatory Power');  hold on; 
    title(['p = ',num2str(signrank(LFP_table.tmr_power(:,1), LFP_table.tmr_power(:,2)),3)]);
    
    %=== Look at spike-tgd LFP
    figure('units','normalized','outerposition',[0.3 0.2 0.15 0.3]);
    hold on;
    data = vertcat(LFP_table.spkTgdLfP_all{:});
    plotWinterval_AF_v0(LFP_table.spkTgdLfP_time{1,1},mean(data),mean(data)-std(data)./sqrt(size(data,1)),mean(data)+std(data)./sqrt(size(data,1)),'k');
    data = vertcat(LFP_table.spkTgdLfP_flight{:});
    plotWinterval_AF_v0(LFP_table.spkTgdLfP_time{1,1},mean(data),mean(data)-std(data)./sqrt(size(data,1)),mean(data)+std(data)./sqrt(size(data,1)),'r');
    xlim('tight');  xlabel('Time from spike (s)');  ylabel('LFP (uV)'); title('Spike Triggered LFP (All sessions)');    legend('All spikes','','Flight');
    
    %=== Look at LFP spectrogram around takeoff (EXAMPLES)
    t_around_flight = LFP(1).interval{1, 1};
    LFP_around_flight_all = [];
    WBT_around_flight_all = [];
    for i=1:size(LFP,1)
        LFP_around_flight_all = [LFP_around_flight_all; LFP(i).LFPATflight{:}];
        WBT_around_flight_all = [WBT_around_flight_all; LFP(i).accATflight{:}];
    end
    f_num_Hi = size(LFP_around_flight_all,1); 
   
    figure('units','normalized','outerposition',[0 0 1 1]);
    tiledlayout(7,10,'TileSpacing','compact');
    rng(3)
    for i=1:min(f_num_Hi,35)
        rand_idx = randi(size(LFP_around_flight_all,1));
        nexttile;   plot(t_around_flight,normalize(LFP_around_flight_all(rand_idx,:)),'r');
        hold on;    plot(t_around_flight,normalize(WBT_around_flight_all(rand_idx,:))+5,'k');
        plot([0 0],ylim,'k--'); yticks([]);
        if i==min(f_num_Hi,70),xlabel('Time (s)');else,xticks([]); end
        if i==1,ylabel('Flight aligned LFP');end
        nexttile;
        [PS_LFP,freq_SG] = cwt(LFP_around_flight_all(rand_idx,:),LFP(1).Fs_Hi);
        imagesc([t_around_flight(1) t_around_flight(end)],[freq_SG(1),freq_SG(end)],imgaussfilt(abs(PS_LFP),[.5 100])); shading interp;  colormap(hot); 
        hold on;    plot(xlim,[8 8],'w--'); plot([0 0],[0.1 20],'w--'); 
        set(gca, 'YScale', 'log','YDir', 'normal','TickLength',[0 0]);   ylim([0.5 20]);
        yticks([1 5 10 20]); ylabel('Freq (Hz)');
    end
    %colorbar
    
    CWTL_around_flight_all = [];
    CWTW_around_flight_all = [];
    for i=1:f_num_Hi
        [PS_LFP,freq_SG] = cwt(LFP_around_flight_all(i,:),LFP(1).Fs_Hi);
        [PS_WBT,freq_SG] = cwt(WBT_around_flight_all(i,:),LFP(1).Fs_Hi);
        CWTL_around_flight_all = cat(3,CWTL_around_flight_all,abs(PS_LFP));
        CWTW_around_flight_all = cat(3,CWTW_around_flight_all,abs(PS_WBT));
    end
    
    figure('units','normalized','outerposition',[.2 .3 .15 .5]);
    tiledlayout(2,1,'TileSpacing','compact');
    nexttile;
    imagesc([t_around_flight(1) t_around_flight(end)],[freq_SG(1),freq_SG(end)],imgaussfilt(squeeze(mean(CWTL_around_flight_all,3)),[.1 0.1])); shading interp;  colormap(redblue); 
    hold on;    plot(xlim,[8 8],'w--');     plot([0 0],[0.1 20],'w--');
    set(gca, 'YScale', 'log','YDir', 'normal','TickLength',[0 0]);   ylim([1 20]);
    yticks([1 5 10 20]); ylabel('Freq (Hz)');   xticks([]); title('LFP');
    nexttile;
    imagesc([t_around_flight(1) t_around_flight(end)],[freq_SG(1),freq_SG(end)],imgaussfilt(squeeze(mean(CWTW_around_flight_all,3)),[.1 0.1])); shading interp;  colormap(redblue); 
    hold on;    plot(xlim,[8 8],'w--');     plot([0 0],[0.1 20],'w--');
    set(gca, 'YScale', 'log','YDir', 'normal','TickLength',[0 0]);   ylim([1 20]);
    yticks([1 5 10 20]); ylabel('Freq (Hz)');   xlabel('Time from takeoff (s)');    title('Wingbeat');
    
end
