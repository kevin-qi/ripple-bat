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

%% LOOK AT DELTA
for hide=1
    
    %=== Quantify theta and delta power during flight and rest
    tht_power_flt = arrayfun(@(s) mean(s.avg_tht_power( s.bflying_Hi)), LFP);
    tht_power_rst = arrayfun(@(s) mean(s.avg_tht_power(~s.bflying_Hi)), LFP);
    dlt_power_flt = arrayfun(@(s) mean(s.avg_dlt_power( s.bflying_Hi)), LFP);
    dlt_power_rst = arrayfun(@(s) mean(s.avg_dlt_power(~s.bflying_Hi)), LFP);
    dlt_power_zscored = arrayfun(@(s) zscore(s.avg_dlt_power), LFP, 'UniformOutput', false);
    dlt_power_lwrSD = mean(cellfun(@(x) sum(x<-1)/numel(x), dlt_power_zscored));
    dlt_power_hgrSD = mean(cellfun(@(x) sum(x>+1)/numel(x), dlt_power_zscored));
    delta_ratio = tht_power_rst./tht_power_flt;
    
    mean_and_sem_AF_v0(arrayfun(@(s) mean(s.tht_bouts_dur), LFP));
    
    %=== Show theta and delta power during flight and rest
    figure('units','normalized','outerposition',[.3 .1 .2 .3]);
    tiledlayout(1,2,'TileSpacing','tight');
    nexttile;   plot_distr_AF_v0(tht_power_flt, tht_power_rst, {'Flight', 'Rest'}, 'SEM', 'Theta Power');
    title(['p = ',num2str(signrank(tht_power_flt,tht_power_rst),3)]);
    nexttile;   plot_distr_AF_v0(dlt_power_flt, dlt_power_rst, {'Flight', 'Rest'}, 'SEM', 'Delta Power');
    title(['p = ',num2str(signrank(dlt_power_flt,dlt_power_rst),3)]);
    
    %=== Show the ratio
    figure('units','normalized','outerposition',[.3 .1 .1 .4]);
    boxchart(delta_ratio); ylabel('Delta ratio (flight/rest)'); ylim([0 4]);
    title(['p = ',num2str(signrank(delta_ratio-1),3)]);
    mean_and_sem_AF_v0(delta_ratio);
    
    %=== Show delta power around theta bouts vs control
    dlt_around_bout = arrayfun(@(s) mean(s.dlt_around_bout,2)', LFP, 'UniformOutput', false);   dlt_around_bout = vertcat(dlt_around_bout{:});
    dlt_around_ctrl = arrayfun(@(s) mean(s.dlt_around_ctrl,2)', LFP, 'UniformOutput', false);   dlt_around_ctrl = vertcat(dlt_around_ctrl{:});
    tht_around_bout = arrayfun(@(s) mean(s.tht_around_bout,2)', LFP, 'UniformOutput', false);   tht_around_bout = vertcat(tht_around_bout{:});
    tht_around_ctrl = arrayfun(@(s) mean(s.tht_around_ctrl,2)', LFP, 'UniformOutput', false);   tht_around_ctrl = vertcat(tht_around_ctrl{:});
    tme_around_bout = LFP(1).tme_around_bout;
    
    %=== Normalization
    dlt_norm_f = mean(dlt_around_ctrl,2);
    tht_norm_f = mean(tht_around_ctrl,2);
    dlt_around_bout = dlt_around_bout./dlt_norm_f;
    dlt_around_ctrl = dlt_around_ctrl./dlt_norm_f;
    tht_around_bout = tht_around_bout./tht_norm_f;
    tht_around_ctrl = tht_around_ctrl./tht_norm_f;
    
%     %=== Normalization
%     dlt_around_bout = normalize(dlt_around_bout,2,'medianiqr');
%     dlt_around_ctrl = normalize(dlt_around_ctrl,2,'medianiqr');
%     tht_around_bout = normalize(tht_around_bout,2,'medianiqr');
%     tht_around_ctrl = normalize(tht_around_ctrl,2,'medianiqr');
    
    %=== Show delta and theta power around theta bouts vs control
    figure('units','normalized','outerposition',[.3 .1 .35 .3]);
    tiledlayout(1,3,'TileSpacing','tight');
    nexttile;
    %data_tmp = normalize(dlt_around_bout,2,'range');
    data_tmp = dlt_around_bout;
    plotWinterval_AF_v0(tme_around_bout,mean(data_tmp,'omitnan'),mean(data_tmp,'omitnan')-std(data_tmp,[],'omitnan')./sqrt(size(data_tmp,1)),mean(data_tmp,'omitnan')+std(data_tmp,[],'omitnan')./sqrt(size(data_tmp,1)),'b');  hold on;
    %data_tmp = normalize(dlt_around_ctrl,2,'range');
    data_tmp = dlt_around_ctrl;
    plotWinterval_AF_v0(tme_around_bout,mean(data_tmp,'omitnan'),mean(data_tmp,'omitnan')-std(data_tmp,[],'omitnan')./sqrt(size(data_tmp,1)),mean(data_tmp,'omitnan')+std(data_tmp,[],'omitnan')./sqrt(size(data_tmp,1)),'k'); 
    plot([0 0], ylim, 'k--'); xlabel('Time from Theta Bout Onset (s)'); ylabel('Delta Power');
    nexttile;
    data_tmp = tht_around_bout;
    plotWinterval_AF_v0(tme_around_bout,mean(data_tmp,'omitnan'),mean(data_tmp,'omitnan')-std(data_tmp,[],'omitnan')./sqrt(size(data_tmp,1)),mean(data_tmp,'omitnan')+std(data_tmp,[],'omitnan')./sqrt(size(data_tmp,1)),'r');  hold on;
    %data_tmp = normalize(dlt_around_ctrl,2,'range');
    data_tmp = tht_around_ctrl;
    plotWinterval_AF_v0(tme_around_bout,mean(data_tmp,'omitnan'),mean(data_tmp,'omitnan')-std(data_tmp,[],'omitnan')./sqrt(size(data_tmp,1)),mean(data_tmp,'omitnan')+std(data_tmp,[],'omitnan')./sqrt(size(data_tmp,1)),'k'); 
    plot([0 0], ylim, 'k--'); xlabel('Time from Theta Bout Onset (s)'); ylabel('Delta Power');
    nexttile;
    data_tmp = tht_around_bout;
    plotWinterval_AF_v0(tme_around_bout,mean(data_tmp,'omitnan'),mean(data_tmp,'omitnan')-std(data_tmp,[],'omitnan')./sqrt(size(data_tmp,1)),mean(data_tmp,'omitnan')+std(data_tmp,[],'omitnan')./sqrt(size(data_tmp,1)),'r');  hold on;
    %data_tmp = normalize(dlt_around_ctrl,2,'range');
    data_tmp = dlt_around_bout;
    plotWinterval_AF_v0(tme_around_bout,mean(data_tmp,'omitnan'),mean(data_tmp,'omitnan')-std(data_tmp,[],'omitnan')./sqrt(size(data_tmp,1)),mean(data_tmp,'omitnan')+std(data_tmp,[],'omitnan')./sqrt(size(data_tmp,1)),'b'); 
    plot([0 0], ylim, 'k--'); xlabel('Time from Theta Bout Onset (s)'); ylabel('Delta Power');
    
end
