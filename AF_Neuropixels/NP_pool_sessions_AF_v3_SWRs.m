%% POOL DATA RELATED TO SWRs
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
    Folder = cd;    FileList = dir(fullfile(Folder, '**', 'Analyzed_NPs_*'));   MS = [];
    for nc = 1:length(FileList)
        cd(FileList(nc).folder);
        load(FileList(nc).name);
        if ~isempty(MSC) && any(cellfun(@(row) isequal(row, MSC.unique_ID(1,:)), num2cell(sessions2include, 2)))
            MS = [MS; MSC];
        end
        disp([num2str(length(FileList)-nc),' remaining sessions to load...']);
    end
    cd(Folder); clearvars -except MS
    
    %=== Adjust the orientation of the Sharp-wave
    MS(2).ave_SWR_trace = -MS(2).ave_SWR_trace;
    MS(3).ave_SWR_trace = -MS(3).ave_SWR_trace;
    
    %=== If normalizing population firing by number of neurons
    normalize_by_neurons = 1;
    if normalize_by_neurons
        for i=1:size(MS,1)
            MS(i).PSTH_Spk =  MS(i).PSTH_Spk/MS(i).n_cells;
            MS(i).PSTH_Spk_CTRl =  MS(i).PSTH_Spk_CTRl/MS(i).n_cells;
            MS(i).PSTH_Spk_finer =  MS(i).PSTH_Spk_finer/MS(i).n_cells;
            MS(i).PSTH_Spk_CTRl_finer =  MS(i).PSTH_Spk_CTRl_finer/MS(i).n_cells;
        end
    end
    
end

%% LOOK AT SPIKING AND THETA BOUTS AROUND SWRs
for hide=1
    
%     bar(MS(1).time_bins(1:end-1)+diff(MS(1).time_bins)/2  ,mean(vertcat(MS.N_SWR_per_time)./sum(vertcat(MS.N_SWR_per_time),2)))
    
    %=== Plot all the average SWRs
    figure('units','normalized','outerposition',[.3 .3 .4 .3]);
    tiledlayout(1,3,'TileSpacing','tight');
    nexttile;
    plot(MS(1).ave_SWR_time,horzcat(MS.ave_SWR_trace),'Color',[.7 .7 .7]);  hold on;
    plot(MS(1).ave_SWR_time,mean(horzcat(MS.ave_SWR_trace),2),'Color',[0 0 0],'LineWidth',2);
    xlabel('Time (s)'); ylabel('LFP (uV)');
    nexttile;
    plot(MS(1).ave_SWR_time,normalize(horzcat(MS.ave_SWR_trace),1,'scale'),'Color',[.7 .7 .7]);  hold on;
    plot(MS(1).ave_SWR_time,mean(normalize(horzcat(MS.ave_SWR_trace),1,'scale'),2),'Color',[0 0 0],'LineWidth',2);
    xlabel('Time (s)'); ylabel('LFP (Norm.)');
    nexttile;
    data_tmp = normalize(horzcat(MS.ave_SWR_trace),1,'scale')';
    plotWinterval_AF_v0(MS(1).ave_SWR_time,mean(data_tmp,'omitnan'),mean(data_tmp,'omitnan')-std(data_tmp,[],'omitnan')./sqrt(size(data_tmp,1)),mean(data_tmp,'omitnan')+std(data_tmp,[],'omitnan')./sqrt(size(data_tmp,1)),'r');
    xlabel('Time (s)'); ylabel('LFP (Norm.)');
    sgtitle('Sharp Wave Ripples');
    
    %=== Plot all the average Firing Rates
    figure('units','normalized','outerposition',[.3 .3 .4 .3]);
    tiledlayout(1,2,'TileSpacing','tight');
    nexttile;
    data_tmp = vertcat(MS.PSTH_Spk_finer);
    plotWinterval_AF_v0(MS(1).PSTH_Spk_ctrs_finer,mean(data_tmp,'omitnan'),mean(data_tmp,'omitnan')-std(data_tmp,[],'omitnan')./sqrt(size(data_tmp,1)),mean(data_tmp,'omitnan')+std(data_tmp,[],'omitnan')./sqrt(size(data_tmp,1)),'r');    hold on;
    data_tmp = vertcat(MS.PSTH_Spk_CTRl_finer);
    plotWinterval_AF_v0(MS(1).PSTH_Spk_ctrs_finer,mean(data_tmp,'omitnan'),mean(data_tmp,'omitnan')-std(data_tmp,[],'omitnan')./sqrt(size(data_tmp,1)),mean(data_tmp,'omitnan')+std(data_tmp,[],'omitnan')./sqrt(size(data_tmp,1)),'k');
    xlabel('Time(s)');  ylabel('Population Firing Rate (Hz)');
    nexttile;
    data_tmp = normalize(vertcat(MS.PSTH_Spk_finer),2,'scale');
    plotWinterval_AF_v0(MS(1).PSTH_Spk_ctrs_finer,mean(data_tmp,'omitnan'),mean(data_tmp,'omitnan')-std(data_tmp,[],'omitnan')./sqrt(size(data_tmp,1)),mean(data_tmp,'omitnan')+std(data_tmp,[],'omitnan')./sqrt(size(data_tmp,1)),'r');    hold on;
    data_tmp = normalize(vertcat(MS.PSTH_Spk_CTRl_finer),2,'scale');
    plotWinterval_AF_v0(MS(1).PSTH_Spk_ctrs_finer,mean(data_tmp,'omitnan'),mean(data_tmp,'omitnan')-std(data_tmp,[],'omitnan')./sqrt(size(data_tmp,1)),mean(data_tmp,'omitnan')+std(data_tmp,[],'omitnan')./sqrt(size(data_tmp,1)),'k');
    xlabel('Time(s)');  ylabel('Population Firing Rate (Norm.)');
    
    %=== Look at theta bouts around SWRs
    figure('units','normalized','outerposition',[.3 .3 .3 .3]);
    tiledlayout(1,2,'TileSpacing','tight');
    nexttile;
    data_tmp = normalize(vertcat(MS.PSTH_ThtBts),2,'scale');
    plotWinterval_AF_v0(MS(1).PSTH_ThtBts_ctrs,mean(data_tmp,'omitnan'),mean(data_tmp,'omitnan')-std(data_tmp,[],'omitnan')./sqrt(size(data_tmp,1)),mean(data_tmp,'omitnan')+std(data_tmp,[],'omitnan')./sqrt(size(data_tmp,1)),'r');
    hold on;
    data_tmp = normalize(vertcat(MS.PSTH_ThtBts_CTRl),2,'scale');
    plotWinterval_AF_v0(MS(1).PSTH_ThtBts_ctrs,mean(data_tmp,'omitnan'),mean(data_tmp,'omitnan')-std(data_tmp,[],'omitnan')./sqrt(size(data_tmp,1)),mean(data_tmp,'omitnan')+std(data_tmp,[],'omitnan')./sqrt(size(data_tmp,1)),'k');
    nexttile;
    data_tmp = vertcat(MS.PSTH_ThtBts);
    plotWinterval_AF_v0(MS(1).PSTH_ThtBts_ctrs,mean(data_tmp,'omitnan'),mean(data_tmp,'omitnan')-std(data_tmp,[],'omitnan')./sqrt(size(data_tmp,1)),mean(data_tmp,'omitnan')+std(data_tmp,[],'omitnan')./sqrt(size(data_tmp,1)),'r');
    hold on;
    data_tmp = vertcat(MS.PSTH_ThtBts_CTRl);
    plotWinterval_AF_v0(MS(1).PSTH_ThtBts_ctrs,mean(data_tmp,'omitnan'),mean(data_tmp,'omitnan')-std(data_tmp,[],'omitnan')./sqrt(size(data_tmp,1)),mean(data_tmp,'omitnan')+std(data_tmp,[],'omitnan')./sqrt(size(data_tmp,1)),'k');
    
    
    %=== Quantify significance in time windows leading to the SWR
    bin_edges = linspace(MS(1).PSTH_Spk_ctrs(1),MS(1).PSTH_Spk_ctrs(end),10);
    bin_interval = mean(diff(bin_edges)/2);
    bin_centers = bin_edges(1:end-1)+bin_interval;
    time_fs = 1/mean(diff(MS(1).PSTH_Spk_ctrs));
    smoothed_trace = movmean(vertcat(MS.PSTH_Spk)',round(time_fs*bin_interval/1),1);
    PSTH_Spk_discr = interp1(MS(1).PSTH_Spk_ctrs,smoothed_trace,bin_centers);
    baseline = PSTH_Spk_discr(1,:);
    
    %=== Show firing rate, SWR and significance
    figure('units','normalized','outerposition',[.3 .3 .2 .5]);
    tiledlayout(2,1,'TileSpacing','tight');
    nexttile;
    data_tmp = vertcat(MS.PSTH_Spk);
    plotWinterval_AF_v0(MS(1).PSTH_Spk_ctrs,mean(data_tmp,'omitnan'),mean(data_tmp,'omitnan')-std(data_tmp,[],'omitnan')./sqrt(size(data_tmp,1)),mean(data_tmp,'omitnan')+std(data_tmp,[],'omitnan')./sqrt(size(data_tmp,1)),'k');  hold on;   
    ave_text_lvl = mean(ylim);  ylimmax = ylim;
    plot(repmat(bin_edges,2,1),repmat(ylim,numel(bin_edges),1)','k');   ylim([0 ylimmax(2)]);
    for i=1:size(PSTH_Spk_discr,1)
        p = signrank(PSTH_Spk_discr(1,:),PSTH_Spk_discr(i,:));
        if p > .05
            text(bin_centers(i),ave_text_lvl,'ns','HorizontalAlignment','center');
        elseif p<.05 && p>.01
            text(bin_centers(i),ave_text_lvl,'*','HorizontalAlignment','center');
        elseif p<.01 && p>.001
            text(bin_centers(i),ave_text_lvl,'**','HorizontalAlignment','center');
        elseif p<.001
            text(bin_centers(i),ave_text_lvl,'***','HorizontalAlignment','center');
        end
    end
    xticks([]);  ylabel('Population Firing Rate (Hz/neuron)');
    nexttile;
    data_tmp = normalize(horzcat(MS.ave_SWR_trace),1,'scale')';
    plotWinterval_AF_v0(MS(1).ave_SWR_time,mean(data_tmp,'omitnan'),mean(data_tmp,'omitnan')-std(data_tmp,[],'omitnan')./sqrt(size(data_tmp,1)),mean(data_tmp,'omitnan')+std(data_tmp,[],'omitnan')./sqrt(size(data_tmp,1)),'r');
    xlabel('Time (s)'); ylabel('LFP (Norm.)');
    sgtitle('Sharp Wave Ripples');
    
    %=== Similar quantification for theta bouts
    bin_edges = linspace(MS(1).PSTH_ThtBts_ctrs(1),MS(1).PSTH_ThtBts_ctrs(end),10);
    bin_interval = mean(diff(bin_edges)/2);
    bin_centers = bin_edges(1:end-1)+bin_interval;
    time_fs = 1/mean(diff(MS(1).PSTH_ThtBts_ctrs));
    smoothed_trace = movmean(vertcat(MS.PSTH_ThtBts)',round(time_fs*bin_interval/2),1);
    PSTH_ThtBts_discr = interp1(MS(1).PSTH_ThtBts_ctrs,smoothed_trace,bin_centers);
    smoothed_trace = movmean(vertcat(MS.PSTH_ThtBts_CTRl)',round(time_fs*bin_interval/2),1);
    PSTH_ThtBts_CTRl_discr = interp1(MS(1).PSTH_ThtBts_ctrs,smoothed_trace,bin_centers);
    
    figure('units','normalized','outerposition',[.3 .3 .17 .4]);
    nexttile;
    data_tmp = vertcat(MS.PSTH_ThtBts);
    plotWinterval_AF_v0(MS(1).PSTH_ThtBts_ctrs,mean(data_tmp,'omitnan'),mean(data_tmp,'omitnan')-std(data_tmp,[],'omitnan')./sqrt(size(data_tmp,1)),mean(data_tmp,'omitnan')+std(data_tmp,[],'omitnan')./sqrt(size(data_tmp,1)),'r');  hold on;
    data_tmp = vertcat(MS.PSTH_ThtBts_CTRl);
    plotWinterval_AF_v0(MS(1).PSTH_ThtBts_ctrs,mean(data_tmp,'omitnan'),mean(data_tmp,'omitnan')-std(data_tmp,[],'omitnan')./sqrt(size(data_tmp,1)),mean(data_tmp,'omitnan')+std(data_tmp,[],'omitnan')./sqrt(size(data_tmp,1)),'k');  hold on;
    ylim_fix = ylim;
    ave_text_lvl = mean(ylim);
    plot(repmat(bin_edges,2,1),repmat(ylim,numel(bin_edges),1)','k')
    for i=1:size(PSTH_ThtBts_CTRl_discr,1)
        p = signrank(PSTH_ThtBts_discr(i,:),PSTH_ThtBts_CTRl_discr(i,:));
        if p > .05
            text(bin_centers(i),ave_text_lvl,'ns','HorizontalAlignment','center');
        elseif p<.05 && p>.01
            text(bin_centers(i),ave_text_lvl,'*','HorizontalAlignment','center');
        elseif p<.01 && p>.001
            text(bin_centers(i),ave_text_lvl,'**','HorizontalAlignment','center');
        elseif p<.001
            text(bin_centers(i),ave_text_lvl,'***','HorizontalAlignment','center');
        end
    end
    xlabel('Time relative to SWR (s)');  ylabel('Theta Bouts Rate (Hz)');
    
    %=== Add bar plot for average theta bout rate at the moment 
    PSTH_ThtBts = vertcat(MS.PSTH_ThtBts);
    
    mean(PSTH_ThtBts,'all')
    std(PSTH_ThtBts,[],'all')/sqrt(size(PSTH_ThtBts,1))
    
    PSTH_ThtBts_CTRl = vertcat(MS.PSTH_ThtBts_CTRl);
    signrank(mean(PSTH_ThtBts(:,MS(1).PSTH_ThtBts_ctrs<0),2), mean(PSTH_ThtBts_CTRl(:,MS(1).PSTH_ThtBts_ctrs<0),2))
    data_tmp = [mean(PSTH_ThtBts(:,MS(1).PSTH_ThtBts_ctrs<0),2),mean(PSTH_ThtBts_CTRl(:,MS(1).PSTH_ThtBts_ctrs<0),2)];
    figure('units','normalized','outerposition',[.3 .3 .1 .4]);
    boxchart(data_tmp); ylim(ylim_fix);
    
    mean(mean(PSTH_ThtBts_CTRl(:,MS(1).PSTH_ThtBts_ctrs<0),2),'all')
    std(mean(PSTH_ThtBts_CTRl(:,MS(1).PSTH_ThtBts_ctrs<0),2),[],'all')/sqrt(size(PSTH_ThtBts,1))
    
end