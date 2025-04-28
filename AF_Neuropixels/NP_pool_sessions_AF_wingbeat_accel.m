%% POOL REPLAYS AND QUANTIFY INTERESTING FEATURES

%=== Load Data
for hide=1
    clr;
    
    %=== Default settings
    set(groot,'defaultAxesFontSize',12);
    set(groot,'defaultHistogramEdgeColor','none','defaultHistogramFaceAlpha',0.5);
    
    %=== Params for cross-correlation Replay SWR
    min_n = 5;                      % Minimum number of active cells
    min_f = 0.2;                    % Minimum fraction of active cells
    min_p = 0.05;                   % Minimum p value
    min_C = 0.2;                    % Minimum Correlation
    
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
        'Dataset_2','14611','231212';...
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
    Folder = cd;    FileList = dir(fullfile(Folder, '**', 'Analyzed_NPs_*'));  MS = [];
    for nc = 1:length(FileList)
        cd(FileList(nc).folder);
        load(FileList(nc).name);
        if ~isempty(MSC) && any(cellfun(@(row) isequal(row, MSC.unique_ID(1,:)), num2cell(sessions2include, 2)))
            
            %=== Accumulate the replay table and single units
            MS = [MS; MSC];
            
        end
        disp([num2str(length(FileList)-nc),' remaining sessions to load...']);
    end
    cd(Folder); clearvars -except MS
    
end

%% WINGBEAT ANALYSIS
for hide=1
    
    %=== Convert to table and add ID
    MS_table = struct2table(MS);
    n_sessions = size(MS_table,1);
    [unique_IDs,~,MS_table.bat_id] = unique(MS_table.unique_ID(:,2));
    Fs_Hi = 500;
    m_lc = 1.8; q_lc = -0.5;    % For classifying flights into straight or loops
    %bat_color = hsv(numel(unique_IDs));
    bat_color = lines(numel(unique_IDs));

%     %=== Scatter of wingbeat frequencies, color-coded by bat
%     figure('units','normalized','outerposition',[.3 .1 .1 .4]);
%     for i=1:max(MS_table.bat_id)
%         scatter(i*ones(sum(MS_table.bat_id==i),1),MS_table.wbf_all(MS_table.bat_id==i),'filled');    hold on;
%     end
%     xlim([0 7]);    ylim([0 10]);   xlabel('Bat ID');   ylabel('Average Wingbeat Frequency (Hz)');
    
    %=== Calculate average wingbeat frequency for each flight
    fg = [];    fl = [];    medc = [];    r2 = [];    maxc = [];    fc = [];    wbt_freq_rsp = [];  col = [];   batID = []; vel = [];   wbt_absv_rsp =[];
    ses_id = [];    cls_id = [];
    for i=1:n_sessions
        
        f_clus = MS_table.f_clus(i,1);
        col_tmp = bat_color(MS_table.bat_id(i),:);
        
        for jj = 1:f_clus.N
            
            if all(~isnan(f_clus.a_abs_NP{1,jj})) && f_clus.id(jj)>1    %=== Exclude unclustered
                
                %=== Accumulate values
                fg = [fg;   f_clus.wb_frequency_FFT(jj)];
                fl = [fl;   f_clus.wb_frequency_HLB(jj)];
                medc = [medc;   f_clus.med_curv2(jj)];
                r2 = [r2;   f_clus.Rsq1(jj)];
                col = [col; col_tmp];
                batID = [batID; MS_table.bat_id(i)];
                vel = [vel; median(f_clus.v_abs_Hi{1,jj})];
                ses_id = [ses_id;   str2double(MS_table.unique_ID(i,3))];
                cls_id = [cls_id;   f_clus.id(jj)];
                
                %=== Assign flight class
                fc = [fc; f_clus.med_curv2(jj)>q_lc+m_lc*f_clus.Rsq1(jj)];
                
                %=== Look at average wingbeat frequency profile
                int_rsp = 20:80;
                wbt_freq_rsp = [wbt_freq_rsp; interp1(linspace(1,100,numel(f_clus.wbt_freq{1,jj})),f_clus.wbt_freq{1,jj},int_rsp,'linear',9)];
                wbt_absv_rsp = [wbt_absv_rsp; interp1(linspace(1,100,numel(f_clus.v_abs_Hi{1,jj})),f_clus.v_abs_Hi{1,jj},int_rsp,'linear',9)];
                
            end
        end
        
    end
    
    %=== Show distributions
    figure('units','normalized','outerposition',[.3 .3 .3 .3]);
    tiledlayout(1,3,'TileSpacing','tight');
    nexttile;   scatter(fg,fl,10,'filled','MarkerFaceAlpha',0.3);    hold on; plot([6 10],[6 10],'k');  xlim([6 10]); axis('equal');  xlabel('WB Frequency FFT (Hz)');    ylabel('WB Frequency Hilbert (Hz)'); 
    [corr_val,p_val] = corr(fg,fl,'type','Pearson');  % Spearman correlation (several outliers)
    title([corr_val,p_val]);
    nexttile;   histogram(fg,[6:0.1:10],'Normalization','probability'); xlabel('Wingbeat Frequency FFT (Hz)');  ylabel('Fraction');  
    nexttile;   histogram(fl,[6:0.1:10],'Normalization','probability'); xlabel('Wingbeat Frequency HLB (Hz)');  ylabel('Fraction');  
    
    %=== Show a few features
    figure('units','normalized','outerposition',[0 .3 1 .3]);
    tiledlayout(1,8,'TileSpacing','tight');
    nexttile;   scatter(fg,fl,10,col,'filled','MarkerFaceAlpha',0.3);    hold on; plot([6 10],[6 10],'k');  xlim('tight');  axis('equal');  xlabel('WB Frequency FFT (Hz)');    ylabel('WB Frequency Hilbert (Hz)');
    nexttile;   gscatter(r2,medc,fc); hold on; plot([0 1],q_lc+m_lc*[0 1],'k--');   legend('off');  
    xlabel('R2');   ylabel('Median Curvature @middle flight (1/m)');    ylim([0 prctile(medc,99)]);  
    nexttile;   histogram(fg,[6:0.1:10],'Normalization','probability'); xlabel('Wingbeat Frequency (Hz)');  ylabel('Fraction'); title(std(fg)/mean(fg));   
    nexttile;   plot_distr_AF_v0(fg(find(fc)), fg((find(~fc))), {'Loops', 'Straight'}, 'SEM', 'Accelerometer'); %ylim_vals = ylim;   ylim([0 ylim_vals(2)]);
    nexttile;   data_tmp = wbt_freq_rsp(find(fc),:);
    plotWinterval_AF_v0(int_rsp,mean(data_tmp,'omitnan'),mean(data_tmp,'omitnan')-std(data_tmp,[],'omitnan')./sqrt(size(data_tmp,1)),mean(data_tmp,'omitnan')+std(data_tmp,[],'omitnan')./sqrt(size(data_tmp,1)),'b');    hold on;
    data_tmp = wbt_freq_rsp(find(~fc),:);
    plotWinterval_AF_v0(int_rsp,mean(data_tmp,'omitnan'),mean(data_tmp,'omitnan')-std(data_tmp,[],'omitnan')./sqrt(size(data_tmp,1)),mean(data_tmp,'omitnan')+std(data_tmp,[],'omitnan')./sqrt(size(data_tmp,1)),'r');    
    nexttile;   data_tmp = wbt_absv_rsp(find(fc),:);
    plotWinterval_AF_v0(int_rsp,mean(data_tmp,'omitnan'),mean(data_tmp,'omitnan')-std(data_tmp,[],'omitnan')./sqrt(size(data_tmp,1)),mean(data_tmp,'omitnan')+std(data_tmp,[],'omitnan')./sqrt(size(data_tmp,1)),'b');    hold on;
    data_tmp = wbt_absv_rsp(find(~fc),:);
    plotWinterval_AF_v0(int_rsp,mean(data_tmp,'omitnan'),mean(data_tmp,'omitnan')-std(data_tmp,[],'omitnan')./sqrt(size(data_tmp,1)),mean(data_tmp,'omitnan')+std(data_tmp,[],'omitnan')./sqrt(size(data_tmp,1)),'r');    
    nexttile;  scatter(fg,vel,15,col,'filled','MarkerFaceAlpha',0.3);   xlim([7 9]);   ylim(prctile(vel,[5 99]));   ylabel('Average Speed (ms-1)'); xlabel('Wingbeat Frequency (Hz)');  axis('square');    
    [corr_val,p_val] = corr(zscore(fg),zscore(vel),'type','Spearman');  % Spearman correlation (several outliers)
    title([corr_val,p_val]);
    
    %=== Show histograms for all bats
    for i=1:numel(unique_IDs),mean_tmp(i) = mean(fg(batID==i));end
    [~,tmp_idx] = sort(mean_tmp);
    figure('units','normalized','outerposition',[.1 .1 .2 .7]);
    tiledlayout(numel(unique_IDs),2,'TileSpacing','none');
    for i=1:numel(unique_IDs)
        nexttile; histogram(fg(batID==tmp_idx(i)),[6:0.1:10],'Normalization','probability','FaceColor',bat_color(tmp_idx(i),:),'EdgeColor','none'); xlim('tight'); 
        if i<numel(unique_IDs),xticks([]);
        else,xlabel('Wingbeat Frequency (Hz)');end
        nexttile; histogram(vel(batID==tmp_idx(i)),[0:0.1:5],'Normalization','probability','FaceColor',bat_color(tmp_idx(i),:),'EdgeColor','none'); xlim('tight'); 
        if i<numel(unique_IDs),xticks([]);
        else,xlabel('Velocity (ms-1)');end
    end
    
    %=== Compare GLMs with wingbeat eplained by bat and cluster
    GLM_table = table(fg,batID,ses_id,cls_id);
    [~,~,GLM_table.cls_abs_id] = unique(GLM_table(:,2:end));
    GLM_table.batID = categorical(GLM_table.batID);             % Ensure batID is categorical
    GLM_table.cls_abs_id = categorical(GLM_table.cls_abs_id);   % Ensure cls_abs_id is categorical
    GLM_table.fg_z = zscore(GLM_table.fg);                        % zscore Wingbeat frequency
    
    % Fit a stepwise GLM starting with an intercept-only model
    stepwise_model = stepwiseglm(GLM_table(:,[1,2,5]), 'constant', ...
        'Distribution', 'normal', ... % Assuming normal residuals
        'Upper', 'interactions'); % Allows only linear terms
    
    % Fit a stepwise GLM starting with a bat ID only model
    stepwise_model = stepwiseglm(GLM_table(:,[1,2,5]), 'fg ~ batID ', ...
        'Distribution', 'normal', ... % Assuming normal residuals
        'Upper', 'interactions'); % Allows only linear terms
    
    %=== Look at coefficient of variation
    FG_summ = groupsummary(GLM_table,'cls_abs_id',{'mean','std'},{'fg'});
    FG_summ.cv = FG_summ.std_fg./FG_summ.mean_fg;
    cv_all = std(fg)/mean(fg);
    figure('units','normalized','outerposition',[.3 .3 .15 .3]);
    histogram(FG_summ.cv,[0:0.005:0.1]); hold on; plot(cv_all*[1 1],ylim);    xlabel('Coefficient of variation');
    title(1-sum(FG_summ.cv<cv_all)/numel(FG_summ.cv));
    
    %=== Show some example accelerometer profiles across flights of the same cluster
    figure('units','normalized','outerposition',[0 0 1 1]);
    tiledlayout(8,10,'TileSpacing','tight');
    int_rsp_1 = 1:0.5:100;
    for i=1:n_sessions
         f_clus = MS_table.f_clus(i,1);
         for j= setdiff(unique(f_clus.id),1)
             acc_mtx = [];   
             clus_idx = find(f_clus.id==j);
             for jj=1:numel(clus_idx)
                a_trace = f_clus.a_flt_NP{1,clus_idx(jj)};
                acc_mtx = [acc_mtx; interp1(linspace(1,100,numel(a_trace)),a_trace,int_rsp_1,'linear')];
             end
             nexttile;  imagesc(acc_mtx);   xticks([]); yticks([numel(clus_idx)]);  colormap('inferno');
         end
    end
    
    % Show some example accelerometer profiles across flights of the same cluster (NEW)
    figure('units','normalized','outerposition',[0 0 1 1]);
    tiledlayout(8,10,'TileSpacing','loose');
    int_rsp_1 = 1:0.5:100;
    for i=1:n_sessions
         f_clus = MS_table.f_clus(i,1);
         for j= setdiff(unique(f_clus.id),1)
             a_traces = {};   
             clus_idx = find(f_clus.id==j);
             for jj=1:numel(clus_idx)
                a_traces = [a_traces; {f_clus.a_flt_NP{1,clus_idx(jj)}}];
             end
             acc_mtx = shift_wbt_AF_v0(a_traces,500,0.5);
             nexttile;  imagesc([0 size(acc_mtx,1)./500],[],acc_mtx');   
             xticks([0 size(acc_mtx,1)./500]);    yticks([numel(clus_idx)]);  colormap('inferno');  title(i);
         end
    end
end
