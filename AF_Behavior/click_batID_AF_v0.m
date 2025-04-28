function click_batID_AF_v0(folder_name)

%% LOAD DATA
for hide=1
    %=== TT data and related
    TTlFile = dir('TTL_*');                             load([TTlFile.folder '/' TTlFile.name]);        % Synchronization TTLs
    C3DFile = dir('*Bat_cluster.mat');                  load([C3DFile.folder '/' C3DFile.name]);        % C3d File
    load('imp_bat.mat');                                                                                % Tag of the implanted bat
    
    %=== Behavioral data and related
    BHV1_file = dir(fullfile(fileparts(TTlFile.folder),'Ext_Beh*','Extracted_Behavior*'));                              % Preprocessed Behavioral Data
    load([BHV1_file.folder,'\',BHV1_file.name],'a','a_abs','angle','bat_clr','bat_nms','bat_pair_nms','bat_pairs',...
        'batdate','bflying','Fs','f_num','Group_name','n_tags','r','r_lim','t','T','v','v_abs','v_th','wBeats');        % Processed Behavioral Data
    BHV2_file = dir(fullfile(fileparts(TTlFile.folder),'Ext_Beh*','Analysis*/Analyzed_Behavior*'));
    load([BHV2_file.folder,'\',BHV2_file.name]);
    
    %=== Echolocation Clicks
    if ~isempty(dir(fullfile(fileparts(TTlFile.folder),'Detected_clicks','Detected_Clicks_*')))
        Clck_file = dir(fullfile(fileparts(TTlFile.folder),'Detected_clicks','Detected_Clicks_*'));                     % File with times of the detected echolocation clicks
        load([Clck_file.folder,'\',Clck_file.name]);
    else
        Detected_Clicks =[];
        Detected_Clicks.times = [];
    end
end

%% PARAMETERS AND DEFINITIONS
for hide=1
    
    %=== Parameters and options
    options.savefigures = 1;                                                                                        % Save Figures
    fig_count = 1;                                                                                                  % Id of the first figure
   
    %=== Create analysis folder for storing the results
    if ~exist('folder_name')
        analysis_directory=fullfile(pwd,['CK&BHv_Analysis_',datestr(now, 'yymmdd_HHMM')]);
        if ~exist(analysis_directory,'dir');mkdir(analysis_directory);end
    else
        current_directory = pwd;
        analysis_directory = [replace(current_directory,'Ephys_Restrictive_Raw',folder_name),'\CK&BHv_Analysis_',datestr(now, 'yymmdd_HHMM')];
        if ~exist(analysis_directory,'dir');mkdir(analysis_directory);end
    end
end

%% LOOK AT ECHOLOCATION CLICKS (ALL FLIGHTS)
for hide=1
    
    %=== Extract all flights with 1 bat max flying
    cond = [FLIGHTS.maxN ==1];
    F_temp = FLIGHTS(all(cond,2),:);
    
    figure('units','normalized','outerposition',[0.2 0.3 0.2 0.6]);
    tiledlayout(4,3,'TileSpacing','tight');
    nexttile([3,1]);   [PSTH_mn,PSTH_cs] = Raster_AF_v3(Detected_Clicks.times,F_temp.t1,F_temp.id,bat_clr(F_temp.id,:),[-1 1],'A-Tko',1,1);
    nexttile([3,1]);   Raster_AF_v3(Detected_Clicks.times,F_temp.t2,F_temp.id,bat_clr(F_temp.id,:),[-1 1],'A-Lnd',1,1);
    nexttile([3,1]);   Raster_TimeWrap_AF_v1(Detected_Clicks.times,F_temp.t1,F_temp.t2,F_temp.id,bat_clr(F_temp.id,:),1,'T-Wrp',1);
    Click_rise = rescale(mean(PSTH_mn,1),0,1)';
    nexttile;          plot(PSTH_cs,Click_rise,'k','LineWidth',2);   hold on;    plot([0 0],ylim,'k--'); hold off;
    nexttile;          plot(PSTH_cs,Click_rise,'k','LineWidth',2);   hold on;    plot([0 0],ylim,'k--');
    plot(xlim,[.1 .1],'r--');   plot(PSTH_cs(knnsearch(Click_rise,0.1))*[1 1],ylim,'r--');  hold off;
    title(['T_r = ', num2str(PSTH_cs(knnsearch(Click_rise,0.1)),2), ' s']);
    nexttile;   plot(PSTH_cs,normalize(PSTH_mn,2,'range',[0 1]),'LineWidth',2); hold on;    plot([0 0],ylim,'k--'); hold off;
    sgtitle([batdate,' ',Group_name]);
    fig_count = saveFig(analysis_directory,batdate,fig_count,options.savefigures);
    
    %=== Quantify bat-specific aspects of echolocation
    bat_click_times = {};   bat_click_power = {};
    for i=1:n_tags
        FC_temp = F_temp(F_temp.id==i,:);
        [~,valid_click_times] = count_spikes_AF_v0(Detected_Clicks.times,t,[FC_temp.smp1 FC_temp.smp2]);
        click_ids = ismember(Detected_Clicks.times,valid_click_times);
        bat_click_times{i,1} = Detected_Clicks.times(click_ids);
        bat_click_power{i,1} = Detected_Clicks.power(click_ids,:);
    end
    
    bat_click_ISI = cellfun(@diff,bat_click_times,'UniformOutput',false);
    bat_click_MeanPow = cellfun(@(x) mean(x,1),bat_click_power,'UniformOutput',false);
    
    figure('units','normalized','outerposition',[0.2 0.3 0.2 0.6]);
    tiledlayout(n_tags,2,'TileSpacing','tight');
    for i=1:n_tags
        nexttile;   histogram(bat_click_ISI{i,1},[0:0.005:0.3],'FaceColor',bat_clr(i,:));
        nexttile;   plot(bat_click_MeanPow{i,1},'Color',bat_clr(i,:));
    end
    
end


end