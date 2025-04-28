%% POOL REPLAYS AND QUANTIFY INTERESTING FEATURES

%=== Load Data
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
        'Dataset_2','14611','231208';...
        'Dataset_2','14611','231209';...
        'Dataset_2','14611','231210';...
        'Dataset_2','14611','231211';...
        'Dataset_2','14611','231212';...
        'Dataset_3','00000','240402';...
        'Dataset_3','00000','240403';...
        'Dataset_3','00000','240404';...
        'Dataset_3','00000','240405';...
        %'Dataset_3','00000','240406';...
        %'Dataset_3','00000','240407';...
        'Dataset_3','14543','240419';...
        'Dataset_3','14543','240420';...
        'Dataset_3','14543','240421';...
        %'Dataset_3','14543','240422';...
        'Dataset_3','29959','240402';...
        'Dataset_3','29959','240403';...
        'Dataset_3','29959','240404';...
        'Dataset_3','29959','240405';...
        %'Dataset_3','29959','240406';...
        %'Dataset_3','29959','240407';...
        ...
        };
    
    %=== Load data and aggregate them in the Multi_Day structure
    Folder = cd;    FileList = dir(fullfile(Folder, '**', 'Analyzed_NPs_*'));   SB = []; ids = [];
    for nc = 1:length(FileList)
        cd(FileList(nc).folder);
        load(FileList(nc).name);
        if ~isempty(SWP_table) && any(cellfun(@(row) isequal(row, SWP_table.unique_ID(1,:)), num2cell(sessions2include, 2)))
            SB = [SB; SWP_table];
        end
        disp([num2str(length(FileList)-nc),' remaining sessions to load...']);
    end
    cd(Folder); clearvars -except SB
end

%% Look at theta sweeps


%=== Params
max_rms = 2;                  % Max RMS
bin_size_1D = 0.15;
min_mtm = 0;                % Min Momentum
max_avs = inf;              % Max Average Posterior Spread
n_res_pp = 100;
min_common_bins = 17;

%=== Add a few features
[~,~,SB.groupID] =  unique([string(SB.unique_ID),string(SB.clus_id)],'rows');       % Id for each cluster from each session
SB.GoodClus = SB.avg_rms<max_rms;

%=== Select good events
SWB_table = SB(SB.GoodClus,:);
n_clus = size(SWB_table,1);

%=== Show average decoding error at Wmax
deATwmax=[];
wingbeatATwmax = [];
posteriorATwmax = [];
momentum = [];
posteriorATwmax_cut = [];
avspread = [];
sdATwmax = [];

for i=1:n_clus
    
    pp_matrix = SWB_table.posteriorATwmax{i,1};
    [~,max_tmp] = max(pp_matrix,[],1);
    
    avDevsweep = mean(abs(squeeze(max_tmp)-round(size(pp_matrix,1)/2)),1)*bin_size_1D;
    avspread = [avspread, mean(sqrt(squeeze(std(pp_matrix,[],1))),1)*bin_size_1D];
    momentum = [momentum avDevsweep];
    deATwmax_tmp =  SWB_table.deATwmax{i,1};
    deATwmax = [deATwmax, deATwmax_tmp];
    wingbeatATwmax = [wingbeatATwmax, SWB_table.wingbeatATwmax{i,1}];
    posteriorATwmax = cat(3,posteriorATwmax,imresize(pp_matrix,[n_res_pp size(pp_matrix,2)]));
    
    posteriorATwmax_cut = cat(3,posteriorATwmax_cut,pp_matrix(1:min_common_bins,:,:));
    
    sdATwmax = [sdATwmax, zscore(SWB_table.sdATwmax{i,1})];

    t =  SWB_table.wingbeattime{1,1};
    t_Hi = linspace(t(1),t(end),size(SWB_table.wingbeatATwmax{i,1},1));
    
end

%=== Plot average decoding error across the wave-cycle
avg_DE = mean(deATwmax(:,momentum>min_mtm & avspread<max_avs),2);
figure('units','normalized','outerposition',[.2 .2 .15 .3]);
tiledlayout(2,1,'TileSpacing','tight');
plot(t,avg_DE,'k'); hold on;    plot(t_Hi,normalize(mean(wingbeatATwmax,2),'range',[min(avg_DE) max(avg_DE)]),'r');
xlabel('Time (s)'); ylabel('Decoding error (m)');   legend('Decoding error','Wingbeat');    xlim('tight');  title(['Min Momentum: ', num2str(min_mtm,2)]);

%=== Plot average decoding error  
figure('units','normalized','outerposition',[.4 .1 0.1 .4]);
tiledlayout(2,1,'TileSpacing','tight');
data_tmp = deATwmax(:,momentum>min_mtm & avspread<max_avs)';
ax(1) = nexttile; plotWinterval_AF_v0(t,mean(data_tmp,'omitnan'),mean(data_tmp,'omitnan')-std(data_tmp,[],'omitnan')./sqrt(size(data_tmp,1)),mean(data_tmp,'omitnan')+std(data_tmp,[],'omitnan')./sqrt(size(data_tmp,1)),'r');
ylabel('Avg. Decoding Error (m)');   xticks([]); 
data_tmp = wingbeatATwmax';
ax(2) = nexttile; plotWinterval_AF_v0(t_Hi,mean(data_tmp,'omitnan'),mean(data_tmp,'omitnan')-std(data_tmp,[],'omitnan')./sqrt(size(data_tmp,1)),mean(data_tmp,'omitnan')+std(data_tmp,[],'omitnan')./sqrt(size(data_tmp,1)),'r');
ylabel('Avg. Acceleration');  
 xlim('tight');  xlabel('Time relative to wingbeat (s)');
linkaxes(ax,'x');

%=== Plot average spike density
figure('units','normalized','outerposition',[.6 .1 0.1 .4]);
tiledlayout(2,1,'TileSpacing','tight');
data_tmp = sdATwmax(:,momentum>min_mtm & avspread<max_avs)';
ax(1) = nexttile; plotWinterval_AF_v0(t,mean(data_tmp,'omitnan'),mean(data_tmp,'omitnan')-std(data_tmp,[],'omitnan')./sqrt(size(data_tmp,1)),mean(data_tmp,'omitnan')+std(data_tmp,[],'omitnan')./sqrt(size(data_tmp,1)),'r');
ylabel('Avg. Spike Density (m)');   xticks([]); 
data_tmp = wingbeatATwmax';
ax(2) = nexttile; plotWinterval_AF_v0(t_Hi,mean(data_tmp,'omitnan'),mean(data_tmp,'omitnan')-std(data_tmp,[],'omitnan')./sqrt(size(data_tmp,1)),mean(data_tmp,'omitnan')+std(data_tmp,[],'omitnan')./sqrt(size(data_tmp,1)),'r');
ylabel('Avg. Acceleration');  
 xlim('tight');  xlabel('Time relative to wingbeat (s)');
linkaxes(ax,'x');

%=== Wing-beat triggered decoding
figure('units','normalized','outerposition',[.5 .2 .1 .4]);
tiledlayout(4,1,'TileSpacing','tight');
medPostAtWmax = mean(posteriorATwmax(:,:,momentum>min_mtm & avspread<max_avs),3,'omitnan');
zx(1) = nexttile([3,1]);   imagesc([t(1) t(end)],[bin_size_1D*n_res_pp*0.5*[-1 1]],medPostAtWmax,prctile(medPostAtWmax, [1 100],'all')');    ylim([-2.5 2.5]);
colormap(flipud(gray));     title(['Min Momentum: ', num2str(min_mtm,2)]);
ylabel('Decoded m'); set(gca,'YDir','normal');  xticks([]); hold on;    refline(0,0); 
zx(2) = nexttile;   plot(t_Hi,mean(wingbeatATwmax,2)); xlabel('Time relative to wingbeat (s)');    ylim('tight');
linkaxes(zx,'x');   xlim('tight'); 

%=== Wing-beat triggered decoding
figure('units','normalized','outerposition',[.7 .2 .1 .4]);
tiledlayout(4,1,'TileSpacing','tight');
medPostAtWmax = mean(posteriorATwmax_cut(:,:,momentum>min_mtm & avspread<max_avs),3,'omitnan');
zx(1) = nexttile([3,1]);   imagesc([t(1) t(end)],[bin_size_1D*min_common_bins*0.5*[-1 1]],medPostAtWmax,prctile(medPostAtWmax, [0 100],'all')');    
colormap(flipud(gray));     title(['Min Momentum: ', num2str(min_mtm,2)]);
ylabel('Decoded m'); set(gca,'YDir','normal');  xticks([]); hold on;    refline(0,0); 
zx(2) = nexttile;   plot(t_Hi,mean(wingbeatATwmax,2)); xlabel('Time relative to wingbeat (s)');    ylim('tight');
linkaxes(zx,'x');   xlim('tight'); 

