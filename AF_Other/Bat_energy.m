%% === Load data and aggregate them in the Multi_Day structure
clr;    Folder = cd;    FileList = dir(fullfile(Folder, '**', 'Extracted_Behavior_*'));
if ~isempty(dir(fullfile(Folder,'Dataset_name.mat'))),load('Dataset_name.mat');disp(['DATSET: === ', Dataset_name, ' ===']);end
bat_dist_cell = cell(length(FileList),1);
for iii = 1:length(FileList)
    cd(FileList(iii).folder);
    load(FileList(iii).name);
    
    %=== Calculate bat distance
    n_pairs = length(bat_pairs);                                                                                    % Number of bat pairs
    bat_dist = zeros(T,n_pairs);
    for i = 1:n_pairs,    bat_dist(:,i) = vecnorm(r(:,:,bat_pairs(i,1))-r(:,:,bat_pairs(i,2)),2,2); end
    
    bat_dist_cell{iii,:} = bat_dist;
    disp([num2str(length(FileList)-iii),' remaining sessions to load...']);
end
cd(Folder); clearvars -except bat_dist_cell 

%% ===
n_sessions = size(bat_dist_cell,1);
T_all = cellfun(@(x) size(x,1),bat_dist_cell);
T_max = max(T_all);
T_min = min(T_all);
E_all = NaN(T_max,1);
E_all_min = zeros(T_min,1);

figure('units','normalized','outerposition',[.3 .3 .15 .4]);
nexttile;
for i = 1:n_sessions
    E = [];
    bat_dist = bat_dist_cell{i,1};
    E = sum(bat_dist,2).^1/2;
    E = sum(-1./bat_dist,2);
    E_all(1:numel(E),i) = E;
    E_all_min = smoothdata(E(1:T_min),'movmedian',1e4);
    plot(smoothdata(E,'movmedian',100),'--','LineWidth',.5,'Color',[.7 .7 .7]);   hold on;
end
plot(mean(E_all_min,2),'k','LineWidth',3);
xlabel('Sample into session');
ylabel('Energy');

