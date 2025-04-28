%% Analyze_Collective_Behavior_AF_v2
% Script for the analysis of collective behavior (April 2022)
% Requires Processed data from Ciholas recordings and Group name (D or F)

%=== Load data
extracted_BHVfile = dir(fullfile(cd, 'Extracted_Behavior_*'));          load(extracted_BHVfile.name);
disp(['Processing ', extracted_BHVfile.name, '...']);

%=== Parameters
n_ds = 100;                                                                                                     % Downsampling factor for clustering positions
n_min = 10*Fs;                                                                                                  % Minimum number of samples to be considered a cluster (at least 10s)
n_rep = 1000;                                                                                                   % Number of repetitions for the d_th shuffling
n_pairs = length(bat_pairs);                                                                                    % Number of bat pairs
edges_dist = 10.^linspace(-3,1,100);                                                                            % Edges for distance histograms (from 1 mm to 10 m)
r_lim = [-2.9 2.9; -2.6 2.6; 0 2.30];                                                                           % Room boundaries
edges_d = {r_lim(1,1):(r_lim(1,2)-r_lim(1,1))/20:r_lim(1,2) r_lim(2,1):(r_lim(2,2)-r_lim(2,1))/20:r_lim(2,2)};  % Edges for density histogram
bat_ids = [1:n_tags]';                                                                                          % Bat identities
d_th = 0.27;                                                                                                    % Threshold distance for various classifications  
for i = 1:n_tags; for j = 1:3; custom_map(:,j,i) = linspace(1,bat_clr(i,j))'; end; end                          % Custom graded colormap
if Group_name == 'D', r_fd = [-0.5,0.05,0.45];                                                                  % Feeder position  for D bats
else, r_fd = [2.77,0.82,1.75; 2.79,-0.99,1.64; 2.78,1.29,0.84; 2.78,-1.43,0.80]; end                            % Feeder positions for F bats
bat_nms_OBJ = ['Bt1'; 'Bt2'; 'Bt3'; 'Bt4'; 'OBJ';];                                                             % Modify Bat Names to include OBJ tag
bat_nms(1:n_tags-1,:) = bat_nms_OBJ(1:n_tags-1,:);  bat_nms(end,:) = bat_nms_OBJ(end,:);                        % Modify Bat Names to include OBJ tag
bat_pair_nms = [bat_nms(bat_pairs(:,1),:), '-'.*ones(length(bat_pairs),1), bat_nms(bat_pairs(:,2),:)];          % Bat Pairs Names
d = [.6 .9];                                                                                                    % Threshold distance for classifying close vs far {Optimal = [.6 .9]}

%=== Parameters
options.save_data = 1;                                                                                          % Save Data
options.clusterFl = 1;                                                                                          % Perform Flight Clustering
options.savefigures = 1;  fig_count = 1;                                                                        % Save Figures

%=== Create analysis folder for storing the results
if options.save_data
    analysis_directory=fullfile(pwd,['Analysis_',datestr(now, 'yymmdd_HHMM')]);
    if ~exist(analysis_directory,'dir')
        mkdir(analysis_directory);
    end
end

%% Define a corrective vector for excluding flight tails with v>v_th but not classified as 'flight'

b_corr = ones(size(bflying));
for i=1:n_tags-1
    if f_num(i)>0
        v_diff = diff(v_abs(:,i)>v_th);
        transitions = find(v_diff);
        if v_diff(transitions(1))   ==-1,b_corr(1:transitions(1),i) = 0;    end
        if v_diff(transitions(end)) == 1,b_corr(transitions(end):end,i) = 0;end
    end
end

%=== Define true rest vector
true_rest = all(bflying==0,2) & all(b_corr,2);

%% Calculate interbat-distances and percentage of time with 0,1,2,... bats flying

%=== bat_dist(:,i): distance between bat_pairs(i,:)
bat_dist = zeros(T,n_pairs);
for i = 1:n_pairs,    bat_dist(:,i) = vecnorm(r(:,:,bat_pairs(i,1))-r(:,:,bat_pairs(i,2)),2,2); end

%=== Bat distances at take-off and landing
f_nn = cell(n_tags,2);  f_bd = cell(n_tags,2);
for i=1:n_tags
    %=== All distances
    for j = 1:f_num(i)
        f_bd{i,1}(j,:)= pdist2(squeeze(r(f_smp{i,1}(j,1),:,:))',squeeze(r(f_smp{i,1}(j,1),:,i)));
        f_bd{i,2}(j,:)= pdist2(squeeze(r(f_smp{i,2}(j,1),:,:))',squeeze(r(f_smp{i,2}(j,1),:,i)));
    end
    %=== NN distances
    f_nn{i,1} = min(bat_dist(f_smp{i,1},any(bat_pairs==i,2)),[],2);
    f_nn{i,2} = min(bat_dist(f_smp{i,2},any(bat_pairs==i,2)),[],2);
end

%=== Quantify percentages of time with 0,1,2,3...bats flying
perc_time = zeros(n_tags+1,1);
for i=1:n_tags+1
    perc_time(i,:) = nnz(sum(bflying,2)== i-1)/T;
end

%% Cluster position data

%=== Pool positions from all the bats and cluster them
X = []; for i = 1:n_tags-1; X = [X;r(~bflying(:,i)& b_corr(:,i),:,i)]; end
figure('units','normalized','outerposition',[0.3 0.3 0.6 0.4]);
tiledlayout(1,4,'TileSpacing','tight');
[id_all,centroid_all] = Cluster3D_AF_v1(X,n_min,0.2,100);
fig_count = saveFig(analysis_directory,batdate,fig_count,options.savefigures);

%=== Correct feeder position(s) 
if Group_name == 'D'
    [~,fd_idx] = min(vecnorm(r_fd-centroid_all,2,2));
    r_fd(1,:) = centroid_all(fd_idx,:);                                                  
else
    for i=1:4
    [feeder_distance,fd_idx] = min(vecnorm(r_fd(i,:)-centroid_all,2,2));
    if feeder_distance<0.2,r_fd(i,:) =  centroid_all(fd_idx,:);end  % Do not correct if further than 20 cm
    end
end

%=== Single bat clusters
r_clus_id = NaN(T,n_tags);  centroid = cell(n_tags,1);  samplesIn = cell(n_tags,1); p_val_clus = zeros(n_tags,1); perc_clus = zeros(n_tags,1);  CH_idx = zeros(n_tags,19);
figure('units','normalized','outerposition',[0.3 0 0.33 1]);
tiledlayout(n_tags,4,'TileSpacing','tight');
for i = 1:n_tags
    X = r(~bflying(:,i)& b_corr(:,i),:,i);                                                              % Consider only stationary epochs
    [id,centroid{i,1},samplesIn{i,1},p_val_clus(i,:),CH_idx(i,:)] = Cluster3D_AF_v1(X,n_min,0.2,n_ds);  % Perform clustering
    perc_clus(i,:) = nnz(id)/length(id);                                                                % fraction of points belonging to a cluster
    r_clus_id(~bflying(:,i)& b_corr(:,i),i)=id;                                                         % Assign cluster id across entire session (nan when bat is flying)
end
fig_count = saveFig(analysis_directory,batdate,fig_count,options.savefigures);

%% Analysis of context (possible states) and dimensionality

%=== Define state matrix as a T x n_tags vector 
%=== A state at time t is defined as the set of positional clusters for each bat
%=== When a bat is flying (or not in a cluster), its state is NaN
state= NaN(T,n_tags-1);
id_all_temp = id_all;
for i = 1:n_tags-1 
    valid_samples = nnz(~bflying(:,i)& b_corr(:,i));
    state(~bflying(:,i)& b_corr(:,i),i) = id_all_temp(1:valid_samples); 
    id_all_temp = id_all_temp(valid_samples+1:end);
end

%=== Define context at rest (only clustered positions) 
state_rest = state(~any(isnan(state),2) & ~any(state == 0,2),:);
[state_templates,state_reference,state_categories] = unique(state_rest,'rows');
state_density = zeros(size(state_templates,1),1);
for i = 1:size(state_templates,1)
    state_density(i)  = nnz(all(state_rest == state_templates(i,:),2))/size(state_rest,1);
end

%===State matrix 
for s = 2:max(id_all)
    state_matrix(s-1,:)= sum(state_rest==s,1)/size(state_rest,1);
end

%=== Calculate the entropy(observed and maximal)
E_obs = -sum(state_density.*log2(state_density));
E_max = log2(size(state_templates,1));

%% Low dimensional representation of states

r_emb = [];
% %=== Look at epochs when all bats are resting
% r_rest = r(true_rest,:,1:n_tags-1);
% 
% %=== Define positional samples to look at (t-SNE works better if they are independent (roughly one sample in between take offs), 3s DEFAULT)
% r_emb = state_embed_AF_v0(r_rest,3*Fs,'sammon',2,0);
% fig_count = saveFig(analysis_directory,batdate,fig_count,options.savefigures);
% 
% %=== Look at preferred positions of the embedded points
% sts_2look = state(true_rest,:);
% best_num = 6;                   % number of best states to look at
% [unique_states,~,~] = unique(sts_2look,'rows');
% for i = 1:size(unique_states,1)
%     state_occ(i) = nnz(all(sts_2look == unique_states(i,:),2));
% end
% [state_occ,sorted_ids] = sort(state_occ,'descend');
% common_states = unique_states(sorted_ids(1:best_num),:);
% 
% %=== Redefine state id for the most common states
% state_id = zeros(size(sts_2look,1),1);
% for i=1:best_num
%     state_id(all(sts_2look == common_states(i,:),2)) = i; 
% end
% 
% %=== Show t-sne embedding and meaning of some of the clustered states
% figure('units','normalized','outerposition',[0.1 0.3 0.4 0.4]);
% tiledlayout(2,5,'TileSpacing','compact');
% nexttile([2,2]);   gscatter(r_emb(:,1),r_emb(:,2),state_id,[0 0 0; jet(best_num)],[],[3; 10*ones(best_num,1)]);    axis off; sgtitle('Most populated clusters');
% tlordr = [3 4 5 8 9 10];
% for j = 1:best_num
%     nexttile(tlordr(j));    
%     for i = 1:n_tags-1
%     scatter3(r_rest(state_id == j,1,i),r_rest(state_id == j,2,i),r_rest(state_id == j,3,i),10,bat_clr(i,:),'filled');
%         hold on;
%     end
%     hold off;   xlim(r_lim(1,:)); ylim(r_lim(2,:)); zlim(r_lim(3,:));   grid on;    
%     set(gca,'xticklabel',[]);   set(gca,'yticklabel',[]);   set(gca,'zticklabel',[]); 
%     title([num2str(j), ' Most common state']);
% end
% fig_count = saveFig(analysis_directory,batdate,fig_count,options.savefigures);

%% Permutation-invariant embedding
%=== Takes time!!!
while 0
    r_emb = state_embed_AF_v0(r_rest,3*Fs,'sammon',2,1);
    fig_count = saveFig(analysis_directory,batdate,fig_count,options.savefigures);
    
    %=== Cluster
    perm_id = dbscan(r_emb,0.5,10);
    gscatter(r_emb(:,1),r_emb(:,2),perm_id);
    
    %=== Calculate occupancy
    [G_c,G_id] = groupcounts(perm_id);
    [G_c,st_i] = sort(G_c,'descend');
    G_id = G_id(st_i);
    
    %=== Redefine state id for the most common states
    best_num = 3;
    state_id = zeros(size(perm_id,1),1);
    for i=1:best_num
        state_id(perm_id == G_id(i)) = i;
    end
    
    %=== Show t-sne embedding and meaning of some of the clustered states
    figure('units','normalized','outerposition',[0.1 0.3 0.7 0.52]);
    tiledlayout(3,11,'TileSpacing','tight');
    nexttile(1,[3,3]);   scatter(r_emb(:,1),r_emb(:,2),5,'k','filled');                                                  axis off; sgtitle('2-D Embedding (Identity-Invariant)');
    nexttile(4,[3,3]);   gscatter(r_emb(:,1),r_emb(:,2),state_id,[0 0 0; jet(best_num)],[],[3; 10*ones(best_num,1)]);    axis off; sgtitle('Most populated clusters');
    tlordr = [7:11, 18:22, 29:33];
    for j = 1:15
        nexttile(tlordr(j));
        lookat = r_rest(state_id == ceil(j/5),:,:);
        n = randi(size(lookat,1));
        for i = 1:n_tags,plot3([0; lookat(n,1,i)],[0; lookat(n,2,i)],[0; lookat(n,3,i)],'Color',bat_clr(i,:),'LineWidth',3); hold on; end
        scatter3(lookat(n,1,:),lookat(n,2,:),lookat(n,3,:),40,bat_clr,'filled');    hold off;
        xlim(r_lim(1,:)); ylim(r_lim(2,:)); zlim(r_lim(3,:));   view(0,90); title([num2str(ceil(j/5)),' most common state']);
    end
    fig_count = saveFig(analysis_directory,batdate,fig_count,options.savefigures);
end

%% Association Index (see Niehof & Morley 2012. Determining the significance of associations between two series of discrete events: bootstrap methods)

h = 3*Fs;   % window size
L = 1000;   % number of time lags
lag = 20;   % maximal lag
u = linspace(-lag*Fs,lag*Fs,L);

figure('units','normalized','outerposition',[0.3 0 0.1 1]);
tiledlayout(ceil(n_pairs/2),2,'TileSpacing','tight');
c=zeros(n_pairs,L);
for i = 1:n_pairs
    for j = 1:L
        for a = f_smp{bat_pairs(i,1),1}'
            c(i,j) = c(i,j)+nnz(abs(f_smp{bat_pairs(i,2),1}-a+u(j))<h/2);
        end
    end
    nexttile;   plot(linspace(-lag,lag,L),c(i,:));    
    grid on;
    xlabel('Time lag (s)'); ylabel('AI'); 
    title([bat_nms(bat_pairs(i,1),:) '-' bat_nms(bat_pairs(i,2),:)]);
end

fig_count = saveFig(analysis_directory,batdate,fig_count,options.savefigures);

%% Cluster flights

% alpha_clus = 0.3; 'Points',7  no Frechet
% alpha_clus = 1.1; 'Points',9  Frechet

alpha_clus = 0.6;         %Parameter for flight clustering (historically 1.1)
f_cls = cell(n_tags,1);     
for i = 1:n_tags
    if options.clusterFl && i~= n_tags
        f_clus(i) = FlightClus_AF_v2(squeeze(r(:,:,i)),bflying(:,i),'Alpha',alpha_clus,'Frechet',1,'Points',7);
        f_cls{i,1} = f_clus(i).sorted_ids;
        fig_count = saveFig(analysis_directory,batdate,fig_count,options.savefigures);
    else   
        f_clus(i) = f_clus(i-1);
        f_cls{i,1} = 0.*f_smp{i,1}';  
    end
end

%% Create FLIGHTS table features associated with flights

all_tko = t(vertcat(f_smp{:,1}));
all_lnd = t(vertcat(f_smp{:,2}));

warning('off','all');
c = 1;  FLIGHTS = table();
for i = 1:n_tags
    for j=1:f_num(i)
        FLIGHTS.id(c,:) = i;                                                        % id of the bat flying
        FLIGHTS.smp1(c,:) = f_smp{i,1}(j,1);                                        % sample takeoff
        FLIGHTS.smp2(c,:) = f_smp{i,2}(j,1);                                        % sample landing
        FLIGHTS.t1(c,:) = t(f_smp{i,1}(j,1));                                       % time takeoff
        FLIGHTS.t2(c,:) = t(f_smp{i,2}(j,1));                                       % time landing
        FLIGHTS.maxN(c,:) = max(sum(bflying(f_smp{i,1}(j,1):f_smp{i,2}(j,1),:),2)); % max number of bats flying during that flight
        FLIGHTS.r1(c,:) = r(f_smp{i,1}(j,1),:,i);                                   % position at takeoff
        FLIGHTS.r2(c,:) = r(f_smp{i,2}(j,1),:,i);                                   % position at landing
        FLIGHTS.pclus1(c,:) = r_clus_id(f_smp{i,1}(j,1)-1,i);                       % positional cluster at take off
        FLIGHTS.pclus2(c,:) = r_clus_id(f_smp{i,2}(j,1)+1,i);                       % positional cluster at landing
        FLIGHTS.fclus(c,:) = f_cls{i,1}(1,j);                                       % flight cluster
        FLIGHTS.state1(c,:) = {state(f_smp{i,1}(j,1)-1,:)};                         % State 1 sample before take off
        FLIGHTS.state2(c,:) = {state(f_smp{i,2}(j,1)+1,:)};                         % State 1 sample after landing
        
        FLIGHTS.odst1(c,:) = squeeze(vecnorm(r(f_smp{i,1}(j,1),:,i)-r(f_smp{i,1}(j,1),:,bat_ids(bat_ids~=i))))';    % Distance from the other bats or object at takeoff (ordered: bat1, bat2, object)
        FLIGHTS.odst2(c,:) = squeeze(vecnorm(r(f_smp{i,2}(j,1),:,i)-r(f_smp{i,2}(j,1),:,bat_ids(bat_ids~=i))))';    % Distance from the other bats or object at landing (ordered: bat1, bat2, object)
        FLIGHTS.nno_dist(c,:) = vecnorm(r(f_smp{i,2}(j,1),:,i)-r(f_smp{i,2}(j,1),:,end));                           % Distance from the object at landing
        FLIGHTS.nnb_dist(c,:) = min(vecnorm(r(f_smp{i,2}(j,1),:,i)-r(f_smp{i,2}(j,1),:,bat_ids(bat_ids~=i & bat_ids~=n_tags)),2,2));  % Nn bat distance
        FLIGHTS.pos(c,:) = {r(f_smp{i,2}(j,1)+[-Fs:0.5*Fs],:,i)};                                                   % Position around landing [-1s,0.5s]
        FLIGHTS.hdn(c,:) = {angle(f_smp{i,2}(j,1)+[-Fs:0.5*Fs],i)};                                                 % Heading angle around landing [-1s,0.5s]
        FLIGHTS.acc(c,:) = {a_abs(f_smp{i,2}(j,1)+[-Fs:0.5*Fs],i)};                                                 % Accelerometer around landing [-1s,0.5s]
        
        cond_f = vecnorm(r(f_smp{i,2}(j,1),:,i)-r_fd(1,:),2,2)<d_th;
        cond_s = any(vecnorm(r(f_smp{i,2}(j,1),:,i)-r(f_smp{i,2}(j,1),:,bat_ids(bat_ids~=i)),2,2)<d_th);
        if      cond_f && ~cond_s
            FLIGHTS.class(c,:) = 'f';                                               % exclusive feeding flight
        elseif ~cond_f && cond_s
            FLIGHTS.class(c,:) = 's';                                               % exclusive social flight
        elseif  cond_f && cond_s
            FLIGHTS.class(c,:) = 'b';                                               % feeding + social flight
        else
            FLIGHTS.class(c,:) = 'o';                                               % other flight
        end
        
        %=== Calculate time to previous landing and time to next takeoff
        if ~isempty(min(t(f_smp{i,1}(j,1))-all_lnd(all_lnd<t(f_smp{i,1}(j,1)))))
            FLIGHTS.t_prev_lnd(c,:) = min(t(f_smp{i,1}(j,1))-all_lnd(all_lnd<t(f_smp{i,1}(j,1))));
        else
            FLIGHTS.t_prev_lnd(c,:) = Inf;
        end
        if ~isempty(min(-t(f_smp{i,2}(j,1))+all_tko(all_tko>t(f_smp{i,2}(j,1)))))
            FLIGHTS.t_next_tko(c,:) = min(-t(f_smp{i,2}(j,1))+all_tko(all_tko>t(f_smp{i,2}(j,1))));
        else
            FLIGHTS.t_next_tko(c,:) = Inf;
        end
        c = c+1;
    end
end
warning('on','all');

%=== Define configurations at landing
FLIGHTS.config = FLIGHTS.t1*nan;
tmp_stt = discretize(FLIGHTS.odst2,[0 d inf]);                            % 1 present, 3 absent, 2 undetermined
tmp_stt(tmp_stt==3) = 0;    tmp_stt(tmp_stt==2) = nan;                    % 1 present, 0 absent, Nan undetermined
for i=0:7
    config_tmp = str2num(dec2bin(i,3)')';
    FLIGHTS.config(all(tmp_stt == config_tmp,2)) = i+1;
end
% [0, 0, 0],2)) = 1;                          % Landing on empty spot
% [0, 0, 1],2)) = 2;                          % Landing on object
% [0, 1, 0],2)) = 3;                          % Landing on 2nd bat
% [0, 1, 1],2)) = 4;                          % Landing on 2nd bat + object
% [1, 0, 0],2)) = 5;                          % Landing on 1st bat
% [1, 0, 1],2)) = 6;                          % Landing on 1st bat + object
% [1, 1, 0],2)) = 7;                          % Landing on both bats
% [1, 1, 1],2)) = 8;                          % Landing on both bats + object

%% Look at the social nonsocial landings and at the obejct

%=== Look at the distances from the object on the perch and feeder sides
figure('units','normalized','outerposition',[0.4 0.3 0.3 0.3]);
tiledlayout(2,4,'TileSpacing','tight');
FC_temp = FLIGHTS(vecnorm(FLIGHTS.r2 - [-2.5 2.2 2.0],2,2)<1,:);
nexttile;   scatter3(FC_temp.r2(:,1),FC_temp.r2(:,2),FC_temp.r2(:,3),5,'filled','MarkerFaceColor','k');
axis equal; xlim(r_lim(1,:)); ylim(r_lim(2,:)); zlim(r_lim(3,:));   view(0,90);
set(gca,'XTickLabels',[],'YTickLabels',[],'ZTickLabels',[]);
nexttile([1,3]);   bar(FC_temp.nno_dist);  hold on;    plot(xlim,[.6 .6],'r--');      ylabel('Distance (m)');
FC_temp = FLIGHTS(vecnorm(FLIGHTS.r2 - [+2.7 0.9 1.7],2,2)<1,:);
nexttile;   scatter3(FC_temp.r2(:,1),FC_temp.r2(:,2),FC_temp.r2(:,3),5,'filled','MarkerFaceColor','k');
axis equal; xlim(r_lim(1,:)); ylim(r_lim(2,:)); zlim(r_lim(3,:));   view(0,90);
set(gca,'XTickLabels',[],'YTickLabels',[],'ZTickLabels',[]);
nexttile([1,3]);   bar(FC_temp.nno_dist);  hold on;    plot(xlim,[.6 .6],'r--');   ylabel('Distance (m)');
sgtitle('Landing Close to object');
fig_count = saveFig(analysis_directory,batdate,fig_count,options.savefigures); 

%=== Look at trajectories depending on the landing configuration (feeder only)
for i=1:n_tags-1
    
    figure('units','normalized','outerposition',[0.2 0.3 0.43 0.33]);
    tiledlayout(2,7,'TileSpacing','tight');
    
    % Select flights of a specific bat to the feeder
    FC_temp = FLIGHTS(FLIGHTS.id == i & any(FLIGHTS.class == ['b','f'],2),:);   
    FCB_temp = sortrows(FC_temp,'nnb_dist','descend');    % Sort by landing distance from bat
    FCO_temp = sortrows(FC_temp,'nno_dist','descend');    % Sort by landing distance from object
    
    %=== Plots (bat)
    ds_temp = FCB_temp.nnb_dist;    id_temp = discretize(FCB_temp.nnb_dist,[0 d inf]);  cls_tmp = knnsearch(FCB_temp.nnb_dist,d(1));    far_tmp = knnsearch(FCB_temp.nnb_dist,d(2));
    if nnz(id_temp==1) && nnz(id_temp==3)
        nexttile;   plot(ds_temp ,-0.5+[1:size(FCB_temp,1)],'Color','k','LineWidth',2);   ylim([0 size(FCB_temp,1)]);
        set(gca, 'YDir','reverse','YGrid','on','YTickLabel',[]);    title('To Bat');    xlim([0 7]);    xlabel('Distance at landing (m)');
        nexttile;   hold on;
        cellfun(@(x) plot3(x(:,1),x(:,2),x(:,3),'k'),FCB_temp.pos(id_temp==1,:));   % Close
        cellfun(@(x) plot3(x(:,1),x(:,2),x(:,3),'r'),FCB_temp.pos(id_temp==3,:));   % Far
        axis equal; xlim(r_lim(1,:)); ylim(r_lim(2,:)); zlim(r_lim(3,:));   view(0,90); xticks([]); yticks([]);
        nexttile;   imagesc([-1, 0.5],[],cell2mat(cellfun(@(x) x(:,1)',FCB_temp.pos,'UniformOutput',false)));   hold on;
        plot([0 0],ylim,'k--'); plot(xlim,cls_tmp*[1 1],'r'); plot(xlim,far_tmp*[1 1],'k'); hold off;   title('x (m)'); xlabel('Time (s)');
        nexttile;   imagesc([-1, 0.5],[],cell2mat(cellfun(@(x) x(:,2)',FCB_temp.pos,'UniformOutput',false)));   hold on;
        plot([0 0],ylim,'k--'); plot(xlim,cls_tmp*[1 1],'r'); plot(xlim,far_tmp*[1 1],'k'); hold off;   title('y (m)'); xlabel('Time (s)');
        nexttile;   imagesc([-1, 0.5],[],cell2mat(cellfun(@(x) x(:,3)',FCB_temp.pos,'UniformOutput',false)));   hold on;
        plot([0 0],ylim,'k--'); plot(xlim,cls_tmp*[1 1],'r'); plot(xlim,far_tmp*[1 1],'k'); hold off;   title('z (m)'); xlabel('Time (s)');
        nexttile;   imagesc([-1, 0.5],[],cell2mat(cellfun(@(x) x(:,1)',FCB_temp.acc,'UniformOutput',false)));   hold on;
        plot([0 0],ylim,'k--'); plot(xlim,cls_tmp*[1 1],'r'); plot(xlim,far_tmp*[1 1],'k'); hold off;   title('a (g)'); xlabel('Time (s)');
        nexttile;   hold on;
        cellfun(@(x) plot([1:numel(x)]./Fs-1,x,'Color',[.1 .1 .1 .1]),FCB_temp.hdn(id_temp==1,:));   % Close
        cellfun(@(x) plot([1:numel(x)]./Fs-1,x,'Color',[.9 .1 .1 .1]),FCB_temp.hdn(id_temp==3,:));   % Far
        plot([1:size(FCB_temp.hdn{1},1)]./Fs-1,median([FCB_temp.hdn{id_temp==1,:}],2),'k','LineWidth',2)
        plot([1:size(FCB_temp.hdn{1},1)]./Fs-1,median([FCB_temp.hdn{id_temp==3,:}],2),'r','LineWidth',2)
        xlim([-1 0]);   ylim([5.6 6.4]);    title('Angle'); xlabel('Time (s)');
    else
        nexttile;   nexttile;   nexttile;   nexttile;
        nexttile;   nexttile;   nexttile;
    end
    
    %=== Plots (object)
    ds_temp = FCO_temp.nno_dist;    id_temp = discretize(FCO_temp.nno_dist,[0 d inf]);  cls_tmp = knnsearch(FCO_temp.nno_dist,d(1));    far_tmp = knnsearch(FCO_temp.nno_dist,d(2));
    if nnz(id_temp==1) && nnz(id_temp==3)
        nexttile;   plot(ds_temp ,-0.5+[1:size(FCO_temp,1)],'Color','k','LineWidth',2);   ylim([0 size(FCO_temp,1)]);
        set(gca, 'YDir','reverse','YGrid','on','YTickLabel',[]);    title('To Object');    xlim([0 7]);    xlabel('Distance at landing (m)');
        nexttile;   hold on;
        cellfun(@(x) plot3(x(:,1),x(:,2),x(:,3),'k'),FCO_temp.pos(id_temp==1,:));   % Close
        cellfun(@(x) plot3(x(:,1),x(:,2),x(:,3),'b'),FCO_temp.pos(id_temp==3,:));   % Far
        axis equal; xlim(r_lim(1,:)); ylim(r_lim(2,:)); zlim(r_lim(3,:));   view(0,90); xticks([]); yticks([]);
        nexttile;   imagesc([-1, 0.5],[],cell2mat(cellfun(@(x) x(:,1)',FCO_temp.pos,'UniformOutput',false)));   hold on;
        plot([0 0],ylim,'k--'); plot(xlim,cls_tmp*[1 1],'r'); plot(xlim,far_tmp*[1 1],'k'); hold off;   title('x (m)'); xlabel('Time (s)');
        nexttile;   imagesc([-1, 0.5],[],cell2mat(cellfun(@(x) x(:,2)',FCO_temp.pos,'UniformOutput',false)));   hold on;
        plot([0 0],ylim,'k--'); plot(xlim,cls_tmp*[1 1],'r'); plot(xlim,far_tmp*[1 1],'k'); hold off;   title('y (m)'); xlabel('Time (s)');
        nexttile;   imagesc([-1, 0.5],[],cell2mat(cellfun(@(x) x(:,3)',FCO_temp.pos,'UniformOutput',false)));   hold on;
        plot([0 0],ylim,'k--'); plot(xlim,cls_tmp*[1 1],'r'); plot(xlim,far_tmp*[1 1],'k'); hold off;   title('z (m)'); xlabel('Time (s)');
        nexttile;   imagesc([-1, 0.5],[],cell2mat(cellfun(@(x) x(:,1)',FCO_temp.acc,'UniformOutput',false)));   hold on;
        plot([0 0],ylim,'k--'); plot(xlim,cls_tmp*[1 1],'r'); plot(xlim,far_tmp*[1 1],'k'); hold off;   title('a (g)'); xlabel('Time (s)');
        nexttile;   hold on;
        cellfun(@(x) plot([1:numel(x)]./Fs-1,x,'Color',[.1 .1 .1 .1]),FCO_temp.hdn(id_temp==1,:));   % Close
        cellfun(@(x) plot([1:numel(x)]./Fs-1,x,'Color',[.1 .1 .9 .1]),FCO_temp.hdn(id_temp==3,:));   % Far
        plot([1:size(FCO_temp.hdn{1},1)]./Fs-1,median([FCO_temp.hdn{id_temp==1,:}],2),'k','LineWidth',2)
        plot([1:size(FCO_temp.hdn{1},1)]./Fs-1,median([FCO_temp.hdn{id_temp==3,:}],2),'b','LineWidth',2)
        xlim([-1 0]);   ylim([5.6 6.4]);    title('Angle'); xlabel('Time (s)');
    else
        nexttile;   nexttile;   nexttile;   nexttile;
        nexttile;   nexttile;   nexttile;
    end
    
    sgtitle(bat_nms_OBJ(i,:));
    fig_count = saveFig(analysis_directory,batdate,fig_count,options.savefigures);  
    
end

%=== Loop across each bat and each positional cluster at landing
for i=1:n_tags-1
        
    FC_temp = FLIGHTS(FLIGHTS.id==i,:);
    nbcls = max(FC_temp.pclus2)+1;
    figure('units','normalized','outerposition',[0.4 0.15 0.15 min([0.1*nbcls,0.8])]);
    tiledlayout(nbcls,3,'TileSpacing','tight');
    for pp =0:nbcls-1
        FP_temp = FC_temp(FC_temp.pclus2==pp,:);
        nexttile;   hold on;
        cellfun(@(x) plot3(x(:,1),x(:,2),x(:,3),'k'),FP_temp.pos(FP_temp.nnb_dist>0.9,:));
        cellfun(@(x) plot3(x(:,1),x(:,2),x(:,3),'r'),FP_temp.pos(FP_temp.nnb_dist<0.6,:));
        axis equal; xlim(r_lim(1,:)); ylim(r_lim(2,:)); zlim(r_lim(3,:));   view(0,90); xticks([]); yticks([]);
        title(['NS = ',num2str(nnz(FP_temp.nnb_dist>0.9)),', S = ', num2str(nnz(FP_temp.nnb_dist<0.6))]);
        
        nexttile;   hold on;
        cellfun(@(x) plot3(x(:,1),x(:,2),x(:,3),'b'),FP_temp.pos(FP_temp.nno_dist>0.9,:));
        cellfun(@(x) plot3(x(:,1),x(:,2),x(:,3),'g'),FP_temp.pos(FP_temp.nno_dist<0.6,:));
        axis equal; xlim(r_lim(1,:)); ylim(r_lim(2,:)); zlim(r_lim(3,:));   view(0,90); xticks([]); yticks([]);
        title(['NO = ',num2str(nnz(FP_temp.nno_dist>0.9)),', O = ', num2str(nnz(FP_temp.nno_dist<0.6))]);
        
        nexttile;   histogram(FP_temp.config,-0.5+[1:9]);
    end
    sgtitle(bat_nms_OBJ(i,:));
    fig_count = saveFig(analysis_directory,batdate,fig_count,options.savefigures);   
end


%% Relationships between bats: Proximity, Landing and Coupling Index

%=== Initialization of matrices
%=== Proximity and coupling indexes and associated p-values
PI = zeros(n_rep+1,n_pairs);    PI_p = zeros(n_pairs,1);    PI_c = zeros(n_pairs,1);
CI = zeros(n_rep+1,n_pairs);    CI_p = zeros(n_pairs,1);    CI_c = zeros(n_pairs,1);
%=== Landing index (forward and inverse) and associated p-values
LI_1 = zeros(n_rep+1,n_pairs);  LI_2 = zeros(n_rep+1,n_pairs); LI_p = zeros(n_pairs,2);
%=== Shuffling procedure
for i = 1:n_pairs
    %=== Observed values
    PI(1,i) = nnz(bat_dist(:,i)<d_th)/T;                                    % PI: fraction of frames < d_th
    CI(1,i) = sum(bflying(:,bat_pairs(i,1)).*bflying(:,bat_pairs(i,2)),1)/T;% CI: fraction of frames the bats fly together
    landing_temp = f_bd{bat_pairs(i,1),2}<d_th;                             % Temporary matrix with all the landings<d_th for bat_pairs(i,1)
    LI_1(1,i) = nnz(landing_temp(:,bat_pairs(i,2)))/f_num(bat_pairs(i,1));  % LI(1,i,1) = fraction of flights for bat_pairs(i,1) landing less than d_th from bat_pairs(i,2)
    landing_temp = f_bd{bat_pairs(i,2),2}<d_th;                             % Temporary matrix with all the landings<d_th for bat_pairs(i,2)
    LI_2(1,i) = nnz(landing_temp(:,bat_pairs(i,1)))/f_num(bat_pairs(i,2));  % LI(1,i,2) = fraction of flights for bat_pairs(i,2) landing less than d_th from bat_pairs(i,1)
    frames_to_shift = randi([60*Fs T-60*Fs],1,n_rep);                       % Random number of frames between 1 min and session time-1min
    %=== Values after shuffling position and velocity
    parfor n = 1:n_rep
        r_shfl = circshift(r(:,:,bat_pairs(i,2)),frames_to_shift(n),1);     % Shift the position of the second bat in the couple
        PI(n+1,i) = nnz(vecnorm(r(:,:,bat_pairs(i,1))-r_shfl,2,2)<d_th)/T;  % Recalculate the PI
        LI_1(n+1,i) = nnz(vecnorm(r(f_smp{bat_pairs(i,1),2},:,bat_pairs(i,1))-r_shfl(f_smp{bat_pairs(i,1),2},:),2,2)<d_th)/f_num(bat_pairs(i,1)); % Calculate LI as the fraction of flights for bat_pairs(i,1) landing less than d_th from r_shfl
        r_shfl = circshift(r(:,:,bat_pairs(i,1)),frames_to_shift(n),1);     % Shift the position of the first bat in the couple
        LI_2(n+1,i) = nnz(vecnorm(r(f_smp{bat_pairs(i,2),2},:,bat_pairs(i,2))-r_shfl(f_smp{bat_pairs(i,2),2},:),2,2)<d_th)/f_num(bat_pairs(i,2)); % Calculate LI as the fraction of flights for bat_pairs(i,2) landing less than d_th from r_shfl
        b_shfl = circshift(bflying(:,bat_pairs(i,2)),frames_to_shift(n),1); % Shift the flight vector of the second bat in the couple
        CI(n+1,i) = sum(bflying(:,bat_pairs(i,1)).*b_shfl,1)/T;             % Recalculate the CI
    end
    PI_p(i,:) = nnz(PI(2:end,i)>PI(1,i))/n_rep;                                 % Calculate the p value as the fraction of PI > observed PI
    PI_c(i,:) = PI(1,i)-mean(PI(2:end,i));                                      % Calculate the corrected proximity index (by subtracting the mean of the shuffled
    CI_p(i,:) = nnz(CI(2:end,i)>CI(1,i))/n_rep;                                 % Calculate the p value as the fraction of CI > observed CI
    CI_c(i,:) = CI(1,i)-mean(CI(2:end,i));
    LI_p(i,1) = nnz(LI_1(2:end,i)>LI_1(1,i))/n_rep;                             % Calculate the p value as the fraction of LI > observed LI
    LI_p(i,2) = nnz(LI_2(2:end,i)>LI_2(1,i))/n_rep;                             % Calculate the p value as the fraction of LI > observed LI
end

LI = cat(3,LI_1,LI_2);

%% Relationships between bats: Preferred Locations

%=== Extract first 3 preferred locations for each bat
pref_location = NaN(n_tags,3,3);    pref_time_loc = NaN(n_tags,1,3);
for i=1:n_tags
    [~,indx] = sort(samplesIn{i,1},'descend','MissingPlacement','last');                     % Sort centroids occupancy
    num_pref = min(3,numel(indx));                                                           % Check how many centroids survive
    pref_location(i,:,1:num_pref) = centroid{i,1}(indx(1:num_pref),:)';                      % Assign the first 3 preferred locations
    pref_time_loc(i,:,1:num_pref) = samplesIn{i,1}(indx(1:num_pref),:)'./T;                  % Fraction of time in the preferred locations
end

%% Relationships between bats: Coordinated flight

%=== Probability of flight (one bat vs two bats)
p_fly1 = sum(bflying,1)./T; p_fly2 = zeros(n_pairs,2);
for i = 1:length(bat_pairs)
    p_fly2(i,1) = sum(bflying(:,bat_pairs(i,1)).*bflying(:,bat_pairs(i,2)),1)/T;
    p_fly2(i,2) = p_fly1(bat_pairs(i,1))*p_fly1(bat_pairs(i,2));
end

%=== Calculate the difference in the number of flights across the session
figure('units','normalized','outerposition',[0.3 0.1 0.2 0.9]);
tiledlayout(1,1,'TileSpacing','tight');
cum_flights = [zeros(1,n_tags); cumsum(max(diff(bflying),0))];
fligh_diff = zeros(T,n_pairs);
nexttile;
for i = 1:n_pairs
    fligh_diff(:,i) = cum_flights(:,bat_pairs(i,1))-cum_flights(:,bat_pairs(i,2));
    plot(t,fligh_diff(:,i),'LineWidth',2);    hold on;
    text(t(T),fligh_diff(T,i),[bat_nms(bat_pairs(i,1),:) '-' bat_nms(bat_pairs(i,2),:)]);
end
rectangle('Position',[0 -10 t(T) 20],'FaceColor',[0 0 0 0.2],'EdgeColor','none');
rectangle('Position',[0 -20 t(T) 40],'FaceColor',[0 0 0 0.1],'EdgeColor','none');
hold off;
xlabel('Time (s)'); ylabel('Difference in the number of flights');
fig_count = saveFig(analysis_directory,batdate,fig_count,options.savefigures);

%% Analysis of the accelerometer Signal

figure('units','normalized','outerposition',[0.5 0 0.1 1]);
tiledlayout(n_tags,1,'TileSpacing','none');
for i=1:n_tags
    %=== Calculate FFT and PSD of the accelerometer signal
    extended_flight = movmax(bflying(:,i),[1*Fs 1*Fs]);
    trace = a_abs(find(extended_flight),i);
    n = 2^nextpow2(numel(trace));
    Y = fft(trace,n);
    P = abs(Y/n).^2;
    f = Fs*(0:n/2)/n;
    PSD = P(1:n/2+1);
    sm_PSD = smoothdata(PSD,'movmedian',n/Fs);
    %=== Find the peak corresponding to wingbeat
    [~,lo_f] = min(abs(f-6));
    [~,hi_f] = min(abs(f-10));
    [~,loc] = findpeaks(log(sm_PSD(lo_f:hi_f)),n/Fs,'MinPeakProminence',0.1);
    wb_F(i) = 6+max(loc);
    %=== Plot the spectrum and wingbeat harmonics
    nexttile;
    plot(f,PSD,'.k','LineWidth',2);     xlim([1 50]);    set(gca, 'YScale', 'log');
    hold on;
    plot(f,sm_PSD,'LineWidth',4,'Color',bat_clr(i,:));
    for j = 1:5,plot(wb_F(i)*[j j],ylim,'--k');end
    hold off;   xlabel('Hz');   yticks([]);
end
sgtitle('Accelerometer Spectrum');
fig_count = saveFig(analysis_directory,batdate,fig_count,options.savefigures);

%% Assign features associated with each bat and each pair

for i=1:n_tags   
    %=== Number of flights,p-value and percentage of clustered points, flight probability, number of rewards
    bat(i).f_num = f_num(i);
    bat(i).p_val_clus = p_val_clus(i,:);
    bat(i).CH_idx = CH_idx(i,:);
    bat(i).perc_clus = perc_clus(i,:);
    bat(i).p_fly1 = p_fly1(:,i);
    bat(i).wb_F = wb_F(:,i);
    bat(i).r_num = nnz(FLIGHTS.class(FLIGHTS.id==i)== 'f')+nnz(FLIGHTS.class(FLIGHTS.id==i)== 'b');
    
    %=== Exploration ratio see Obenhaus et al., "Functional network topography of the medial entorhinal cortex." PNAS 2022)
    occupancy = histcounts2(r(~bflying(:,i)& b_corr(:,i),1,i),r(~bflying(:,i)& b_corr(:,i),2,i),cell2mat(edges_d(1,1)),cell2mat(edges_d(1,2)));
    perimeter = [occupancy(1,:),occupancy(end,:),occupancy(:,1)',occupancy(:,end)'];
    bat(i).exp_ratio_all = nnz(occupancy>Fs*5)/numel(occupancy);
    bat(i).exp_ratio_per = nnz(perimeter>Fs*5)/numel(perimeter);
    bat(i).spatial_map = occupancy./T;

    %=== Fraction of time spent close to the feeder
    for j = 1:size(r_fd,1),condition(:,j) = vecnorm(r(:,:,i)-r_fd(j,:),2,2)<d_th;end
    bat(i).feedfraction = nnz(any(condition,2))/T;
    
    %=== Fraction of feeder,social,both and other flights
    bat(i).f_flights = (nnz(FLIGHTS.class(FLIGHTS.id==i)== 'f')+nnz(FLIGHTS.class(FLIGHTS.id==i)== 'b'))/f_num(i);    % fraction of feeder flights, including social + feeder
    bat(i).s_flights = (nnz(FLIGHTS.class(FLIGHTS.id==i)== 's')+nnz(FLIGHTS.class(FLIGHTS.id==i)== 'b'))/f_num(i);    % fraction of social flights, including social + feeder    
    bat(i).b_flights = nnz(FLIGHTS.class(FLIGHTS.id==i)== 'b')/f_num(i);    
    bat(i).o_flights = nnz(FLIGHTS.class(FLIGHTS.id==i)== 'o')/f_num(i);
end

for i=1:n_pairs
    %=== Proximity/Landing index and p value
    pair(i).PI = PI(1,i);
    pair(i).PI_pval = PI_p(i);
    pair(i).PI_c = PI_c(i);
    pair(i).CI = CI(1,i);
    pair(i).CI_pval = CI_p(i);
    pair(i).CI_c = CI_c(i);
    pair(i).LI = LI(1,i,:);
    pair(i).LI_pval = LI_p(i,:);
    
    %=== Probability of simultaneous flight (observed vs calculated)
    pair(i).p_fly2 = p_fly2(i,:);
    
    %=== Maximal and median difference in the number of flights
    pair(i).max_flight_diff = max(abs(fligh_diff(:,i)));
    pair(i).med_flight_diff = median(fligh_diff(:,i));
    pair(i).sig_flight_diff = nnz(fligh_diff(:,i)>0)/T;
end

%% Save data

if options.save_data
    save([analysis_directory,'/Analyzed_Behavior_', batdate, '.mat'],...
            'bat','bat_dist','centroid','centroid_all','d_th','f_bd','f_clus','f_nn','f_smp','FLIGHTS','Group_name','n_ds','n_min','n_rep',...
            'options','pair','perc_time','pref_location','pref_time_loc','r_clus_id','r_fd','r_emb','state','samplesIn','T','true_rest');           
end

%% Plot additional figures

%=== Show Density Histograms heat-map
figure('units','normalized','outerposition',[0 0.25 1 0.35]);
for i=1:n_tags
    sbpt = subplot(1,n_tags,i);
    hist3(r(:,1:2,i),'edges',edges_d,'CdataMode','auto','edgecolor','none','FaceColor','interp');
    xlabel('x');
    xlim(r_lim(1,:)); ylim(r_lim(2,:));   title(bat_nms(i,:),['Feeds: ' ,num2str(nnz(FLIGHTS.class(FLIGHTS.id==i)== 'f')+nnz(FLIGHTS.class(FLIGHTS.id==i)== 'b'))]);
    view(90,90);  colormap(sbpt,custom_map(:,:,i)); % Change color scheme
    axis square;
end
fig_count = saveFig(analysis_directory,batdate,fig_count,options.savefigures);

%=== FIGURE: Scatter plot all bats
figure();       set(gcf, 'units','normalized','outerposition',[0 0.25 1 0.5]);
for i=1:n_tags
    subplot(131);  plot3(r(:,1,i),r(:,2,i),r(:,3,i),'k');  xlim(r_lim(1,:)); ylim(r_lim(2,:)); zlim(r_lim(3,:));  title('3D view');                 hold on;  axis equal;
    subplot(132);  plot3(r(:,1,i),r(:,2,i),r(:,3,i),'k');  xlim(r_lim(1,:)); ylim(r_lim(2,:)); zlim(r_lim(3,:));  title('Top view');   view(0,90);  hold on;  axis equal;
    subplot(133);  plot3(r(:,1,i),r(:,2,i),r(:,3,i),'k');  xlim(r_lim(1,:)); ylim(r_lim(2,:)); zlim(r_lim(3,:));  title('Door view');  view(-30,0); hold on;  axis equal;
end
hold off;
fig_count = saveFig(analysis_directory,batdate,fig_count,options.savefigures);

%=== Show flights and positions during stationary epochs (cut the first 10s and last 3s)
figure('units','normalized','outerposition',[0.3 0.2 0.35 0.35]);
tiledlayout(2,n_tags,'TileSpacing','tight');
for i = 1:n_tags
    nexttile(i);
    plot3(r(:,1,i),r(:,2,i),r(:,3,i),'-','Color', bat_clr(i,:));  xlim(r_lim(1,:)); ylim(r_lim(2,:)); zlim(r_lim(3,:));  view(90,90);   axis equal
    nexttile(n_tags+i);
    scatter3(r_qt(10*Fs:end-3*Fs,1,i),r_qt(10*Fs:end-3*Fs,2,i),r_qt(10*Fs:end-3*Fs,3,i),5,'filled','MarkerFaceColor', bat_clr(i,:));
    xlim(r_lim(1,:)); ylim(r_lim(2,:)); zlim(r_lim(3,:));
    set(gca,'XTickLabels',[],'YTickLabels',[],'ZTickLabels',[]); 
    xlabel('x','fontsize',14);ylabel('y','fontsize',14);zlabel('z','fontsize',14);
end
fig_count = saveFig(analysis_directory,batdate,fig_count,options.savefigures);

%=== Show all the clusters
figure;     
str = string(1:size(centroid_all,1));
textscatter3(centroid_all(:,1),centroid_all(:,2),centroid_all(:,3),str);    hold on;
scatter3(r_fd(:,1),r_fd(:,2),r_fd(:,3),50,'filled','MarkerFaceColor', 'r'); hold off;
fig_count = saveFig(analysis_directory,batdate,fig_count,options.savefigures);

%=== Plot all centroids, with dimensions proportional to occupancy
figure('units','normalized','outerposition',[0.3 0.3 0.2 0.3]);
for i = 1:n_tags
    scatter3(centroid{i,1}(:,1),centroid{i,1}(:,2),centroid{i,1}(:,3),samplesIn{i,1}./T*200,'filled');
    hold on;
end
axis equal; xlim(r_lim(1,:)); ylim(r_lim(2,:)); zlim(r_lim(3,:));   grid on;    hold off;
fig_count = saveFig(analysis_directory,batdate,fig_count,options.savefigures);

%=== Plot Density of states
figure('units','normalized','outerposition',[0 0.25 1 0.45]);
tiledlayout(1,10,'TileSpacing','tight');
nexttile([1,9]);
bar(state_density); hold on;    plot(xlim,size(state_templates,1)/size(state_rest,1)*[1 1]);    hold off;
xlabel('Possible states of the system');    ylabel('Fraction of total time spent into state');
title(['Fraction of explored states of all possible states: ', num2str(size(unique(state_rest,'rows'),1)/(max(id_all)-1)^n_tags,2) newline...
    'Fraction of significant states of all explored states:' num2str(nnz(state_density>size(state_templates,1)/size(state_rest,1))/size(state_templates,1),2) newline...
    'Entropy (only oberved states) = ' num2str(E_obs,2) ,'  Entropy (max) = ' num2str(E_max,2)]);
nexttile; imagesc(state_matrix);  colormap(flip(gray)); xlabel('Bat Id');   ylabel(' Positional Cluster');  title('Average Observed States');   
co = colorbar;   co.Label.String = 'Fraction of rest time';
fig_count = saveFig(analysis_directory,batdate,fig_count,options.savefigures);

%=== Plot Probabilities of flying together
figure('units','normalized','outerposition',[0.3 0.3 0.2 0.33]);
scatter(p_fly2(:,2),p_fly2(:,1),'k','filled');
xlabel('P coupled fly (independent)');  ylabel('P coupled fly (real)');
text(p_fly2(:,2),p_fly2(:,1),bat_pair_nms,'VerticalAlignment','bottom','HorizontalAlignment','right')
axis square;    xlim(ylim);     hold on;    h = refline(1,0);   hold off;   h.Color = 'k';  h.LineStyle = '--';
fig_count = saveFig(analysis_directory,batdate,fig_count,options.savefigures);

%=== FIG: Proximity index histograms
figure('units','normalized','outerposition',[0.3 0.1 0.05 0.9]);
tiledlayout(n_pairs,1,'TileSpacing','tight');
for i = 1:n_pairs
    ax(i) = nexttile;       histogram(PI(:,i),'edgecolor','none','FaceColor','k');
    yticks([]);  xlim([0 max(PI(:))]);
    y1=get(gca,'ylim');  hold on; plot([PI(1,i) PI(1,i)],y1); hold off;
    title([bat_nms(bat_pairs(i,1),:) '-' bat_nms(bat_pairs(i,2),:) ': ' num2str(PI_p(i),3)]);
end
xlabel('Proximity Index');
fig_count = saveFig(analysis_directory,batdate,fig_count,options.savefigures);

%=== FIG: Coupling index histograms
figure('units','normalized','outerposition',[0.3 0.1 0.05 0.9]);
tiledlayout(n_pairs,1,'TileSpacing','tight');
for i = 1:n_pairs
    ax(i) = nexttile;       histogram(CI(:,i),'edgecolor','none','FaceColor','k');
    yticks([]);  xlim([0 max(CI(:))]);
    y1=get(gca,'ylim');  hold on; plot([CI(1,i) CI(1,i)],y1); hold off;
    title([bat_nms(bat_pairs(i,1),:) '-' bat_nms(bat_pairs(i,2),:) ': ' num2str(CI_p(i),3)]);
end
xlabel('Coupling Index');
fig_count = saveFig(analysis_directory,batdate,fig_count,options.savefigures);

%=== FIG: Landing index histograms
figure();   set(gcf, 'units','normalized','outerposition',[0.3 0.1 0.05 0.9]);
tiledlayout(n_pairs,2,'TileSpacing','tight');
for i = 1:n_pairs
    ax(i) = nexttile;       histogram(LI(:,i,1),'edgecolor','none','FaceColor','k');
    yticks([]);  xlim([0 max(LI(:))]);
    y1=get(gca,'ylim');  hold on; plot([LI(1,i,1) LI(1,i,1)],y1); hold off;
    title([bat_nms(bat_pairs(i,1),:) '-' bat_nms(bat_pairs(i,2),:) ': ' num2str(LI_p(i,1),3)]);
    ax(i+1) = nexttile;     histogram(LI(:,i,2),'edgecolor','none','FaceColor','k');
    yticks([]);  xlim([0 max(LI(:))]);
    y1=get(gca,'ylim');  hold on; plot([LI(1,i,2) LI(1,i,2)],y1); hold off;
    title([bat_nms(bat_pairs(i,2),:) '-' bat_nms(bat_pairs(i,1),:) ': ' num2str(LI_p(i,2),3)]);
end
xlabel('Proximity Index');
fig_count = saveFig(analysis_directory,batdate,fig_count,options.savefigures);

%=== Network graph with PI
A = zeros(n_tags);
for i = 1:n_pairs
    if PI_p(i)<1e-3,     A(bat_pairs(i,1),bat_pairs(i,2)) = 0.7;
    elseif PI_p(i)<1e-2, A(bat_pairs(i,1),bat_pairs(i,2)) = 0.5;
    elseif PI_p(i)<5e-2, A(bat_pairs(i,1),bat_pairs(i,2)) = 0.3;
    else,                A(bat_pairs(i,1),bat_pairs(i,2)) = 0.01;    
    end
end
G = graph(A,cellstr(bat_nms),'upper');
figure();
G.Edges.LWidths = 20*G.Edges.Weight;    p_NG = plot(G); title('PI Network'); 
set(p_NG,'LineWidth',G.Edges.LWidths,'MarkerSize',10,'NodeLabelColor',bat_clr,'NodeColor',bat_clr,'NodeFontSize',15,'NodeFontWeight','bold','EdgeColor',0.5*[1 1 1]);
fig_count = saveFig(analysis_directory,batdate,fig_count,options.savefigures);

%=== Network graph with CI
A = zeros(n_tags);
for i = 1:n_pairs
    if CI_p(i)<1e-3,     A(bat_pairs(i,1),bat_pairs(i,2)) = 0.7;
    elseif CI_p(i)<1e-2, A(bat_pairs(i,1),bat_pairs(i,2)) = 0.5;
    elseif CI_p(i)<5e-2, A(bat_pairs(i,1),bat_pairs(i,2)) = 0.3;
    else,                A(bat_pairs(i,1),bat_pairs(i,2)) = 0.01;    
    end
end
G = graph(A,cellstr(bat_nms),'upper');
figure();
G.Edges.LWidths = 20*G.Edges.Weight;    p_NG = plot(G); title('CI Network'); 
set(p_NG,'LineWidth',G.Edges.LWidths,'MarkerSize',10,'NodeLabelColor',bat_clr,'NodeColor',bat_clr,'NodeFontSize',15,'NodeFontWeight','bold','EdgeColor',0.5*[1 1 1]);
fig_count = saveFig(analysis_directory,batdate,fig_count,options.savefigures);

%=== Network graph with LI
A = zeros(n_tags);
for i = 1:n_pairs
    if LI_p(i,1)<1e-3,     A(bat_pairs(i,1),bat_pairs(i,2)) = 0.7;
    elseif LI_p(i,1)<1e-2, A(bat_pairs(i,1),bat_pairs(i,2)) = 0.5;
    elseif LI_p(i,1)<5e-2, A(bat_pairs(i,1),bat_pairs(i,2)) = 0.3;
    else,                  A(bat_pairs(i,1),bat_pairs(i,2)) = 0.01;    
    end
    if LI_p(i,2)<1e-3,     A(bat_pairs(i,2),bat_pairs(i,1)) = 0.7;
    elseif LI_p(i,2)<1e-2, A(bat_pairs(i,2),bat_pairs(i,1)) = 0.5;
    elseif LI_p(i,2)<5e-2, A(bat_pairs(i,2),bat_pairs(i,1)) = 0.3;
    else,                  A(bat_pairs(i,2),bat_pairs(i,1)) = 0.01;    
    end
end
G = digraph(A,cellstr(bat_nms));
figure();
G.Edges.LWidths = 20*G.Edges.Weight;    p_NG = plot(G); title('LI Network'); 
set(p_NG,'LineWidth',G.Edges.LWidths,'MarkerSize',10,'NodeLabelColor',bat_clr,'NodeColor',bat_clr,'ArrowSize',10,'NodeFontSize',15,'NodeFontWeight','bold','EdgeColor',0.5*[1 1 1]);
fig_count = saveFig(analysis_directory,batdate,fig_count,options.savefigures);

%=== Coordinated flight
figure();   set(gcf, 'units','normalized','outerposition',[0.3 0.3 0.2 0.3]);
for i = 1:n_tags
    plot(t(f_smp{i,1})/60,1:f_num(i),'Color', bat_clr(i,:),'LineWidth',3); hold on;
end
hold off;   xlabel('Time into session (min)');  ylabel('Flight number');    axis square;
fig_count = saveFig(analysis_directory,batdate,fig_count,options.savefigures);
