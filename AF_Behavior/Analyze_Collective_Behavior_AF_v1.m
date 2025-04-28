%% Analyze_Collective_Behavior_AF_v1
% Script for the analysis of collective behavior (March 2022)
% Requires Processed data from Ciholas recordings and Group name (D or F)

%=== Load data
extracted_BHVfile = dir(fullfile(cd, 'Extracted_Behavior_*'));          load(extracted_BHVfile.name);
disp(['Processing ', extracted_BHVfile.name, '...']);

%=== Parameters
n_ds = 100;                                                                                                     % downsampling factor for clustering positions
n_min = 10*Fs;                                                                                                  % minimum number of samples to be considered a cluster (at least 10s)
n_rep = 1000;                                                                                                   % number of repetitions for the d_th shuffling
n_pairs = length(bat_pairs);                                                                                    % number of bat pairs
edges_dist = 10.^linspace(-3,1,100);                                                                            % edges for distance histograms (from 1 mm to 10 m)
edges_d = {r_lim(1,1):(r_lim(1,2)-r_lim(1,1))/20:r_lim(1,2) r_lim(2,1):(r_lim(2,2)-r_lim(2,1))/20:r_lim(2,2)};  % edges for density histogram
bat_ids = [1:n_tags]';                                                                                          % bat identities
d_th = 0.27;                                                                                                    % threshold distance for various classifications  
for i = 1:n_tags; for j = 1:3; custom_map(:,j,i) = linspace(1,bat_clr(i,j))'; end; end                          % custom graded colormap
if Group_name == 'D', r_fd = [-0.5,0.05,0.45];                                                                  % Feeder position  for D bats
else, r_fd = [2.77,0.82,1.75; 2.79,-0.99,1.64; 2.78,1.29,0.84; 2.78,-1.43,0.80]; end                            % Feeder positions for F bats
options.save_data = 1;
options.clusterFl = 1;                                                                                          % Performing flight clustering or not

%=== Create analysis folder for storing the results
if options.save_data
    analysis_directory=fullfile(pwd,['Analysis_',datestr(now, 'yymmdd_HHMM')]);
    if ~exist(analysis_directory,'dir')
        mkdir(analysis_directory);
    end
end

%=== Correction vector for excluding tails with v>v_th but no formal 'flight'
b_corr = ones(size(bflying));
for i=1:n_tags
    if f_num(i)>0
        v_diff = diff(v_abs(:,i)>v_th);
        transitions = find(v_diff);
        if v_diff(transitions(1))   ==-1,b_corr(1:transitions(1),i) = 0;    end
        if v_diff(transitions(end)) == 1,b_corr(transitions(end):end,i) = 0;end
    end
end

%% Interbat distances

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

%=== Show Density Histograms heat-map
figure();   set(gcf, 'units','normalized','outerposition',[0 0.25 1 0.35]);
for i=1:n_tags
    sbpt = subplot(1,n_tags,i);
    hist3(r(:,1:2,i),'edges',edges_d,'CdataMode','auto','edgecolor','none','FaceColor','interp');
    xlabel('x');
    xlim(r_lim(1,:)); ylim(r_lim(2,:));   title(bat_nms(i,:));  view(90,90);  colormap(sbpt,custom_map(:,:,i)); % Change color scheme
    axis square;
end

%=== Show flights and positions during stationary epochs (cut the first 10s and last 3s)
figure();   set(gcf, 'units','normalized','outerposition',[0.3 0.2 0.35 0.3]);
tiledlayout(2,n_tags,'TileSpacing','tight');
for i = 1:n_tags
    nexttile(i);
    plot3(r(:,1,i),r(:,2,i),r(:,3,i),'-','Color', bat_clr(i,:));  xlim(r_lim(1,:)); ylim(r_lim(2,:)); zlim(r_lim(3,:));  view(90,90);
    nexttile(n_tags+i);
    scatter3(r_qt(10*Fs:end-3*Fs,1,i),r_qt(10*Fs:end-3*Fs,2,i),r_qt(10*Fs:end-3*Fs,3,i),30,'filled','MarkerFaceColor', bat_clr(i,:));
    xlim(r_lim(1,:)); ylim(r_lim(2,:)); zlim(r_lim(3,:));
end

%=== Pool positions from all the bats and cluster them
X = []; for i = 1:n_tags; X = [X;r(~bflying(:,i)& b_corr(:,i),:,i)]; end
figure();   set(gcf, 'units','normalized','outerposition',[0.3 0.2 0.2 0.2]);
[~,centroid_all] = Cluster3D_AF_v1(X,n_min,0.2,100);

%=== Correct feeder position(s) and show all the clusters
if Group_name == 'D'
    [~,fd_idx] = min(vecnorm(r_fd-centroid_all,2,2));
    r_fd(1,:) = centroid_all(fd_idx,:);                                                  
else
    for i=1:4
    [feeder_distance,fd_idx] = min(vecnorm(r_fd(i,:)-centroid_all,2,2));
    if feeder_distance < 0.2        % Do not correct if further than 20 cm
        r_fd(i,:) =  centroid_all(fd_idx,:);
    end
    end
end
figure;     
str = string(1:size(centroid_all,1));
textscatter3(centroid_all(:,1),centroid_all(:,2),centroid_all(:,3),str);    hold on;
scatter3(r_fd(:,1),r_fd(:,2),r_fd(:,3),50,'filled','MarkerFaceColor', 'r'); hold off;

%=== Single bat clusters
r_clus_id = NaN(T,n_tags);  centroid = cell(n_tags,1);  samplesIn = cell(n_tags,1); p_val_clus = zeros(n_tags,1); perc_clus = zeros(n_tags,1);
figure();   set(gcf, 'units','normalized','outerposition',[0.3 0 0.27 1]);
tiledlayout(n_tags,3,'TileSpacing','tight');
for i = 1:n_tags
    X = r(~bflying(:,i)& b_corr(:,i),:,i);                                                  % Consider only stationary epochs
    [id,centroid{i,1},samplesIn{i,1},p_val_clus(i,:)] = Cluster3D_AF_v1(X,n_min,0.2,n_ds);  % Perform clustering
    perc_clus(i,:) = nnz(id)/length(id);                                                    % fraction of points belonging to a cluster
    r_clus_id(~bflying(:,i)& b_corr(:,i),i)=id;                                             % Assign cluster id across entire session (nan when bat is flying)
end

%=== Plot all centroids, with dimensions proportional to occupancy
figure();   set(gcf, 'units','normalized','outerposition',[0.3 0.3 0.2 0.3]);
for i = 1:n_tags
    scatter3(centroid{i,1}(:,1),centroid{i,1}(:,2),centroid{i,1}(:,3),samplesIn{i,1}./T*200,'filled');
    hold on;
end
axis equal; xlim(r_lim(1,:)); ylim(r_lim(2,:)); zlim(r_lim(3,:));   grid on;
hold off;

%% Cluster flights (f_clus = 0 if options.clusterFl =0)

alpha_clus = 1.2;         %Parameter for flight clustering
f_cls = cell(n_tags,1);     
for i = 1:n_tags
    if options.clusterFl 
        f_clus(i) = FlightClus_AF_v2(squeeze(r(:,:,i)),bflying(:,i),'Alpha',alpha_clus,'Frechet',1,'Points',10);
        f_cls{i,1} = f_clus(i).sorted_ids;
    else   
        f_clus = struct([]);
        f_cls{i,1} = 0.*f_smp{i,1};  
    end
end

%% Create FLIGHTS table with several features

warning('off','all');
c = 1;  FLIGHTS = table(); 
for i = 1:n_tags
    for j=1:f_num(i)
        FLIGHTS.id(c,:) = i;                                      % id of the bat flying
        FLIGHTS.smp1(c,:) = f_smp{i,1}(j,1);                      % sample takeoff
        FLIGHTS.smp2(c,:) = f_smp{i,2}(j,1);                      % sample landing
        FLIGHTS.t1(c,:) = t(f_smp{i,1}(j,1));                     % time takeoff
        FLIGHTS.t2(c,:) = t(f_smp{i,2}(j,1));                     % time landing
        FLIGHTS.pclus1(c,:) = r_clus_id(f_smp{i,1}(j,1)-1,i);     % positional cluster at take off
        FLIGHTS.pclus2(c,:) = r_clus_id(f_smp{i,2}(j,1)+1,i);     % positional cluster at landing
        FLIGHTS.fclus(c,:) = f_cls{i,1}(j,1);                    % flight cluster
        cond_f = any(vecnorm(r(f_smp{i,2}(j,1),:,i)-r_fd,2,2)<d_th);
        cond_s = any(vecnorm(r(f_smp{i,2}(j,1),:,i)-r(f_smp{i,2}(j,1),:,bat_ids(bat_ids~=i)),2,2)<d_th);
        if      cond_f && ~cond_s
            FLIGHTS.class(c,:) = 'f';                             % exclusive feeding flight
        elseif ~cond_f && cond_s
            FLIGHTS.class(c,:) = 's';                             % exclusive social flight
        elseif  cond_f && cond_s
            FLIGHTS.class(c,:) = 'b';                             % feeding + social flight
        else
            FLIGHTS.class(c,:) = 'o';                             % other flight
        end
        c = c+1;
    end
end
warning('on','all');

%% Relationships between bats: Proximity Index

%=== Proximity index calculation
PI = zeros(n_rep+1,n_pairs);   
PI_p = zeros(n_pairs,1);
for i = 1:n_pairs
    PI(1,i) = nnz(bat_dist(:,i)<d_th)/T;                                    % PI: fraction of frames < d_th
    frames_to_shift = randi([60*Fs T-60*Fs],1,n_rep);                       % Random number of frames between 1 min and session time-1min
    parfor n = 1:n_rep
        r_shfl = circshift(r(:,:,bat_pairs(i,2)),frames_to_shift(n),1);     % Shift the position of the second bat in the couple
        PI(n+1,i) = nnz(vecnorm(r(:,:,bat_pairs(i,1))-r_shfl,2,2)<d_th)/T;  % Recalculate the PI 
    end
    PI_p(i,:) = nnz(PI(2:end,i)>PI(1,i))/n_rep;                             % Calculate the p value as the fraction of PI > observed PI
end

%=== FIG: Proximity index histograms
figure();   set(gcf, 'units','normalized','outerposition',[0.3 0.1 0.05 0.9]);
tiledlayout(n_pairs,1,'TileSpacing','tight');
for i = 1:n_pairs
    ax(i) = nexttile;       histogram(PI(:,i),'edgecolor','none','FaceColor','k');
    yticks([]);  xlim([0 max(PI(:))]);
    y1=get(gca,'ylim');  hold on; plot([PI(1,i) PI(1,i)],y1); hold off;
    title([bat_nms(bat_pairs(i,1),:) '-' bat_nms(bat_pairs(i,2),:) ': ' num2str(PI_p(i),3)]);
end
xlabel('Proximity Index');

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
G.Edges.LWidths = 20*G.Edges.Weight;    p_NG = plot(G); 
set(p_NG,'LineWidth',G.Edges.LWidths,'MarkerSize',10,'NodeLabelColor',bat_clr,'NodeColor',bat_clr,'NodeFontSize',15,'NodeFontWeight','bold','EdgeColor',0.5*[1 1 1]);


%% Relationships between bats: Landing Index

%=== Landing index calculation
LI = zeros(n_rep+1,n_pairs,2);              % Initialize LI matrix: (1+n_rep,bat pairs in forward order,bat pairs in reverse order) 
LI_p = zeros(n_pairs,2);                    % Initialize p-value matrix
for i = 1:n_pairs
    landing_temp = f_bd{bat_pairs(i,1),2}<d_th;                             % Temporary matrix with all the landings<d_th for bat_pairs(i,1)
    LI(1,i,1) = nnz(landing_temp(:,bat_pairs(i,2)))/f_num(bat_pairs(i,1));  % LI(1,i,1) = fraction of flights for bat_pairs(i,1) landing less than d_th from bat_pairs(i,2)
    landing_temp = f_bd{bat_pairs(i,2),2}<d_th;                             % Temporary matrix with all the landings<d_th for bat_pairs(i,2)
    LI(1,i,2) = nnz(landing_temp(:,bat_pairs(i,1)))/f_num(bat_pairs(i,2));  % LI(1,i,2) = fraction of flights for bat_pairs(i,2) landing less than d_th from bat_pairs(i,1)
    frames_to_shift = randi([60*Fs T-60*Fs],1,n_rep);                       % Random number of frames between 1 min and session time-1min
    for n = 1:n_rep
        r_shfl = circshift(r(:,:,bat_pairs(i,2)),frames_to_shift(n),1);     % Shift the position of the second bat in the couple
        % Calculate LI as the fraction of flights for bat_pairs(i,1) landing less than d_th from r_shfl 
        LI(n+1,i,1) = nnz(vecnorm(r(f_smp{bat_pairs(i,1),2},:,bat_pairs(i,1))-r_shfl(f_smp{bat_pairs(i,1),2},:),2,2)<d_th)/f_num(bat_pairs(i,1));
        r_shfl = circshift(r(:,:,bat_pairs(i,1)),frames_to_shift(n),1);     % Shift the position of the first bat in the couple
        % Calculate LI as the fraction of flights for bat_pairs(i,2) landing less than d_th from r_shfl
        LI(n+1,i,2) = nnz(vecnorm(r(f_smp{bat_pairs(i,2),2},:,bat_pairs(i,2))-r_shfl(f_smp{bat_pairs(i,2),2},:),2,2)<d_th)/f_num(bat_pairs(i,2));
    end
    LI_p(i,1) = nnz(LI(2:end,i,1)>LI(1,i,1))/n_rep;                             % Calculate the p value as the fraction of LI > observed LI
    LI_p(i,2) = nnz(LI(2:end,i,2)>LI(1,i,2))/n_rep;                             % Calculate the p value as the fraction of LI > observed LI
end

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
G.Edges.LWidths = 20*G.Edges.Weight;    p_NG = plot(G); 
set(p_NG,'LineWidth',G.Edges.LWidths,'MarkerSize',10,'NodeLabelColor',bat_clr,'NodeColor',bat_clr,'ArrowSize',10,'NodeFontSize',15,'NodeFontWeight','bold','EdgeColor',0.5*[1 1 1]);

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

%=== Coordinated flight
figure();   set(gcf, 'units','normalized','outerposition',[0.3 0.3 0.2 0.3]);
for i = 1:n_tags
    plot(t(f_smp{i,1})/60,1:f_num(i),'Color', bat_clr(i,:),'LineWidth',3); hold on;
end
hold off;   xlabel('Time into session (min)');  ylabel('Flight number');    axis square;

%=== Calculate the difference in the number of flights across the session
figure();   set(gcf, 'units','normalized','outerposition',[0.3 0.1 0.2 0.9]);
tiledlayout(2,1,'TileSpacing','tight');
cum_flights = [zeros(1,n_tags); cumsum(max(diff(bflying),0))];
kernel = zeros(Fs*600,1);  kernel(1) = -1; kernel(end) = 1;
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
nexttile;
for i = 1:n_pairs
    local_diff = conv(fligh_diff(:,i), kernel, 'same');
    plot([1:length(local_diff)]./Fs,local_diff-50*i,'LineWidth',2);    hold on;
    text(length(local_diff)/Fs,-50*i,[bat_nms(bat_pairs(i,1),:) '-' bat_nms(bat_pairs(i,2),:)]);
    rectangle('Position',[0 -50*i-10 length(local_diff)/Fs 20],'FaceColor',[0 0 0 0.1],'EdgeColor','none');
    rectangle('Position',[0 -50*i-5  length(local_diff)/Fs 10],'FaceColor',[0 0 0 0.2],'EdgeColor','none');
end
hold off;
xlabel('Time (s)'); ylabel('Difference in the number of flights');

%% Group and quantify different features for each bat and for each couple

for i=1:n_tags   
    
    %=== Number of flights,p-value and percentage of clustered points, flight probability
    bat(i).f_num = f_num(i);
    bat(i).p_val_clus = p_val_clus(i,:);
    bat(i).perc_clus = perc_clus(i,:);
    bat(i).p_fly1 = p_fly1(:,i);
    
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
    pair(i).LI = LI(1,i,:);
    pair(i).LI_pval = LI_p(i,:);
    
    %=== Probability of simultaneous flight (observed vs calculated)
    pair(i).p_fly2 = p_fly2(i,:);
    
    %=== Maximal and median difference in the number of flights
    pair(i).max_flight_diff = max(abs(fligh_diff(:,i)));
    pair(i).med_flight_diff = median(fligh_diff(:,i));
    pair(i).sig_flight_diff = nnz(fligh_diff(:,i)>0)/T;

end

%% Save figures and data

if options.save_data
    figHandles = findall(0,'Type','figure');
    for i = 1:numel(figHandles)
        saveas(figHandles(i),[analysis_directory, '/', batdate '_figure' num2str(numel(figHandles)+1-i) '.png']);
    end
    close all;
    save([analysis_directory,'/Analyzed_Behavior_', batdate, '.mat'],...
            'bat','centroid','centroid_all','d_th','f_bd','f_clus','f_nn','f_smp','FLIGHTS','Group_name',...
            'n_ds','n_min','n_rep','pair','perc_time','pref_location','pref_time_loc','r_fd','samplesIn','T');           
end

%% Relationships between bats: Heading difference when flying

% heading_diff = nan(T,length(bat_pairs),n_rep+1);
% for i = 1:length(bat_pairs)
%     frames_to_shift = randi([60*Fs T-60*Fs],1,n_rep);
%     cpled = find(~isnan(angle(:,bat_pairs(i,1))).*~isnan(angle(:,bat_pairs(i,2))));
%     heading_diff(cpled,i,1) = rad2deg(angdiff(angle(cpled,bat_pairs(i,1)),angle(cpled,bat_pairs(i,2))));
%     for n = 1:n_rep
%         angle_shfl = circshift(angle  (:,bat_pairs(i,2)),frames_to_shift(n),1);     % Shift the angle of the second bat in the couple
%         cpled = find(~isnan(angle(:,bat_pairs(i,1))).*~isnan(angle_shfl));
%         heading_diff(cpled,i,1+n) = rad2deg(angdiff(angle(cpled,bat_pairs(i,1)),angle(cpled,bat_pairs(i,2))));
%     end
% end
% 
% figure();   set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
% for i = 1:length(bat_pairs)
%     subplot(2,5,i);
%     if nnz(~isnan(heading_diff(:,i)))>0
%         alpha = deg2rad(heading_diff(:,i,1)); alpha = alpha(~isnan(alpha));
%         alpha_shfl = squeeze(heading_diff(:,i,2:end));  alpha_shfl = deg2rad(alpha_shfl(:)); alpha_shfl = alpha_shfl(~isnan(alpha_shfl));
%         polarhistogram(alpha,20,'edgecolor','none','Normalization','probability');    hold on;
%         polarhistogram(alpha_shfl,20,'edgecolor','none','Normalization','probability');    hold off;
%         title([bat_nms(bat_pairs(i,1),:) '-' bat_nms(bat_pairs(i,2),:) ' (' num2str(0.01*nnz(~isnan(heading_diff(:,i)))) 's)']);
%     end
% end
