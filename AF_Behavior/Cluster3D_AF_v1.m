function [id,centroid,samplesIn,p_val,CH_idx] = Cluster3D_AF_v1(X,n_min,distance,n_ds,varargin)

%% Performs agglomerative hierarchical clustering on positional data
%% Input
% X:                matrix of observations x features (3D coordinates)
% n_min:            minimum number of samples to consitute a cluster
% distance:         linkage distance
% n_ds:             downsampling factor
% figure flag:      ...,'Fig_Flag',0) to avoid plotting figures

%% Output
% id: id of the cluster for each sample (rows), 0 is for unclustered
% centroid: cluster centroid coordinates for each cluster; (NaN,NaN,NaN) for the unclustered
% samplesIn: number of samples spent in that cluster
% p_val: p value of a KS test between the observed distribution of...
% ...distances and the one corresponding to uniform spreading around the room
% CH_idx: Calinski-Harabasz index for 2-10 clusters

%% Parameters and variable arguments

r_lim = [-2.9 2.9; -2.6 2.6; 0 2.30];   % Room boundaries
edges_dist = 10.^linspace(-3,1,100);    % Edges for distance histograms (from 1 mm to 10 m)
fig_flag = 0;
use_dbscan = 0;
if nargin > 1
    nparams=length(varargin);
    for i=1:2:nparams
        switch (varargin{i})
            case 'Fig_Flag'
                fig_flag=varargin{i+1};
        end
    end
end

%% Clustering

%=== Get downsampled X (X_ds) and calculate pairwise distances between points
X_ds = downsample(X,n_ds);                                      % downsample positions                                      
Y = pdist(X_ds,'euclidean');                                    % calculate euclidean distance between points

%=== Cluster points using dbscan or hierarchical clustering
if use_dbscan
    [id_temp,corepts] = dbscan(Y,distance,round(n_min/n_ds),'Distance','precomputed');  % Cluster via DBSCAN algorithm
    id_temp(~corepts) = 0;                                                              % Assign id = 0 to noise points
else
    Z = linkage(Y,'single');                                        % calculate linkage distance
    id_temp = cluster(Z,'Cutoff',distance,'Criterion','distance');  % perform clustering
    [Ns,~,b] = histcounts(id_temp,[1:max(id_temp)+1]-0.5);          % calculate number of points/cluster
    id_temp(Ns(b)<round(n_min/n_ds)) = 0;                           % clusters with less than n_min points are merged into the 0 cluster
end

%=== Reorder clusters
[~,~,id_reord] = unique(id_temp);                               % get unique cluster ids 
if nnz(id_temp)> 0,  id_reord = id_reord - 1;    end            % if cluster 0 is present, shift all cluster ids: 0,1,2...
id_ds = id_reord;  clus_ids = unique(id_ds);                    % assign re-ordered cluster ids

%=== Assign cluster id to the unlabeled points in X
if n_ds>1
    id = NaN(size(X,1),1);                                          % initialize the id vector for all the entries in X
    id(downsample([1:size(id,1)],n_ds)',:)=id_ds;                   % assign the cluster id to the downsampled points
    %id = fillmissing(id,'nearest');                                 % fill the missing ids with the nearest value
    %to_correct = find(diff(id))+1+round([-n_ds/2+1:n_ds/2-2]);      % when there is a change of id, look at [-n_ds/2:n_ds/2] points ...
    %to_correct = to_correct(to_correct>0);                          % ...and assign their ids based on the closest downsampled id
    %idx = knnsearch(X_ds,X(to_correct,:));
    %id(to_correct) = id_ds(idx);
    idx = knnsearch(X_ds,X);
    id = id_ds(idx);
else
    id = id_ds;
end

%=== Calculate cluster centroids and occupancy
centroid = NaN(length(clus_ids),size(X,2));
samplesIn = NaN(length(clus_ids),1);
for i =1:length(clus_ids)
    if clus_ids(i)~=0
    centroid(i,:) = mean(X(id == clus_ids(i),:));
    samplesIn(i,:) = nnz(id == clus_ids(i));
    end
end

%=== Calculate Calinski-Harabasz index for 2-10 clusters
CH_idx = [];
warning('off','all');
for i = 1:2
    eva = evalclusters(X_ds,'kmeans','CalinskiHarabasz','KList',[2:20]);
    CH_idx = [CH_idx; eva.CriterionValues];
end
warning('on','all');
CH_idx = mean(CH_idx);

%% Test the hyphotesis that points are uniformly spread across the available faces

%=== Simulated distribution for points spread across the walls (avoiding floor)
r_rnd = []; for j=1:3 ;  r_rnd(:,j) = r_lim(j,1)+[r_lim(j,2)-r_lim(j,1)]*rand(5000,1);   end
%=== Project on the 6 faces
r_rnd(1    :1e3,1) = r_lim(1,1);    
r_rnd(1e3+1:2e3,1) = r_lim(1,2);
r_rnd(2e3+1:3e3,2) = r_lim(2,1);    
r_rnd(3e3+1:4e3,2) = r_lim(2,2);
r_rnd(4e3+1:5e3,3) = r_lim(3,2);    
%scatter3(r_rnd(:,1),r_rnd(:,2),r_rnd(:,3))
Y_rnd = pdist(r_rnd);

[~,p_val] = kstest2(datasample(Y,min(1e4,length(Y)),'Replace',false),datasample(Y_rnd,min(1e4,length(Y)),'Replace',false));

%% Visualize

if fig_flag
    clus_col = hsv(length(clus_ids));
    %=== Distribution of distances 
    nexttile;
    histogram(Y,edges_dist,'edgecolor','none','FaceAlpha',0.8,'FaceColor',[0,0,0],'Normalization','probability');    set(gca,'xScale','log');    set(gca,'YScale','log');
    yticks([]); line([0.2,0.2],ylim);   xlabel('Within-bat Distance (m)');      grid on;
    %=== Distribution of distances (cdf)
    nexttile;
    histogram(Y,'edgecolor','none','FaceAlpha',0.8,'FaceColor',[0,0,0],'Normalization','cdf');  ylim([0,1]);
    hold on;    histogram(Y_rnd,'edgecolor','none','FaceAlpha',0.2,'FaceColor',[1,0,0],'Normalization','cdf');  hold off;
    yticks([]); line([0.2,0.2],ylim);   xlabel('Within-bat Distance (m)');      grid on;
    title(['p value =  ', num2str(p_val)])
    %=== CH-indexe
    nexttile;
    plot(eva.InspectedK,CH_idx);    xlabel('N Clusters');   ylabel('Calinski-Harabasz index');
    %=== Cluster Visualization
    nexttile;
    for i =1:length(clus_ids)
        if clus_ids(i)~=0
            scatter3(X(id == clus_ids(i),1),X(id == clus_ids(i),2),X(id == clus_ids(i),3),5,'filled','MarkerFaceColor',clus_col(i,:)); hold on;
        else
            scatter3(X(id == clus_ids(i),1),X(id == clus_ids(i),2),X(id == clus_ids(i),3),1,'.k'); hold on;  
        end
        xlim([-2.9 2.9]); ylim([-2.6 2.6]); zlim([0 2.3]);
    end
    hold off;
    title(['Clustered points: ', num2str(nnz(id)/length(id)*100,3) ,' %']);
end

end

