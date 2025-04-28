function id = Cluster3D_AF_v0(X,n_min,distance,varargin)

fig_flag = 1;
if nargin > 1
    nparams=length(varargin);
    for i=1:2:nparams
        switch (varargin{i})
            case 'Fig_Flag'
                fig_flag=varargin{i+1};
        end
    end
end

%===Calculate euclidean distances between points and linkage distance
Y = pdist(X,'euclidean');
Z = linkage(Y,'single');

%===Perform agglomerative hierarchical clustering
id_temp = cluster(Z,'Cutoff',distance,'Criterion','distance');
edges = [1:max(id_temp)+1]-0.5;
[Ns,~,b] = histcounts(id_temp,edges);

%===Clusters with less than n_min points goes on the 'unclustered'
id_temp(Ns(b)<n_min) = max(id_temp)+1;

%===Get unique clusters and plot
[clus_ids,~,id_reord] = unique(id_temp);
clus_col = jet(length(clus_ids));
if fig_flag
    for i =1:length(clus_ids)
        scatter3(X(id_temp == clus_ids(i),1),X(id_temp == clus_ids(i),2),X(id_temp == clus_ids(i),3),'filled','MarkerFaceColor',clus_col(i,:)); hold on;
    end
    title(['d =' num2str(distance) ' m, ' num2str(length(clus_ids)-1) ' clusters']);
    hold off;
end

%===Reassign cluster ids, with nan for the unclustered
id_reord(id_reord == max(id_reord)) = nan;
id = id_reord;

end

