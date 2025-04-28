function flight_clus = FlightClus_AF_v0(r_vector,f_vector,varargin)
%Flight_clus_AF extracts flight clusters from trajectory data, by using
%hierarchical agglomerative clustering (or k-means) on 3D points along the trajectory. 

%INPUTS:
% r_vector = time (rows) x 3D position (columns)
% f_vector = time (rows) x flight (column): 0s (not flying) or 1s (flying)
% varargin = 

%% Parameters and overrides
x2=2.8; x1=-2.8; y2=2.6;  y1=-2.6;  z1=0; z2=2.30;  %Flight volume coordinates
Fs = 100;                                           %sampling frequency (Hz)
ds_clus = 10;                                       %number of 3D-points/flight for clustering (depending of the clustering algorithm, usually >6 points is enough)
pca_features = 0;                                   %if using PCA
Frechet_clus = 0;                                   %if using Frechet-distance for clustering
k_means = 0;                                        %if using k-means
alpha = 3.5;                                        %main clustering parameter
reassign = 1;                                       %re-order clusters
N_min = 3;                                          %min number of flights for being a cluster
Frechet = 0;                                        %if performing Frechet-based calculations
n_plot = 5;                                         %number of clusters/figure
debug = 1;

%User input overrides default params
if nargin > 1
    nparams=length(varargin);
    for i=1:2:nparams
        switch (varargin{i})
            case 'Fs'
                Fs=varargin{i+1};
            case 'Points'
                ds_clus = varargin{i+1};
            case 'N_min'
                N_min = varargin{i+1};
            case 'Alpha'
                alpha = varargin{i+1};
            case 'Frechet'
                Frechet_clus = varargin{i+1};
            case 'PCA'
                pca_features = varargin{i+1};
            case 'Debug'
                debug = varargin{i+1};
        end
    end
end

%Sanity check: a flight is not cut at the beginning or end of the session
if f_vector(1) == 1 || f_vector(end) == 1
    warning('INCOMPLETE FLIGHTS AT START/STOP OF THE SESSION');
end

%Extract flight number, start and stop and calculate velocity
num_flights = nnz(diff(f_vector)>0);
f_start = find(diff(f_vector)>0)+1;
f_stop = find(diff(f_vector)<0);
v = diff(r_vector,1,1).*Fs; v=[zeros(1,3); v];   v_abs = vecnorm(v,2,2);


%% Clustering flights

%Cut out flights, downsample to ds_clus positions per flight
all_flights = NaN(3,max(f_stop-f_start)+1,num_flights);          %3D matrix with all flights (position)
all_flights_vel = NaN(1,max(f_stop-f_start)+1,num_flights);      %3D matrix with all flights (velocity)
all_flights_ds = NaN(3,ds_clus,num_flights);                     %3D matrix with all flights(downsampled position)

for nf = 1 : num_flights
    all_flights(:,1:(f_stop(nf)-f_start(nf)+1),nf) = r_vector(f_start(nf):f_stop(nf),:)';
    all_flights_vel(1,1:(f_stop(nf)-f_start(nf)+1),nf) = v_abs(f_start(nf):f_stop(nf),:)';
    all_flights_ds(:,:,nf) = interp1(linspace(1,100,f_stop(nf)-f_start(nf)+1)',r_vector(f_start(nf):f_stop(nf),:),linspace(1,100,ds_clus)','spline')';
%             %Uncomment if you want to see how the downsampled flights look like
%             plot3(all_flights(1,:,nf),all_flights(2,:,nf),all_flights(3,:,nf),'Color','b');   hold on;
%             plot3(all_flights_ds(1,:,nf),all_flights_ds(2,:,nf),all_flights_ds(3,:,nf),'.','Color','r','MarkerSize',10);  hold off;
%             w = waitforbuttonpress;
end

%Define X matrix of features for clustering (downsampled coordinates, stacked together)
if ~Frechet_clus
    X = reshape(permute(all_flights_ds,[2 1 3]),[],num_flights)';   %so then X = #flights x #features
else
    X = all_flights_ds; 
end

%If dimensionality reduction is needed
if pca_features
    [coeff,score,latent] = pca(X);     X = score(:,1:5);
end

%k-means or hierarchical clustering (with euclidean distance and shortest linkage)
if k_means
    n_clusters = 15;    idx = kmeans(X,n_clusters);
else
    figure();   set(gcf, 'units','normalized','outerposition',[0.5 0.2 0.15 0.7]);
    dist = alpha;                               %linkage distance
    if ~Frechet_clus
        Y = pdist(X,'euclidean');
    else
        pairs = nchoosek(1:num_flights,2);
        for n = 1:length(pairs)
            Y(n) = DiscreteFrechetDist(squeeze(X(:,:,pairs(n,1))),squeeze(X(:,:,pairs(n,2))));
        end
    end
    Z = linkage(Y,'centroid');
    subplot(311); edges_f = 10.^linspace(log10(min(Y)),log10(max(Y)),100);  histogram(Y,edges_f);    set(gca,'XScale','log');    set(gca,'YScale','log');  
    subplot(312); histogram(Y);
    subplot(313); hLines = dendrogram(Z,0);  hold on;    refline(0,dist);    hold off;
    
    idx = cluster(Z,'Cutoff',dist,'Criterion','inconsistent');
    title([num2str(length(unique(idx))) ' clusters']);
    %ylim([0 10]);
end

if debug
    flight_clus = length(unique(idx));
    return
end

%Create structure with flight start stop frames, id of the trajectory
clear flight;
flight.strt_frame = ceil(f_start)';
flight.stop_frame = ceil(f_stop)';
flight.pos = all_flights;
flight.vel = all_flights_vel;
flight.id = idx;
flight.Fs = Fs;
flight.ifd = diff(flight.strt_frame)./flight.Fs;

%Sort structure according to cluster id
clear flight_sorted;
[flight_sorted.id,I] = sort(flight.id);
flight_sorted.strt_frame = flight.strt_frame(I);
flight_sorted.stop_frame = flight.stop_frame(I);
flight_sorted.pos = flight.pos(:,:,I);
flight_sorted.vel = flight.vel(:,:,I);
flight_sorted.Fs = flight.Fs;
flight_sorted.N = size(flight_sorted.id,1);

%Assign isolated clusters to cluster #flights+1
[Ns,b] = histc(flight_sorted.id,unique(flight_sorted.id));
flight_sorted.id(Ns(b)<N_min) = size(all_flights,3)+1;
id_surv_clusters = unique(flight_sorted.id);
n_surv_clusters = size(id_surv_clusters,1);

%Create final structure flight.clus after re-assignment
%unclustered flights are in the last cluster
clear flight_clus;
flight_clus.id = flight_sorted.id;
flight_clus.strt_frame = flight_sorted.strt_frame;
flight_clus.stop_frame = flight_sorted.stop_frame;
flight_clus.pos = flight_sorted.pos;
flight_clus.vel = flight_sorted.vel;
flight_clus.Fs = flight_sorted.Fs;
flight_clus.N = flight_sorted.N;
for jj=1:n_surv_clusters;
    flight_clus.id(flight_sorted.id == id_surv_clusters(jj)) = jj;
end
id_surv_clusters = unique(flight_clus.id);

%Re-assign id for convenience, with cluster 1 of non-clustered flights and
%then clusters reordered in descending order of crowdedness
if reassign
    new_ord = [];   [~,new_ord] = sort(histc(flight_clus.id,id_surv_clusters(1:end-1)),'descend');
    new_ord = [new_ord; id_surv_clusters(end)];
    new_ord = circshift(new_ord,1);
    reassign_matrix =(flight_clus.id == new_ord');
    for jj=1:n_surv_clusters
        flight_clus.id(reassign_matrix(:,jj)) = jj;
    end
end

%Calculate trajectory length, duration in s and interflight (take-off to take-off)
for ii = 1:flight_sorted.N
    trajectory = unique(round(flight_clus.pos(:,~isnan(flight_clus.pos(1,:,ii)),ii)',3),'rows','stable');
    flight_clus.length(ii)= arclength(trajectory(:,1),trajectory(:,2),trajectory(:,3),'s');
    flight_clus.dur(ii) = (flight_clus.stop_frame(ii)-flight_clus.strt_frame(ii))./flight_clus.Fs;
end
flight_clus.ifd = flight.ifd;

%% Visualize
col = hsv(n_surv_clusters); 
figure(); set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
tiledlayout(max([5, round(n_surv_clusters/4)]),12,'TileSpacing','tight');
for jj=1:n_surv_clusters
    id = find(flight_clus.id==jj);
    
    %subplot(3,n_surv_clusters,jj);
    nexttile;  
    plot(r_vector(:,1),r_vector(:,2),':','Color',[0.8 0.8 0.8],'MarkerSize',0.001);
    xlim([x1 x2]); ylim([y1 y2]); 
    xlabel('x');    ylabel('y');    
    hold on;
    avg_take_off = [];
    for ii=1:size(id,1)
        title(['Cluster' num2str(jj) ' (' num2str(size(id,1)) ' flights),' num2str(size(id,1)/num_flights*100,2) '%'])
        plot(flight_clus.pos(1,:,id(ii)),flight_clus.pos(2,:,id(ii)),'-','LineWidth',1.5,'Color', col(jj,:));
        avg_take_off = [avg_take_off flight_clus.pos(:,1,id(ii))];
       
    end
    take_off = mean(avg_take_off,2);
    textscatter(take_off(1),take_off(2),"Take-off");
    hold off;
    
    %subplot(3,n_surv_clusters,n_surv_clusters+jj);
    nexttile;
    plot(r_vector(:,1),r_vector(:,3),':','Color',[0.8 0.8 0.8],'MarkerSize',0.001);
    xlim([x1 x2]); ylim([z1 z2]); 
    xlabel('x');    ylabel('y');    
    hold on;
    avg_take_off = [];
    for ii=1:size(id,1)
        plot(flight_clus.pos(1,:,id(ii)),flight_clus.pos(3,:,id(ii)),'-','LineWidth',1.5,'Color', col(jj,:));
        avg_take_off = [avg_take_off flight_clus.pos(:,1,id(ii))];
       
    end
    take_off = mean(avg_take_off,2);
    textscatter(take_off(1),take_off(3),"Take-off");
    hold off;
    
    %subplot(3,n_surv_clusters,n_surv_clusters*2+jj);
    nexttile;
    histogram(flight_clus.dur(id));
    xlim([0 15]);   xlabel('Duration(s)');  ylabel('Counts');
end

%Plot flights stats
figure();
subplot(1,3,1);     histogram(flight_clus.length);   xlabel('Flight length (m)');        ylabel('Counts');
subplot(1,3,2);     histogram(flight_clus.dur);      xlabel('Fligth duration (s)');      ylabel('Counts');
subplot(1,3,3);     histogram(flight_clus.ifd);      xlabel('Inter-Fligth (s)');         ylabel('Counts');

%Calculate entropy of the session
%[Pr,b] = histc(flight.id,unique(flight.id));    Pr = Pr./flight_clus.N;     Entropy = -sum(Pr.*log2(Pr));

%% Calculate the Frechet distance between flights
if Frechet
    %Calculate distance between all the flights
    flight_pairs = nchoosek(1:flight_clus.N,2);
    F_dist = zeros(size(flight_pairs,1),1);
    for i = 1:size(flight_pairs,1)
        P = flight_clus.pos(:,:,flight_pairs(i,1));     P = P(:,~isnan(P(1,:)))';   P = downsample(P,3);
        Q = flight_clus.pos(:,:,flight_pairs(i,2));     Q = Q(:,~isnan(Q(1,:)))';   Q = downsample(Q,3);
        [F_dist(i),~] = DiscreteFrechetDist(P,Q);
    end
    
    %Intra-cluster Frechet distance
    [Nf,~] = histc(flight_clus.id,unique(flight_clus.id));
    F_dist_intra = nan(max(Nf)*(max(Nf)-1)/2,n_surv_clusters);
    if n_surv_clusters>1
        for jj = 1:n_surv_clusters
            disp(jj);
            flight_pairs = nchoosek(find(flight_clus.id==jj),2);
            for i = 1:size(flight_pairs,1)
                P = flight_clus.pos(:,:,flight_pairs(i,1));     P = P(:,~isnan(P(1,:)))';   P = downsample(P,3);
                Q = flight_clus.pos(:,:,flight_pairs(i,2));     Q = Q(:,~isnan(Q(1,:)))';   Q = downsample(Q,3);
                [F_dist_intra(i,jj),~] = DiscreteFrechetDist(P,Q);
            end
        end
    end
    
    %Extra-cluster Frechet distance
    if n_surv_clusters>1
        flight_pairs = nchoosek(1:flight_clus.N,2);
        idx_del = false(size(flight_pairs,1),1);
        for jj = 2:n_surv_clusters
            to_delete = nchoosek(find(flight_clus.id==jj),2);
            idx_del = idx_del | ismember(flight_pairs,to_delete,'rows') | ismember(flight_pairs,to_delete(:,[2 1]),'rows');
        end
        F_dist_extra = F_dist(~idx_del);
    else
        F_dist_extra = F_dist;
    end
    
    %Stats for different clusters (cv duration, length, Frechet and interflight interval for same cluster)
    Clus = [];
    for jj = 1:n_surv_clusters
        Clus.cv_dur(jj) = std(flight_clus.dur(ismember(flight_clus.id,jj)))/mean(flight_clus.dur(ismember(flight_clus.id,jj)));
        Clus.cv_len(jj) = std(flight_clus.length(ismember(flight_clus.id,jj)))/mean(flight_clus.length(ismember(flight_clus.id,jj)));
        Clus.cv_fre(jj) = mean(F_dist_intra(:,jj),'omitnan')/mean(flight_clus.length(ismember(flight_clus.id,jj)));
        Clus.if_int(jj) = median(diff(flight_clus.strt_frame(ismember(flight_clus.id,jj)))./flight_clus.Fs);
    end
    
    figure();   set(gcf, 'units','normalized','outerposition',[0.1 0.25 0.8 0.5]);
    ax_1 = subplot(1,3,1);
    hfl = histogram(F_dist,'facealpha',.5,'edgecolor','none');
    xlabel('Frechet distance (m)'); ylabel('Counts'); title('All flights');
    
    ax_2 = subplot(1,3,2);
    for jj = 1:n_surv_clusters
        histogram(F_dist_intra(:,jj),hfl.BinEdges,'facecolor',col(jj,:),'facealpha',.5,'edgecolor','none'); hold on;
    end
    hold off;   xlabel('Frechet distance (m)'); title('Intra clusters 1,2,...');
    
    ax_3 = subplot(1,3,3);
    histogram(F_dist_extra,hfl.BinEdges,'facecolor','b','facealpha',.5,'edgecolor','none');   hold on;
    histogram(F_dist_intra(:,2:end),hfl.BinEdges,'facecolor','k','facealpha',.5,'edgecolor','none');   hold off;
    xlabel('Frechet distance (m)');            title('Intra VS Extra');
    
    linkaxes([ax_1,ax_2,ax_3],'y');
    %set([ax_1,ax_2,ax_3], 'YScale', 'log');
    
end
end