%% Script for looking at the collective behavior across days

%=== Load data and aggregate them in the Multi_Day structure
clr;    Folder = cd;    FileList = dir(fullfile(Folder, '**', 'Analyzed_Behavior_*'));
if ~isempty(dir(fullfile(Folder,'Dataset_name.mat'))),load('Dataset_name.mat');disp(['DATSET: === ', Dataset_name, ' ===']);end
Multi_Day = struct([]);
for iii = 1:length(FileList)
    cd(FileList(iii).folder);
    load(FileList(iii).name);
    path = strsplit(FileList(iii).folder,'\');
    Multi_Day(iii).T = T;                         % Total Number of acquired samples
    Multi_Day(iii).Session = path{1,5};           % Session Date
    Multi_Day(iii).r_fd = r_fd;                   % Coordinates of the feeder for that day
    Multi_Day(iii).centroid = centroid;           % Centroid coordinates for the preferred locations
    Multi_Day(iii).samplesIn = samplesIn;         % Number of samples spent in the preferred locations
    Multi_Day(iii).bat = bat;                     % Structure containing several single-bat features
    Multi_Day(iii).pair = pair;                   % Structure containing several paired-bat features
    Multi_Day(iii).pref_location = pref_location; % First 3 preferred locations
    Multi_Day(iii).f_nn = f_nn;                   % Nearest-bat distances at landing and takeoff
    Multi_Day(iii).f_bd = f_bd;                   % Bat distances at landing and takeoff
    Multi_Day(iii).d_th = d_th;                   % Threshold distance for later classifications
    Multi_Day(iii).FLIGHTS = FLIGHTS;             % FLIGHTS table
    Multi_Day(iii).perc_time = perc_time;         % Percentage of time with 0,1,2,... bats flying
    Multi_Day(iii).lt_d = vecnorm(FLIGHTS.r1(2:end,:)-FLIGHTS.r2(1:end-1,:),2,2);     % Get the takeoff to previous landing distance
    Multi_Day(iii).lt_d = Multi_Day(iii).lt_d(~diff(FLIGHTS.id),:);                     % Remove spurious distances across different bat ids
    disp([num2str(length(FileList)-iii),' remaining sessions to load...']);
end
cd(Folder); clearvars -except Multi_Day Group_name
T = struct2table(Multi_Day);                    % Rearrange data in a table

%=== Parameters and additional useful variables
n_sessions = size(Multi_Day,2);                                                                                 % Number of sessions
n_tags = size(Multi_Day(1).centroid,1);                                                                         % Number of bats
bat_ids = [1:n_tags]';                                                                                          % Bat identities
r_lim = [-2.9 2.9; -2.6 2.6; 0 2.30];                                                                           % Room boundaries
diagonal = vecnorm(diff(r_lim,1,2));                                                                            % Room diagonal
edges_d = {r_lim(1,1):(r_lim(1,2)-r_lim(1,1))/20:r_lim(1,2) r_lim(2,1):(r_lim(2,2)-r_lim(2,1))/20:r_lim(2,2)};  % Edges for density histogram
Fs = 100;                                                                                                       % Sampling frequency (Hz) for common time
if Group_name == 'D';bat_nms = ['Dai'; 'Den'; 'Dia'; 'Dor'; 'Dum'; 'Im1'; 'Im2'];                               % Bat names for D-bats
else; bat_nms = ['Far'; 'Fer'; 'Fia'; 'For'; 'Ful'; 'Im1'; 'Im2']; end                                          % Bat names for F-bats
bat_nms = bat_nms(1:n_tags,:);                                                                                  % Bat Names
bat_pairs = nchoosek(1:n_tags,2);                                                                               % Bat Pairs
n_pairs = length(bat_pairs);                                                                                    % Number of bat pairs
bat_pair_nms = [bat_nms(bat_pairs(:,1),:), '-'.*ones(length(bat_pairs),1), bat_nms(bat_pairs(:,2),:)];          % Bat Pairs Names
bat_clr = lines(n_tags);                                                                                        % Bat Colors
for i = 1:n_tags; for j = 1:3; custom_map(:,j,i) = linspace(1,bat_clr(i,j))'; end; end                          % Custom graded colormap(level,RGB,bat)
d_th = unique(T.d_th);                                                                                          % Threshold distance for later classifications

%% Look at the position of the bats and the feeder(s) across sessions

centroids_all = cell2mat(vertcat(T.centroid{:,1}));
samples_all = cell2mat(vertcat(T.samplesIn{:,1}));
view_angle = [-37.5,30;0,90;0,0;90,0];
figure; set(gcf, 'units','normalized','outerposition',[0.1 0.3 0.2 0.43]);
tiledlayout(2,2,'TileSpacing','tight');
for i=1:4
    nexttile;
    scatter3(centroids_all(:,1),centroids_all(:,2),centroids_all(:,3),samples_all./median(T.T)*10,'filled','MarkerFaceColor', 'k');
    view(view_angle(i,:));   
    xlim(r_lim(1,:)); ylim(r_lim(2,:)); zlim(r_lim(3,:));
    axis equal; set(gca,'XTickLabels',[],'YTickLabels',[],'ZTickLabels',[]); 
    xlabel('x','fontsize',16);ylabel('y','fontsize',16);zlabel('z','fontsize',16); 
end

figure;
for j=1:n_sessions
    feeder_pos = Multi_Day(j).r_fd;
    scatter3(feeder_pos(:,1),feeder_pos(:,2),feeder_pos(:,3),50,'filled','MarkerFaceColor', 'r');  hold on;
    textscatter3(feeder_pos(:,1),feeder_pos(:,2),feeder_pos(:,3),string(j*ones(1,size(feeder_pos,1))));
end
hold off;
xlim(r_lim(1,:)); ylim(r_lim(2,:)); zlim(r_lim(3,:));

%% Look at the probability of 0,1,2,... bats flying across sessions

figure('units','normalized','outerposition',[0 0.3 1 0.2]);   tiledlayout(1,n_sessions,'TileSpacing','none');
for j=1:n_sessions
    nexttile;   bar([0:n_tags]',cell2mat(T.perc_time(j,:)),'k');    ylim([0 1]);    xlabel('# Bats flying');    title(['Session ', num2str(j)]);
end

figure;
bar([0:n_tags]',mean(cell2mat(T.perc_time'),2),'k');    ylim([0 1]);    xlabel('# Bats flying');    ylabel('Probability');

%% Are the bats resting across all the available space or do they choose specific locations?

%=== Number of flights/hour
data = reshape([T.bat(:,:).f_num],[],n_tags)'./(T.T/Fs/3600)';
figure('units','normalized','outerposition',[0.1 0.3 0.1 0.3]);
PlotDistr_AF_v0(data,bat_clr,'Number of flights');

%=== Exploration ratio all
data = reshape([T.bat(:,:).exp_ratio_all],[],n_tags)';
figure('units','normalized','outerposition',[0.2 0.3 0.1 0.3]);
PlotDistr_AF_v0(data,bat_clr,'Exploration Ratio (2D)','Y_lim',[0 0.5]);

%=== Exploration ratio perimeter
data = reshape([T.bat(:,:).exp_ratio_per],[],n_tags)';
figure('units','normalized','outerposition',[0.3 0.3 0.1 0.3]);
PlotDistr_AF_v0(data,bat_clr,'Exploration Ratio (Perimeter)','Y_lim',[0 0.5]);

%=== P-values against uniform distribution across the walls
data = reshape([T.bat(:,:).p_val_clus],[],n_tags)';
figure('units','normalized','outerposition',[0.4 0.3 0.1 0.3]);
PlotDistr_AF_v0(data,bat_clr,'P-val against uniform distribution','Y_lim',[0 0.001]);

%=== Fraction clustered samples
data = reshape([T.bat(:,:).perc_clus],[],n_tags)';
figure('units','normalized','outerposition',[0.5 0.3 0.1 0.3]);
PlotDistr_AF_v0(data,bat_clr,'Fraction Clustered Samples','Y_lim',[0.95 1]);

%=== Calinski-Harabasz index for 2-10 clusters
figure('units','normalized','outerposition',[0.2 0.3 0.5 0.3]);
tiledlayout(1,n_tags,'TileSpacing','tight');
data = reshape([T.bat(:,:).CH_idx],[],n_sessions,n_tags);
for i=1:n_tags
    ydata = normalize(data(:,:,i),1,'range');   xdata = [1:size(ydata,1)];
    nexttile;   plot(ydata,'Color',[0 0 0 0.3]);  hold on;
    plot(mean(ydata,2),'Color',bat_clr(i,:));
    ciplot(mean(ydata,2)-std(ydata,[],2)./sqrt(size(ydata,2)),mean(ydata,2)+std(ydata,[],2)./sqrt(size(ydata,2)),xdata,bat_clr(i,:));   alpha(0.3);  hold off;
    ylabel('Calinski-Harabasz index (norm.)');
    xlabel('Number of clusters');
end

%% Do they spend a lot of time on the feeder(s)? Do they fly a lot to and from the feeder?

%=== Calculate amount of samples spent close to feeder vs away from it
t_feed = []; t_else = [];
for j = 1:n_sessions
    if Group_name == 'D'
        for i = 1:n_tags
            t_feed = [t_feed; Multi_Day(j).samplesIn{i,1}(find(vecnorm(Multi_Day(j).centroid{i,:}-Multi_Day(j).r_fd,2,2)<d_th),1)/Fs];
            t_else = [t_else; Multi_Day(j).samplesIn{i,1}(find(vecnorm(Multi_Day(j).centroid{i,:}-Multi_Day(j).r_fd,2,2)>d_th),1)/Fs];
        end
    else
        for f=1:4,for i = 1:n_tags
                t_else = [t_else; Multi_Day(j).samplesIn{i,1}(find(vecnorm(Multi_Day(j).centroid{i,:}-Multi_Day(j).r_fd(f,:),2,2)>d_th),1)/Fs];
                t_feed = [t_feed; Multi_Day(j).samplesIn{i,1}(find(vecnorm(Multi_Day(j).centroid{i,:}-Multi_Day(j).r_fd(f,:),2,2)<d_th),1)/Fs];
            end,      end
    end
end

%=== Test the same distribution hypothesis for t_feed and t_else
[~,p_val] = kstest2(t_else,t_feed);

%=== Plot the cdf for time close to the feeder vs time away from it and
figure; set(gcf, 'units','normalized','outerposition',[0.2 0.7 0.1 0.3]);
bins = 10.^linspace(log10(min([t_else;t_feed])),log10(max([t_else;t_feed])),30);
histogram(t_feed,bins,'edgecolor','none','FaceAlpha',0.5,'FaceColor',[1,1,0],'Normalization','cdf');  hold on;
histogram(t_else,bins,'edgecolor','none','FaceAlpha',0.5,'FaceColor',[0,0,0],'Normalization','cdf');  hold off;
set(gca,'xScale','log');    xlim([bins(1), bins(end)]); ylim([0 1]); axis square;
xlabel('Time spent into location (s)');  ylabel('CDF');  legend('Feeder','Others','Location','northwest');
title(['p = ', num2str(p_val,3)]);

%=== Percentage of time close to the feeder
data = reshape([T.bat(:,:).feedfraction],[],n_tags)';
figure; set(gcf, 'units','normalized','outerposition',[0.1 0.3 0.1 0.3]);
PlotDistr_AF_v0(data,bat_clr,'Fraction Time close to feeder');

%=== Fraction of flights to the feeder
data = reshape([T.bat(:,:).f_flights],[],n_tags)';
figure; set(gcf, 'units','normalized','outerposition',[0.2 0.3 0.1 0.3]);
PlotDistr_AF_v0(data,bat_clr,'Fraction Flights to Feeder','Y_lim',[0 1]);

%=== How does the number of flights to feeder evolve across sessions
f_prc = reshape([T.bat(:,:).f_flights],[],n_tags);
f_num = reshape([T.bat(:,:).f_num],[],n_tags);
figure; set(gcf, 'units','normalized','outerposition',[0.2 0.5 0.1 0.4]);
tiledlayout(2,1,'TileSpacing','none');
nexttile;   plot(f_prc,'LineWidth',3);        ylabel('% Fligths to Feeder');  xticks([]);
nexttile;   plot(f_prc.*f_num,'LineWidth',3); ylabel('# Fligths to Feeder');  xlabel('Session');

%=== Test if the fraction of flights to feeder is = 0.5
for i = 1:n_tags
    if ~kstest(zscore([T.bat(:,i).f_flights]'))
        [~,p] = ttest([T.bat(:,i).f_flights]'-0.5);
        test = 'Paired t test';
    else
        [p,~] = signrank([T.bat(:,i).f_flights]'-0.5);
        test = 'Wilcoxon signed rank test';
    end
    disp(['For bat n.', num2str(i),', flights to feeder are ',...
        num2str(mean([T.bat(:,i).f_flights])/mean(1-[T.bat(:,i).f_flights]'),3),...
        ' times flights to non-feeder locations, p =', num2str(p,3), ' ', test]);
end

%% If only a minor part of their behavior is devoted to feeding, what are they doing?

%=== Rearrange all the landing distances to the nn-bat
nnb_landings = reshape([T.f_nn{:,:}],n_tags,2,n_sessions);  nnb_landings = squeeze(nnb_landings(:,2,:))';

%=== Define the edges and centers for binning distances
edges_dist = 10.^[-2:0.1:log10(diagonal)];
edge_dist_cent = edges_dist(1:end-1)+ diff(edges_dist)/2;

%=== Concatenate all nn-bat distances at landing and count their frequency
nnb_distances = vertcat(nnb_landings{:,:});
[p_distance,~] = histcounts(nnb_distances,edges_dist,'Normalization','probability');

%=== Show the distribution of nn-bat distances
figure; set(gcf, 'units','normalized','outerposition',[0.1 0.3 0.2 0.3]);
tiledlayout(1,2,'TileSpacing','none');
nexttile;   histogram(nnb_distances,edges_dist,'Normalization','probability');  set(gca, 'XScale', 'log');  xlabel('NN-bat distance at landing (m)');   yticks([]);
nexttile;   h = radialHist_AF_v0(p_distance,edge_dist_cent,[0,1],[0 360],[0 0 0]); h.CDataMode = 'auto'; axis equal; xlabel('NN-bat distance at landing (m)');    yticks([]);

%=== Fraction of flights to the feeder
data = reshape([T.bat(:,:).f_flights],[],n_tags)';
figure; set(gcf, 'units','normalized','outerposition',[0.3 0.3 0.1 0.3]);
PlotDistr_AF_v0(data,bat_clr,'Fraction Flights to Feeder','Y_lim',[0 1]);

%=== Fraction of flights to another bat
data = reshape([T.bat(:,:).s_flights],[],n_tags)';
figure; set(gcf, 'units','normalized','outerposition',[0.4 0.3 0.1 0.3]);
PlotDistr_AF_v0(data,bat_clr,'Fraction Flights to Bat','Y_lim',[0 1]);

%=== Fraction of flights to another bat and feeder
data = reshape([T.bat(:,:).b_flights],[],n_tags)';
figure; set(gcf, 'units','normalized','outerposition',[0.5 0.3 0.1 0.3]);
PlotDistr_AF_v0(data,bat_clr,'Fraction Flights to Feeder&Bat','Y_lim',[0 1]);

%=== Fraction of flights to other locations
data = reshape([T.bat(:,:).o_flights],[],n_tags)';
figure; set(gcf, 'units','normalized','outerposition',[0.6 0.3 0.1 0.3]);
PlotDistr_AF_v0(data,bat_clr,'Fraction Flights to Wall','Y_lim',[0 1]);

%% A lot of landings happen close to another bat: is there any bat preference or bat-specific landing distance?

%=== Collect and rearrange all the landing events
all_landings = reshape([T.f_bd{:,:}],n_tags,2,n_sessions);  all_landings = squeeze(all_landings(:,2,:))';

%=== Define the edges and centers for binning distances
edges_dist = 10.^[-2:0.1:log10(diagonal)];
edge_dist_cent = edges_dist(1:end-1)+ diff(edges_dist)/2;

%=== Plot the distribution of distances in a radial histogram (custom plot)
angular_span = 360/(n_tags-1);
for i = 1:n_tags
    all_distances = vertcat(all_landings{:,i});
    oth_bats = bat_ids(bat_ids~=i);
    figure;   set(gcf, 'units','normalized','outerposition',[0.1+(i-1)*0.11 0.6 0.11 0.27]);
    for s = 1:n_tags-1
        [p_distance,~] = histcounts(all_distances(:,oth_bats(s)),edges_dist,'Normalization','probability');
        radialHist_AF_v0(p_distance,edge_dist_cent,[0,1],angular_span*[s-1,s],bat_clr(oth_bats(s),:));              hold on;
        plot3([0 d_th*cos(deg2rad(mean(angular_span*[s-1,s])))],[0 d_th*sin(deg2rad(mean(angular_span*[s-1,s])))],[1 1],'k');
        xticks([]); yticks([]); axis off ;
    end
    hold off;
end

%=== Quantify the percentages of landings on specific bats
landing_table = zeros(n_tags,n_tags,n_sessions);
for i = 1:n_tags
    for j =1:n_sessions
        all_distances = vertcat(all_landings{j,i});
        all_distances = all_distances<d_th;
        oth_bats = bat_ids(bat_ids~=i);
        for s=1:n_tags-1
            landing_table(i,oth_bats(s),j) = nnz(all_distances(:,oth_bats(s)))./size(all_distances,1);
        end
    end
end

% %=== Look at actual histograms, as a sanity check
% figure();   set(gcf, 'units','normalized','outerposition',[0.3 0.3 0.16 0.38]);
% tiledlayout(n_tags,n_tags,'TileSpacing','none');
% edges_dist = 10.^[-3:0.1:0];
% for i = 1:n_tags
%     all_distances = vertcat(all_landings{:,i});
%     oth_bats = bat_ids(bat_ids~=i);
%     for k=1:n_tags
%         nexttile;   histogram(all_distances(:,k),edges_dist,'Normalization','probability','edgecolor','none');    xlim([0 1]);
%     end
% end

%% Bats are often seen close to each other, is this a consequence of shared place preferences? Is there more than that? Proximity,Landing and Coupling indexes

%=== Collect proximity and landing indexes, and the associated p-values
data_Pi = reshape([T.pair(:,:).PI],[],n_pairs);
data_Pv = reshape([T.pair(:,:).PI_pval],[],n_pairs);
data_Li = reshape([T.pair(:,:).LI],[],n_pairs,2);
data_Lv = reshape([T.pair(:,:).LI_pval],2,[],n_pairs);
data_Lv = permute(data_Lv,[2 3 1]);
data_Ci = reshape([T.pair(:,:).CI],[],n_pairs);
data_Cv = reshape([T.pair(:,:).CI_pval],[],n_pairs);

%=== Look at the network graphs across sessions
figure;   set(gcf, 'units','normalized','outerposition',[0.1 0.2 0.8 0.7]);
tiledlayout(6,n_sessions,'TileSpacing','none');

%=== Proximity index (preference)
for j=1:n_sessions
    A = zeros(n_tags);
    for i = 1:n_pairs
        if data_Pv(j,i)<1e-3,     A(bat_pairs(i,1),bat_pairs(i,2)) = 0.7;
        elseif data_Pv(j,i)<1e-2, A(bat_pairs(i,1),bat_pairs(i,2)) = 0.5;
        elseif data_Pv(j,i)<5e-2, A(bat_pairs(i,1),bat_pairs(i,2)) = 0.3;
        else,                A(bat_pairs(i,1),bat_pairs(i,2)) = 0.01;
        end
    end
    G = graph(A,cellstr(bat_nms),'upper');  G.Edges.LWidths = 20*G.Edges.Weight;
    nexttile;   p_NG = plot(G); title(['Session ', num2str(j),': ' ,T.Session{j,1}]);
    set(p_NG,'LineWidth',G.Edges.LWidths,'MarkerSize',10,'NodeLabelColor',bat_clr,'NodeColor',bat_clr,'NodeFontSize',15,'NodeFontWeight','bold','EdgeColor',0.5*[1 1 1]);
end
%=== Landing index (preference)
for j=1:n_sessions
    A = zeros(n_tags);
    for i = 1:n_pairs
        if data_Lv(j,i,1)<1e-3,     A(bat_pairs(i,1),bat_pairs(i,2)) = 0.7;
        elseif data_Lv(j,i,1)<1e-2, A(bat_pairs(i,1),bat_pairs(i,2)) = 0.5;
        elseif data_Lv(j,i,1)<5e-2, A(bat_pairs(i,1),bat_pairs(i,2)) = 0.3;
        else,                       A(bat_pairs(i,1),bat_pairs(i,2)) = 0.01;
        end
        if data_Lv(j,i,2)<1e-3,     A(bat_pairs(i,2),bat_pairs(i,1)) = 0.7;
        elseif data_Lv(j,i,2)<1e-2, A(bat_pairs(i,2),bat_pairs(i,1)) = 0.5;
        elseif data_Lv(j,i,2)<5e-2, A(bat_pairs(i,2),bat_pairs(i,1)) = 0.3;
        else,                  A(bat_pairs(i,2),bat_pairs(i,1)) = 0.01;
        end
    end
    G = digraph(A,cellstr(bat_nms));  G.Edges.LWidths = 20*G.Edges.Weight;
    nexttile;   p_NG = plot(G);
    set(p_NG,'LineWidth',G.Edges.LWidths,'MarkerSize',10,'NodeLabelColor',bat_clr,'NodeColor',bat_clr,'ArrowSize',10,'NodeFontSize',15,'NodeFontWeight','bold','EdgeColor',0.5*[1 1 0]);
end
%=== Coupling index (preference)
for j=1:n_sessions
    A = zeros(n_tags);
    for i = 1:n_pairs
        if data_Cv(j,i)<1e-3,     A(bat_pairs(i,1),bat_pairs(i,2)) = 0.7;
        elseif data_Cv(j,i)<1e-2, A(bat_pairs(i,1),bat_pairs(i,2)) = 0.5;
        elseif data_Cv(j,i)<5e-2, A(bat_pairs(i,1),bat_pairs(i,2)) = 0.3;
        else,                A(bat_pairs(i,1),bat_pairs(i,2)) = 0.01;
        end
    end
    G = graph(A,cellstr(bat_nms),'upper');  G.Edges.LWidths = 20*G.Edges.Weight;
    nexttile;   p_NG = plot(G); title(['Session ', num2str(j),': ' ,T.Session{j,1}]);
    set(p_NG,'LineWidth',G.Edges.LWidths,'MarkerSize',10,'NodeLabelColor',bat_clr,'NodeColor',bat_clr,'NodeFontSize',15,'NodeFontWeight','bold','EdgeColor',0.5*[1 0.8 0.2]);
end
%=== Proximity index (avoidance)
for j=1:n_sessions
    A = zeros(n_tags);
    for i = 1:n_pairs
        if data_Pv(j,i)>1-1e-3,     A(bat_pairs(i,1),bat_pairs(i,2)) = 0.7;
        elseif data_Pv(j,i)>1-1e-2, A(bat_pairs(i,1),bat_pairs(i,2)) = 0.5;
        elseif data_Pv(j,i)>1-5e-2, A(bat_pairs(i,1),bat_pairs(i,2)) = 0.3;
        else,                A(bat_pairs(i,1),bat_pairs(i,2)) = 0.01;
        end
    end
    G = graph(A,cellstr(bat_nms),'upper');  G.Edges.LWidths = 20*G.Edges.Weight;
    nexttile;   p_NG = plot(G);
    set(p_NG,'LineWidth',G.Edges.LWidths,'MarkerSize',10,'NodeLabelColor',bat_clr,'NodeColor',bat_clr,'NodeFontSize',15,'NodeFontWeight','bold','EdgeColor',0.5*[1 0 0]);
end
%=== Landing landing (avoidance)
for j=1:n_sessions
    A = zeros(n_tags);
    for i = 1:n_pairs
        if data_Lv(j,i,1)>1-1e-3,     A(bat_pairs(i,1),bat_pairs(i,2)) = 0.7;
        elseif data_Lv(j,i,1)>1-1e-2, A(bat_pairs(i,1),bat_pairs(i,2)) = 0.5;
        elseif data_Lv(j,i,1)>1-5e-2, A(bat_pairs(i,1),bat_pairs(i,2)) = 0.3;
        else,                       A(bat_pairs(i,1),bat_pairs(i,2)) = 0.01;
        end
        if data_Lv(j,i,2)>1-1e-3,     A(bat_pairs(i,2),bat_pairs(i,1)) = 0.7;
        elseif data_Lv(j,i,2)>1-1e-2, A(bat_pairs(i,2),bat_pairs(i,1)) = 0.5;
        elseif data_Lv(j,i,2)>1-5e-2, A(bat_pairs(i,2),bat_pairs(i,1)) = 0.3;
        else,                  A(bat_pairs(i,2),bat_pairs(i,1)) = 0.01;
        end
    end
    G = digraph(A,cellstr(bat_nms));  G.Edges.LWidths = 20*G.Edges.Weight;
    nexttile;   p_NG = plot(G);
    set(p_NG,'LineWidth',G.Edges.LWidths,'MarkerSize',10,'NodeLabelColor',bat_clr,'NodeColor',bat_clr,'ArrowSize',10,'NodeFontSize',15,'NodeFontWeight','bold','EdgeColor',0.5*[0.8 0.3 0]);
end
%=== Coupling index (avoidance)
for j=1:n_sessions
    A = zeros(n_tags);
    for i = 1:n_pairs
        if data_Cv(j,i)>1-1e-3,     A(bat_pairs(i,1),bat_pairs(i,2)) = 0.7;
        elseif data_Cv(j,i)>1-1e-2, A(bat_pairs(i,1),bat_pairs(i,2)) = 0.5;
        elseif data_Cv(j,i)>1-5e-2, A(bat_pairs(i,1),bat_pairs(i,2)) = 0.3;
        else,                A(bat_pairs(i,1),bat_pairs(i,2)) = 0.01;
        end
    end
    G = graph(A,cellstr(bat_nms),'upper');  G.Edges.LWidths = 20*G.Edges.Weight;
    nexttile;   p_NG = plot(G);
    set(p_NG,'LineWidth',G.Edges.LWidths,'MarkerSize',10,'NodeLabelColor',bat_clr,'NodeColor',bat_clr,'NodeFontSize',15,'NodeFontWeight','bold','EdgeColor',0.5*[1 0 1]);
end

%=== Average network
figure;   set(gcf, 'units','normalized','outerposition',[0.1 0.2 0.1 0.5]);
tiledlayout(3,1,'TileSpacing','none');
ave_Pv = median(data_Pv,1);
A = zeros(n_tags);
for i = 1:n_pairs
    if ave_Pv(:,i)<1e-3,     A(bat_pairs(i,1),bat_pairs(i,2)) = 0.7;
    elseif ave_Pv(:,i)<1e-2, A(bat_pairs(i,1),bat_pairs(i,2)) = 0.5;
    elseif ave_Pv(:,i)<5e-2, A(bat_pairs(i,1),bat_pairs(i,2)) = 0.3;
    else,                A(bat_pairs(i,1),bat_pairs(i,2)) = 0.01;
    end
end
G = graph(A,cellstr(bat_nms),'upper');  G.Edges.LWidths = 20*G.Edges.Weight;
nexttile;   p_NG = plot(G);
set(p_NG,'LineWidth',G.Edges.LWidths,'MarkerSize',10,'NodeLabelColor',bat_clr,'NodeColor',bat_clr,'NodeFontSize',15,'NodeFontWeight','bold','EdgeColor',0.5*[1 1 1]);
ave_Lv = median(data_Lv,1);
A = zeros(n_tags);
for i = 1:n_pairs
    if data_Lv(1,i,1)<1e-3,     A(bat_pairs(i,1),bat_pairs(i,2)) = 0.7;
    elseif data_Lv(1,i,1)<1e-2, A(bat_pairs(i,1),bat_pairs(i,2)) = 0.5;
    elseif data_Lv(1,i,1)<5e-2, A(bat_pairs(i,1),bat_pairs(i,2)) = 0.3;
    else,                       A(bat_pairs(i,1),bat_pairs(i,2)) = 0.01;
    end
    if data_Lv(1,i,2)<1e-3,     A(bat_pairs(i,2),bat_pairs(i,1)) = 0.7;
    elseif data_Lv(1,i,2)<1e-2, A(bat_pairs(i,2),bat_pairs(i,1)) = 0.5;
    elseif data_Lv(1,i,2)<5e-2, A(bat_pairs(i,2),bat_pairs(i,1)) = 0.3;
    else,                  A(bat_pairs(i,2),bat_pairs(i,1)) = 0.01;
    end
end
G = digraph(A,cellstr(bat_nms));  G.Edges.LWidths = 20*G.Edges.Weight;
nexttile;   p_NG = plot(G);
set(p_NG,'LineWidth',G.Edges.LWidths,'MarkerSize',10,'NodeLabelColor',bat_clr,'NodeColor',bat_clr,'ArrowSize',10,'NodeFontSize',15,'NodeFontWeight','bold','EdgeColor',0.5*[1 1 0]);
ave_Cv = median(data_Cv,1);
A = zeros(n_tags);
for i = 1:n_pairs
    if ave_Cv(:,i)<1e-3,     A(bat_pairs(i,1),bat_pairs(i,2)) = 0.7;
    elseif ave_Cv(:,i)<1e-2, A(bat_pairs(i,1),bat_pairs(i,2)) = 0.5;
    elseif ave_Cv(:,i)<5e-2, A(bat_pairs(i,1),bat_pairs(i,2)) = 0.3;
    else,                A(bat_pairs(i,1),bat_pairs(i,2)) = 0.01;
    end
end
G = graph(A,cellstr(bat_nms),'upper');  G.Edges.LWidths = 20*G.Edges.Weight;
nexttile;   p_NG = plot(G);
set(p_NG,'LineWidth',G.Edges.LWidths,'MarkerSize',10,'NodeLabelColor',bat_clr,'NodeColor',bat_clr,'NodeFontSize',15,'NodeFontWeight','bold','EdgeColor',0.5*[1 0.8 0.2]);

% %=== Look at the co-occurrence of proximity and landing indexes
% score_Pv = randn(size(data_Pv))*0.05;    score_Lv = randn(size(data_Lv))*0.05;
% score_Pv(data_Pv<0.001) =  1.0;     score_Lv(data_Lv<0.001) =  1.0;
% score_Pv(data_Pv<0.010) =  1.0;     score_Lv(data_Lv<0.010) =  1.0;
% score_Pv(data_Pv<0.050) =  1.0;     score_Lv(data_Lv<0.050) =  1.0;
% score_Pv(data_Pv>0.999) = -1.0;     score_Lv(data_Lv>0.999) = -1.0;
% score_Pv(data_Pv>0.990) = -1.0;     score_Lv(data_Lv>0.990) = -1.0;
% score_Pv(data_Pv>0.950) = -1.0;     score_Lv(data_Lv>0.950) = -1.0;
% scatter3(reshape(score_Lv(:,:,1),[],1),reshape(score_Lv(:,:,2),[],1),score_Pv(:),[],'filled','MarkerFaceColor','k');
% xlabel('LI1');  ylabel('LI2');  zlabel('PI');
% axis equal;

%% How do preferred locations evolve across days?

figure();   set(gcf, 'units','normalized','outerposition',[0.3 0.3 0.4 0.5]);
tiledlayout(3,n_tags,'TileSpacing','tight');
exp_f = 100;
for i = 1:n_tags
    nexttile;
    for j=1:n_sessions
        scatter3(T.centroid{j,1}{i,1}(:,1),T.centroid{j,1}{i,1}(:,2),j*ones(size(T.centroid{j,1}{i,1}(:,1))),T.samplesIn{j,1}{i,1}./T.T(j)*exp_f,'filled','MarkerFaceColor', bat_clr(i,:));
        hold on;
    end
    hold off;   xlim(r_lim(1,:)); ylim(r_lim(2,:)); grid on;    xlabel('x (m)');    ylabel('y (m)');    zlabel('Session');  xticks([]); yticks([]);
end
for i = 1:n_tags
    nexttile;
    for j=1:n_sessions
        scatter3(T.centroid{j,1}{i,1}(:,1),j*ones(size(T.centroid{j,1}{i,1}(:,1))),T.centroid{j,1}{i,1}(:,3),T.samplesIn{j,1}{i,1}./T.T(j)*exp_f,'filled','MarkerFaceColor', bat_clr(i,:));
        hold on;
    end
    hold off;   xlim(r_lim(1,:)); zlim(r_lim(3,:)); grid on;    xlabel('x (m)');    ylabel('Session');    zlabel('z (m)');  xticks([]); zticks([]);
end
for i = 1:n_tags
    nexttile;
    for j=1:n_sessions
        scatter3(j*ones(size(T.centroid{j,1}{i,1}(:,1))),T.centroid{j,1}{i,1}(:,2),T.centroid{j,1}{i,1}(:,3),T.samplesIn{j,1}{i,1}./T.T(j)*exp_f,'filled','MarkerFaceColor', bat_clr(i,:));
        hold on;
    end
    hold off;    zlim(r_lim(3,:)); ylim(r_lim(2,:)); grid on;   xlabel('Session');    ylabel('y (m)');    zlabel('z (m)');  yticks([]); zticks([]);
end

%% What is the 2D-correlation within and across bats for the 2D-occupancy maps?

%=== Plot the 2D-occupancy maps across bats and days
figure();   set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
tiledlayout(n_tags,n_sessions,'TileSpacing','tight');
for i=1:n_tags, for j=1:n_sessions
        nexttile;   histogram2('XBinEdges',edges_d{1,1},'YBinEdges',edges_d{1,2},'BinCounts',T.bat(j,i).spatial_map,'FaceColor','flat','EdgeColor','none');
        colormap(gca,custom_map(:,:,i));  %view(0,90);
        xticks([]); yticks([]);
        if i==1,title(['Session ', num2str(j)]);end
    end,end

%=== Define Map_corr(i,k,j,l): correlation between bat (i,k), sessions (j,l)
Map_corr = zeros(n_tags,n_tags,n_sessions,n_sessions);  Map_pval = zeros(n_tags,n_tags,n_sessions,n_sessions);
for i=1:n_tags,for j = 1:n_sessions,for k = 1:n_tags,for l = 1:n_sessions
                [Map_corr(i,k,j,l),Map_pval(i,k,j,l)] = corr(reshape([T.bat(j,i).spatial_map],[],1),reshape([T.bat(l,k).spatial_map],[],1));
            end,end,end,end

%=== Plot the correlation between session j+1 and j of a given bat (diagonal) or sessions j+1 from two bats (off-diagonal)
figure('units','normalized','outerposition',[0.3 0.3 0.4 0.5]);
tiledlayout(n_tags,n_tags,'TileSpacing','none');
for i =1:n_tags
    for k=1:n_tags
        nexttile;   data_corr = []; data_pval = [];
        for j = 1:n_sessions-1
            if i==k, data_corr = [data_corr; squeeze(Map_corr(i,k,j+1,j))];
                data_pval = [data_pval; squeeze(Map_pval(i,k,j+1,j))];
            else,    data_corr = [data_corr; squeeze(Map_corr(i,k,j+1,j+1))];
                data_pval = [data_pval; squeeze(Map_corr(i,k,j+1,j+1))];
            end
        end
        %text([1:n_sessions-1]',ones(n_sessions-1,1),char((data_pval<0.05).*42),'HorizontalAlignment','center','FontSize',15);
        hold on;
        if i==k, area(data_corr,'FaceColor',bat_clr(i,:));          hold off;
        elseif i<k,    area(data_corr,'FaceColor','k');             hold off;
        end
        grid on;    grid minor; ylim([0 1]);   yticks([]);    xticks([]);
    end
end

%% Show spatial preferences across days and correlation

figure('units','normalized','outerposition',[0.3 0.3 0.4 0.35]);
tiledlayout(2,n_tags,'TileSpacing','tight');
exp_f = 100;
for i = 1:n_tags
    nexttile;
    for j=1:n_sessions
        scatter3(T.centroid{j,1}{i,1}(:,1),T.centroid{j,1}{i,1}(:,2),j*ones(size(T.centroid{j,1}{i,1}(:,1))),T.samplesIn{j,1}{i,1}./T.T(j)*exp_f,'filled','MarkerFaceColor', bat_clr(i,:));
        hold on;
    end
    hold off;   xlim(r_lim(1,:)); ylim(r_lim(2,:)); grid on;    xlabel('x (m)');    ylabel('y (m)');    zlabel('Session');  xticks([]); yticks([]);
end

%=== Define Map_corr(i,k,j,l): correlation between bat (i,k), sessions (j,l)
Map_corr = zeros(n_tags,n_tags,n_sessions,n_sessions);  Map_pval = zeros(n_tags,n_tags,n_sessions,n_sessions);
for i=1:n_tags,for j = 1:n_sessions,for k = 1:n_tags,for l = 1:n_sessions
                [Map_corr(i,k,j,l),Map_pval(i,k,j,l)] = corr(reshape([T.bat(j,i).spatial_map],[],1),reshape([T.bat(l,k).spatial_map],[],1),'type','Pearson');
end,end,end,end

%=== Plot the correlation between session j+1 and j of a given bat (diagonal) or sessions j+1 from two bats (off-diagonal)
for i =1:n_tags
    nexttile;   data_corr = []; data_pval = [];
    for j = 1:n_sessions-1
        data_corr = [data_corr; squeeze(Map_corr(i,i,j+1,j))];
        data_pval = [data_pval; squeeze(Map_pval(i,i,j+1,j))];
    end
    %text([1:n_sessions-1]',ones(n_sessions-1,1),char((data_pval<0.05).*42),'HorizontalAlignment','center','FontSize',15);
    hold on;
    area(data_corr,'FaceColor',bat_clr(i,:));          hold off;
    ylim([0 1]);   %yticks([]);    xticks([]);
end

%=== Plot the correlation between session j+1 and j of a given bat (diagonal) or sessions j+1 from two bats (off-diagonal)
figure
corr_self = zeros(n_sessions-1,n_tags); corr_othr = zeros(n_sessions-1,n_tags);
for i=1:n_tags
     oth_bats = setdiff(bat_ids,i);
    for j = 1:n_sessions-1
        corr_self(j,i) = squeeze(Map_corr(i,i,j+1,j));
        corr_othr(j,i) = mean(Map_corr(i,oth_bats,j,j));
    end
    hold on;
    %plot(corr_self(:,i),'Color',bat_clr(i,:));      
    %plot(corr_othr(:,i),'Color',bat_clr(i,:),'LineStyle','--');   
end
plot(mean(corr_self,2),'k','LineWidth',3); 
plot(mean(corr_othr,2),'--k','LineWidth',3); 

ylim([0 1]);
hold off;

%% What is the difference in the NN-bat distance at takeoff before feeding and landing after feeding?

% %=== Accumulate distances around feeding times
% c = 1; ifi = [];
% for j=1:n_sessions
%     %=== Fligth table for that day
%     flight_table = sortrows(T.FLIGHTS{j,1},'smp1','ascend');
%     ifi = [ifi; diff(flight_table.t1)];
%     %=== Loop through all the feeding events
%     for k = find(flight_table.class == 'f' | flight_table.class == 'b')'
%         flight_batId = flight_table.id(k);                                                                                      % Id of the bat
%         flight_sample = flight_table.smp1(k);                                                                                   % Takeoff sample
%         f_subset_all = flight_table(flight_table.id == flight_table.id(k),:);                                                   % All Fligths from that bat
%         flight_number = find(f_subset_all.smp1 == flight_sample);                                                               % Flight number for that bat
%         d_tko(c) = T.f_nn{j,1}{flight_batId, 1}(flight_number,1);                                                               % NN-bat at takeoff
%         if flight_number<size(T.f_nn{j,1}{flight_batId, 2},1),  d_lnd(c) = T.f_nn{j,1}{flight_batId, 2}(flight_number+1,1);     % NN-bat at landing of the next flight
%         else,                                                   d_lnd(c) = NaN;end
%         c=c+1;
%     end
% end
%
% %=== Plot the distibution of difference in distance (landing-takeoff) for all flights and close-takeoff flights
% figure('units','normalized','outerposition',[0.3 0.2 0.1 0.7]);
% PlotDistr_AF_v0([d_lnd-d_tko],[1 0 0],'NN Bat distance at takeoff vs landing'); grid on;
% title(['All events, p= ' num2str(signrank(d_tko-d_lnd),3)]);
% figure();   set(gcf, 'units','normalized','outerposition',[0.4 0.2 0.1 0.7]);
% PlotDistr_AF_v0([d_lnd(d_tko<d_th)-d_tko(d_tko<d_th)],[1 0 0],'NN Bat distance at takeoff vs landing'); grid on;
% title(['Close takeoff, p= ' num2str(signrank(d_tko(d_tko<d_th)-d_lnd(d_tko<d_th)),3)]);

%% Do they fly in pairs more than what would be expected by independent events?

cpl_fly = reshape([T.pair(:,:).p_fly2],2,[],n_pairs);
p_max =  max(cpl_fly(:));

figure;
scatter(reshape(cpl_fly(2,:,:),[],1),reshape(cpl_fly(1,:,:),[],1),'k','filled');    hold on;
plot([0 p_max],[0 p_max],'--k');  hold off;
xlabel('P coupled fly (independent)');  ylabel('P coupled fly (real)');
axis equal
xlim([0 p_max]);  ylim([0 p_max]);

figure('units','normalized','outerposition',[0 0.3 1 0.3]);
tiledlayout(1,n_pairs,'TileSpacing','none');
for i=1:n_pairs
    nexttile;   scatter(cpl_fly(2,:,i),cpl_fly(1,:,i),'k','filled');    hold on;
    plot([0 p_max],[0 p_max],'--k');
    patch([0,p_max,0]',[0,p_max,p_max]',bat_clr(bat_pairs(i,1),:),'FaceAlpha',0.3,'LineStyle','none');
    patch([0,p_max,p_max]',[0,p_max,0]',bat_clr(bat_pairs(i,2),:),'FaceAlpha',0.3,'LineStyle','none');  hold off;
    if i==1,xlabel('P coupled fly (independent)');  ylabel('P coupled fly (real)');else,xticks([]);yticks([]);end
    axis equal; xlim([0 p_max]);  ylim([0 p_max]);
end

%% Look at difference in the number of flights

% %=== Maximum difference in the number of flights
% data_max = reshape([T.pair(:,:).max_flight_diff],[],n_pairs);
% figure; set(gcf, 'units','normalized','outerposition',[0.3 0.3 0.1 0.3]);
% PlotDistr_AF_v0(data_max',jet(n_pairs),'Max Flight Difference');
%
% %=== Median difference in the number of flights
% data_med = reshape([T.pair(:,:).med_flight_diff],[],n_pairs);
% figure; set(gcf, 'units','normalized','outerposition',[0.3 0.3 0.1 0.3]);
% PlotDistr_AF_v0(data_med',jet(n_pairs),'Med Flight Difference');
%
% %=== Sign percentage for the difference in the number of flights
% data_sig = reshape([T.pair(:,:).sig_flight_diff],[],n_pairs);
% figure; set(gcf, 'units','normalized','outerposition',[0.3 0.3 0.1 0.3]);
% PlotDistr_AF_v0(data_sig',jet(n_pairs),'Sign Flight Difference');
%
% figure();   set(gcf, 'units','normalized','outerposition',[0.3 0 0.1 1]);
% tiledlayout(n_pairs,1,'TileSpacing','tight');
% for i=1:n_pairs
%     nexttile;
%     scatter(data_med(:,i),data_sig(:,i),'filled');
%     xlim([-90 90]); ylim([0 1]);
%     rectangle('Position',[-10 0 20 1],'FaceColor',[0 0 0 0.2],'EdgeColor','none');
%     title([bat_nms(bat_pairs(i,1),:) '-' bat_nms(bat_pairs(i,2),:)]);
% end

%% Test the relationship between weight and number of flights, time on feeder, etc...

% %=== Test the relationship between total number of flights and weight
% flights_hour = reshape([T.bat(:,:).f_num],[],n_tags)./(T.T/Fs/3600);
% Test_Flights_Weight_AF_v0;
%
% %=== Test the relationship between number of flights to the feeder and weight
% feeds_hour = reshape([T.bat(:,:).f_num],[],n_tags)./(T.T/Fs/3600).*reshape([T.bat(:,:).f_flights],[],n_tags);
% Test_Feeds_Weight_AF_v0;
%
% %=== Test the relationship between time spent close to the feeder and weight


%% Visualize the dataset

Folder = cd;    FileList = dir(fullfile(Folder, '**', 'Analyzed_Behavior_*'));
Multi_Day = struct([]); durations = []; interflight = [];   stop2start = [];
for i = 1:length(FileList)
    cd(FileList(i).folder);
    load(FileList(i).name);
    path = strsplit(FileList(i).folder,'\');
    Multi_Day(i).T = T;
    Multi_Day(i).n_tags = size(pref_location,1);
    Multi_Day(i).group_name = Group_name;
    Multi_Day(i).Dataset = path{1,4};
    Multi_Day(i).f_nn = f_nn;
    Multi_Day(i).d_th = d_th;                   % Threshold distance for later classifications
    for j = 1:size(f_smp,1)
        durations = [durations; (f_smp{j, 2}-f_smp{j, 1})/100];
        interflight = [interflight; diff(f_smp{j, 1})/100];
        stop2start = [stop2start;   diff(f_smp{j, 1})/100-(f_smp{j, 2}(1:end-1)-f_smp{j, 1}(1:end-1))/100];
    end
    Multi_Day(i).lt_d = vecnorm(FLIGHTS.r1(2:end,:)-FLIGHTS.r2(1:end-1,:),2,2);     % Get the takeoff to previous landing distance
    Multi_Day(i).lt_d = Multi_Day(i).lt_d(~diff(FLIGHTS.id),:);                     % Remove spurious distances across different bat ids
    disp([num2str(length(FileList)-i),' remaining sessions to load...']);
end
cd(Folder);
clearvars -except Multi_Day durations interflight stop2start

T = struct2table(Multi_Day);
T.Dataset = categorical(T.Dataset);
n_sessions = size(Multi_Day,2);                                                                                 % Number of sessions
r_lim = [-2.8 2.8; -2.6 2.6; 0 2.30];                                                                           % Room boundaries
diagonal = vecnorm(diff(r_lim,1,2));                                                                            % Room diagonal
d_th = unique(T.d_th);                                                                                          % Threshold distance for later classifications

%%

%=== Look at Takeoff to previous landing distance, duration and interflight
figure('units','normalized','outerposition',[0.3 0.3 0.2 0.6]);
tiledlayout(3,1,'TileSpacing','tight');
lt_d = vertcat(T.lt_d{:,:});
nexttile;   histogram(lt_d,[0:0.01:1],'edgecolor','none','FaceColor','k');          xlabel('Takeoff to previous landing distance (m)');     yticks([]);   set(gca,'fontsize',16);
title(['Median = ', num2str(median(lt_d),2), ' m, 95% of data <', num2str(prctile(lt_d,95),2), ' m']);
nexttile;   histogram(durations,[0:0.1:10],'edgecolor','none','FaceColor','k');     xlabel('Duration (s)');                                 yticks([]);   set(gca, 'fontsize',16);
nexttile;   histogram(interflight,10.^[0:0.01:3],'edgecolor','none','FaceColor','k'); xlabel('Interflight (s)');                            yticks([]);   set(gca,'XScale','log','fontsize',16);
%nexttile;   histogram(interflight,[0:1:100],'edgecolor','none','FaceColor','k');    xlabel('Start 2 start (s)');                            yticks([]);

%=== All sessions
cmap = hsv(14);
cmap = cmap(randperm(size(cmap, 1)), :);
b = bar(T.n_tags);
b.FaceColor = 'flat';
b.CData(find(cell2mat(T.group_name) == 'D'),:) = [.5 0 .5].*ones(length(find(cell2mat(T.group_name) == 'D')),1);
hold on;
gscatter([1:size(T,1)]',7.2*ones(size(T,1),1),T.Dataset,cmap,[],[],'doleg','off');
hold off;       xlabel('Session');  ylabel('Number of bats');

qqplot(log10(interflight));

%=== Session duration
figure; histogram(T.T./(100*3600)); xlabel('Session duration (hrs)');
figure; histogram(T.n_tags); xlabel('Number of bats');

%=== Overall landings
nnb_distances = [];
for j=1:n_sessions
    for i=1:size(T.f_nn{j,1},1)
        nnb_distances = [nnb_distances; T.f_nn{j,1}{i, 2}];
    end
end

%=== Visualize all landings
figure; set(gcf, 'units','normalized','outerposition',[0.3 0.3 0.1 0.4]);
tiledlayout(2,1,'TileSpacing','none');
edges_dist = 10.^(-2:0.01:log10(diagonal));
nexttile;   histogram(nnb_distances,edges_dist,'Normalization','probability','edgecolor','none','FaceColor','k','FaceAlpha',1);  
set(gca, 'XScale', 'log');  xticks([]); hold on;    plot([.6 .6],ylim,'r',[.9 .9],ylim,'b');  hold off; xlim([0.01 10]);                     ylim('tight');  yticks([]); title('Bat distance at landing');
nexttile;   histogram(nnb_distances,edges_dist,'Normalization','cdf','edgecolor','none','FaceColor','k','FaceAlpha',1);  
set(gca, 'XScale', 'log');  hold on;    plot([.6 .6],ylim,'r',[.9 .9],ylim,'b');  hold off; xlim([0.01 10]); xlabel('NN-bat distance (m)');  ylim('tight');  ylabel('Fraction');

%=== Visualize with log scale
edges_dist = 10.^[-2:0.02:log10(diagonal)];
edge_dist_cent = edges_dist(1:end-1)+ diff(edges_dist)/2;
[p_distance,~] = histcounts(nnb_distances,edges_dist,'Normalization','probability');
figure; set(gcf, 'units','normalized','outerposition',[0.1 0.3 0.2 0.3]);
tiledlayout(1,2,'TileSpacing','none');
nexttile;   histogram(nnb_distances,edges_dist,'Normalization','probability','edgecolor','none','FaceColor','k');
set(gca, 'XScale', 'log');  xlabel('NN-bat distance at landing (m)');   yticks([]);
hold on;    y1=get(gca,'ylim');  hold on; plot([d_th d_th],y1); hold off;
nexttile;   h = radialHist_AF_v0(p_distance,edge_dist_cent,[0,0.5],[0 360],[0 0 0]); h.CDataMode = 'auto'; axis equal; xlabel('NN-bat distance at landing (m)');    yticks([]);

%=== Visualize with linear scale
edges_dist = [0.01:0.01:0.50];
edge_dist_cent = edges_dist(1:end-1)+ diff(edges_dist)/2;
[p_distance,~] = histcounts(nnb_distances,edges_dist,'Normalization','probability');
figure('units','normalized','outerposition',[0.1 0.3 0.2 0.3]);
tiledlayout(1,2,'TileSpacing','none');
nexttile;   histogram(nnb_distances,edges_dist,'Normalization','probability','edgecolor','none','FaceColor','k');
xlabel('NN-bat distance at landing (m)');   yticks([]);
hold on;    y1=get(gca,'ylim');  hold on; plot([d_th d_th],y1); hold off;
nexttile;   h = radialHist_AF_v0(p_distance,edge_dist_cent,[0,0.5],[0 360],[0 0 0]); h.CDataMode = 'auto'; axis equal; xlabel('NN-bat distance at landing (m)');    yticks([]);

%=== Visualize CDF
sorted_dist = sort(nnb_distances);
sorted_dist(ceil(end)/2)
plot(sorted_dist,[1:numel(sorted_dist)]/numel(sorted_dist),'k','LineWidth',3);
set(gca,'XScale','Log'); xlim([0.01 10]); 

%=== Fit values below the peak with a Half-Normal
figure;
h = histogram(nnb_distances(nnb_distances<0.5),edges_dist);
[maxcount, whichbin] = max(h.Values);
mu = mean([edges_dist(whichbin),edges_dist(whichbin+1)]);
x = -nnb_distances(nnb_distances<mu)+mu;
pd = fitdist(x,'HalfNormal');
histfit(x,50,'HalfNormal');
title(['Mean = ', num2str(mu,3), ' Std = ', num2str(pd.sigma,3)]);
disp(['Optimal distance threshold: ', num2str(mu+2*pd.sigma,3), ' m']);

%% Concateate all flights in a single table

Folder = cd;    FileList = dir(fullfile(Folder, '**', 'Analyzed_Behavior_*'));
FLIGHTS_all = table; perc_fast = [];
for i = 1:length(FileList)
    cd(FileList(i).folder);
    load(FileList(i).name);
    path = strsplit(FileList(i).folder,'\');
    FLIGHTS.g_name = repmat(Group_name,size(FLIGHTS,1),1);
    FLIGHTS.session = repmat(path{1,5},size(FLIGHTS,1),1);
    FLIGHTS.T = T*ones(size(FLIGHTS,1),1);
    FLIGHTS.nn = cell2mat(f_nn(:,2));
    FLIGHTS_all = [FLIGHTS_all; FLIGHTS];
    
    %=== Quantify percentage of flights during which there is another bat taking off from the landing spot
    FLIGHTS = sortrows(FLIGHTS,'t1');
    counter = 0;
    for j=1:size(FLIGHTS,1)
        candidate_t1 = FLIGHTS.t1(FLIGHTS.t1>FLIGHTS.t1(j) & FLIGHTS.t1<FLIGHTS.t2(j));
        candidate_r1 = FLIGHTS.r1(FLIGHTS.t1>FLIGHTS.t1(j) & FLIGHTS.t1<FLIGHTS.t2(j),:);
        for nn=1:numel(candidate_t1)
            if vecnorm(candidate_r1(nn,:)-FLIGHTS.r2(j))<0.6
                counter = counter+1;
            end 
        end
    end
    perc_fast = [perc_fast;counter/size(FLIGHTS,1)];
    
    disp([num2str(length(FileList)-i),' remaining sessions to load...']);
end
cd(Folder);
clearvars -except FLIGHTS_all perc_fast

bat_clr = lines(5);
r_lim = [-2.9 2.9; -2.6 2.6; 0 2.30];                                                                           % Room boundaries


%=== Show all landing spots (FOR FIGURE 1)
[centroids_all,subsample] = datasample(FLIGHTS_all.r2,1e4,'Replace',false);
category_all = FLIGHTS_all.class(subsample);
view_angle = [-37.5,30;0,90;0,0;90,0];
figure; set(gcf, 'units','normalized','outerposition',[0.1 0.3 0.1 0.21]);
cond = ~(category_all=='f' | category_all=='b');
i=1;
scatter3(centroids_all(cond,1),centroids_all(cond,2),centroids_all(cond,3),3,'MarkerFaceColor', 'k','MarkerEdgeColor','none','MarkerFaceAlpha',.2);
view(view_angle(i,:));
xlim(r_lim(1,:)); ylim(r_lim(2,:)); zlim(r_lim(3,:));
axis equal; set(gca,'XTickLabels',[],'YTickLabels',[],'ZTickLabels',[]);
xlabel('x','fontsize',16);ylabel('y','fontsize',16);zlabel('z','fontsize',16);
xticks(r_lim(1,:)); yticks(r_lim(2,:)); zticks(r_lim(3,:));
sgtitle(['Landing spots from n = ', num2str(size(centroids_all,1)), ' flights']);

ptCloud = pointCloud(centroids_all(cond,:));
pcshow(ptCloud);
indices = pcbin(ptCloud,round(diff(r_lim,1,2)/0.05)');
occupancyGrid = cellfun(@(c) ~isempty(c), indices);
v = volshow(occupancyGrid);

%=== Keep only D or F bats
D_FLIGHTS = FLIGHTS_all(FLIGHTS_all.g_name == 'D' & FLIGHTS_all.id < 6,:);
F_FLIGHTS = FLIGHTS_all(FLIGHTS_all.g_name == 'F' & FLIGHTS_all.id < 6,:);

%=== Define categorical variables
D_FLIGHTS.session = categorical(cellstr(D_FLIGHTS.session));
F_FLIGHTS.session = categorical(cellstr(F_FLIGHTS.session));
D_FLIGHTS.id = categorical(D_FLIGHTS.id);
F_FLIGHTS.id = categorical(F_FLIGHTS.id);

D_tmp = varfun(@mean,D_FLIGHTS,'InputVariables','T','GroupingVariables','session');
F_tmp = varfun(@mean,F_FLIGHTS,'InputVariables','T','GroupingVariables','session');

%=== Calculate Flights number in 5 min chunks
D_Flights_per_5min = varfun(@(x) histcounts(x,[0:300:3300]),D_FLIGHTS,'InputVariables','t1','GroupingVariables',{'id','session'});
F_Flights_per_5min = varfun(@(x) histcounts(x,[0:300:3300]),F_FLIGHTS,'InputVariables','t1','GroupingVariables',{'id','session'});

%=== Show nn-distance distribution
figure;
edges_dist = 10.^[-2:0.01:log10(8.2)];
histogram(D_FLIGHTS.nn,edges_dist,'edgecolor','none','FaceColor','k','FaceAlpha',1);  hold on;
plot([.5 .5],ylim,'r'); hold off;
set(gca, 'XScale', 'log');  xlabel('NN-bat distance at landing (m)');   yticks([]);
set(gca,'fontsize',16); box off;
figure;
area(sort(D_FLIGHTS.nn),[1:numel(D_FLIGHTS.nn)]./numel(D_FLIGHTS.nn),'FaceColor','k');  hold on;
plot([.5 .5],ylim,'r'); hold off;
set(gca, 'XScale', 'log');  xlabel('NN-bat distance at landing (m)');   
set(gca,'fontsize',16); box off;

%=== Show all landing spots
centroids_all = D_FLIGHTS.r2;
view_angle = [-37.5,30;0,90;0,0;90,0];
figure; set(gcf, 'units','normalized','outerposition',[0.1 0.3 0.2 0.43]);
tiledlayout(2,2,'TileSpacing','tight');
for i=1:4
    nexttile;
    scatter3(centroids_all(:,1),centroids_all(:,2),centroids_all(:,3),1,'filled','MarkerFaceColor', 'k');
    view(view_angle(i,:));   
    xlim(r_lim(1,:)); ylim(r_lim(2,:)); zlim(r_lim(3,:));
    axis equal; set(gca,'XTickLabels',[],'YTickLabels',[],'ZTickLabels',[]); 
    xlabel('x','fontsize',16);ylabel('y','fontsize',16);zlabel('z','fontsize',16); 
end
sgtitle([num2str(size(D_tmp,1)),' Sessions, Tot Rec Time: ', num2str(sum(D_tmp.mean_T/100)/(3600*24),2),' days']);

%=== Plot
figure('units','normalized','outerposition',[.3 0 .1 1]);
for i =1:5
    data = D_Flights_per_5min.Fun_t1(D_Flights_per_5min.id == num2str(i),:);
    nexttile;   plot([5:5:55],data','.k');    hold on;    plot([5:5:55],mean(data),'Color',bat_clr(i,:),'LineWidth',5);   hold off;
    xlabel('Time (min)');   ylabel('Number of flights/5min');
end
sgtitle('D bats');

figure('units','normalized','outerposition',[.3 0 .1 1]);
for i =1:5
    data = F_Flights_per_5min.Fun_t1(F_Flights_per_5min.id == num2str(i),:);
    nexttile;   plot([5:5:55],data','.k');    hold on;    plot([5:5:55],mean(data),'Color',bat_clr(i,:),'LineWidth',5);   hold off;
    xlabel('Time (min)');   ylabel('Number of flights/5min');
end
sgtitle('F bats');

%% === Test for Exponential disribution of the interflight intervals (after 10 min)
D_Flights_exp = varfun(@(x) custom(x),D_FLIGHTS,'InputVariables','t1','GroupingVariables',{'id','session'});
F_Flights_exp = varfun(@(x) custom(x),F_FLIGHTS,'InputVariables','t1','GroupingVariables',{'id','session'});
disp(nnz(D_Flights_exp.Fun_t1)/numel(D_Flights_exp.Fun_t1)*100);
disp(nnz(F_Flights_exp.Fun_t1)/numel(F_Flights_exp.Fun_t1)*100);

%=== Calculate interflight difference distribution
D_Flights_alldiff = varfun(@(x) histcounts(diff(x(x>600)),[0:1:180],'Normalization','probability'),D_FLIGHTS,'InputVariables','t1','GroupingVariables',{'session'});
F_Flights_alldiff = varfun(@(x) histcounts(diff(x(x>600)),[0:1:180],'Normalization','probability'),F_FLIGHTS,'InputVariables','t1','GroupingVariables',{'session'});

plot(mean(D_Flights_alldiff.Fun_t1));

% https://stats.stackexchange.com/questions/69896/measuring-correlation-of-point-processes

function f = custom(x)
T = 600;
if nnz(x>T)>4
    f = lillietest(diff(x(x>T)),'Distribution','exponential');
else
    f = NaN;
end
end







