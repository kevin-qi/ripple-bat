%% Look at landing/taking off distributions

%---Params and matrices
n_tags = 6;
bat_pairs = nchoosek(1:n_tags,2);                                                                               %Bat Pairs
r = conc.r;
T = length(r);
bflying = conc.b;

%---Calculate bat distances
bat_dist = zeros(T,length(bat_pairs));
for i = 1:length(bat_pairs)
    bat_dist(:,i) = vecnorm(r(:,:,bat_pairs(i,1))-r(:,:,bat_pairs(i,2)),2,2);
end

%---Calculate egocentric coordinates
r_ego = zeros(T,3,n_tags,n_tags);
for i = 1:n_tags
    for j =1:n_tags
        r_ego(:,:,i,j) = r(:,:,j)-r(:,:,i);
    end
end

%---Extract flights takeoff and landing
for i = 1:n_tags
    f_num(i) = nnz(diff(bflying(:,i))>0);
    f_smp{i,1} = find(diff(bflying(:,i))>0)+1;
    f_smp{i,2} = find(diff(bflying(:,i))<0);
end

%---Calculate distances at taking off 
nnbtk = zeros(sum(f_num),1);
c = 1;
for i = 1:n_tags
    i_pairs = any(ismember(bat_pairs,i),2);
    for nn=1:f_num(i)
        nnbtk(c) = min(bat_dist(f_smp{i,1}(nn),i_pairs));
        c = c+1;
    end
end

%---Calculate distances at landing 
nnbld = zeros(sum(f_num),1);
c = 1;
for i = 1:n_tags
    i_pairs = any(ismember(bat_pairs,i),2);
    for nn=1:f_num(i)
        nnbld(c) = min(bat_dist(f_smp{i,2}(nn),i_pairs));
        c = c+1;
    end
end

figure();       set(gcf, 'units','normalized','outerposition',[0 0.25 0.34 0.4]);
tiledlayout(1,2,'TileSpacing','tight','Padding','compact');
edges_t = 10.^[-2:0.05:1];
nexttile; histogram(nnbtk,edges_t,'edgecolor','none','FaceColor','k');  set(gca,'XScale','log');  xlabel('Distance at take off (m)'); ylabel('Counts');
nexttile; histogram(nnbld,edges_t,'edgecolor','none','FaceColor','k');  set(gca,'XScale','log');  xlabel('Distance at landing (m)'); ylabel('Counts');

histogram([nnbtk],edges_t,'edgecolor','none','FaceColor','k','Normalization','cdf');  set(gca,'XScale','log');  xlabel('Distance at take off (m)'); ylabel('Cumulative fraction');

%---FIG: Landing in ego-centric coordinates, grand-histogram
sphere_r = 2;
edges_s = {-sphere_r:(2*sphere_r)/100:sphere_r -sphere_r:(2*sphere_r)/100:sphere_r};                 %Edges for density histogram
l_point = [];
for i=1:n_tags
    for j=1:n_tags
        if j ~= i
            for n = 1:f_num(j)
                l_point = cat(1,l_point,squeeze(r_ego(f_smp{j,1}(n,1),:,i,j)));
            end
            
        end
    end
end
figure();       set(gcf, 'units','normalized','outerposition',[0 0.25 0.34 0.4]);
tiledlayout(1,3,'TileSpacing','tight','Padding','compact');
ax_nms = ['x';'y';'z'];
for i=1:3
    sgtitle('Taking off Spot Distributions');
    ax_idx = setdiff(1:3,i);
    nexttile;   hist3(l_point(:,ax_idx),'edges',edges_s,'CdataMode','auto','edgecolor','none','FaceColor','interp');
    xlabel(ax_nms(ax_idx(1)));    ylabel(ax_nms(ax_idx(2)));    xlim([-sphere_r sphere_r]); ylim([-sphere_r sphere_r]);    axis square;    hold on;    
    line([0,0],[1,-1],'Color','white','LineWidth',0.1); line([-1,1],[0,0],'Color','white','LineWidth',0.1); hold off;    view(0,-90);
    colormap('jet');
end
