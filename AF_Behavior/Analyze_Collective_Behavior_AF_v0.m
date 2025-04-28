%% Script for the analysis of collective behavior (September 2021)
% Requires extracted data from Ciholas recordings

%Load data
extracted_CDPfile = dir(fullfile(cd, '*extracted_*_cdp_*'));          load(extracted_CDPfile.name);
batdate = extracted_CDPfile.name(11:16);

%Parameters and options
n_tags = CDPmtdata.tags;
r_lim = [-2.8 2.8; -2.6 2.6; 0 2.30];                                                                           %Room boundaries
edges_d = {r_lim(1,1):(r_lim(1,2)-r_lim(1,1))/10:r_lim(1,2) r_lim(2,1):(r_lim(2,2)-r_lim(2,1))/10:r_lim(2,2)};  %Edges for density histogram
Fs = 100;                                                                                                       %Sampling frequency (Hz) for common time
bat_nms = ['Dai'; 'Den'; 'Dia'; 'Dor'; 'Dum'; 'Im1'; 'Im2'];    bat_nms = bat_nms(1:n_tags,:);                  %Bat Names
bat_pairs = nchoosek(1:n_tags,2);                                                                               %Bat Pairs
bat_pair_nms = [bat_nms(bat_pairs(:,1),:), '-'.*ones(length(bat_pairs),1), bat_nms(bat_pairs(:,2),:)];          %Bat Pairs Names
bat_clr = lines(n_tags);                                                                                        %Bat Colors
v_th = 0.5;                                                                                                     %Velocity threshold (m/s) for flight segmentation
alpha_clus = 1;                                                                                                 %clustering parameter

options.show_fig = 1;
options.save_data = 1;
options.use_r_corr = 1;
options.savemovie = 0;
options.clusterFl = 0;

%Custom graded colormap
for i = 1:n_tags
    for j = 1:3
        custom_map(:,j,i) = linspace(1,bat_clr(i,j))';
    end
end

%Create analysis folder for storing the results
if options.save_data
    analysis_directory=fullfile(pwd,['Analysis_',datestr(now, 'yymmdd_HHMM')]);
    if ~exist(analysis_directory,'dir')
        mkdir(analysis_directory);
    end
end

%Shift Ciholas Time by 1.05s--> 1st TTL will be aligned to 3s M9 TTL
for i = 1:n_tags
    tag_data{1,i}(:,8) = tag_data{1,i}(:,8)+1.05;
    tag_data_filt{1,i}(:,8) = tag_data_filt{1,i}(:,8)+1.05;
    tag_ac_data{1,i}(:,8) = tag_ac_data{1,i}(:,8)+1.05;
end

%% Calculate evenly sampled kinematic variables (r,v,a)

% As defined below, time 0 corresponds to the first 3s-TTL of the M9 (and so the first camera frame)
t = [CDPmtdata.TTL_times(1):1/Fs:CDPmtdata.TTL_times(end)]';       %Evenly sampled time vector from first to last TTL.
T = length(t);                                                     %Number of time samples
r = zeros(T,3,n_tags);                                             %3D position vector (sample,dim,id)

%Interpolate position at evenly spaced time points
for i = 1:n_tags
    r(:,:,i) =  interp1(tag_data_filt{1,i}(:,8), tag_data_filt{1,i}(:,[3:5]),t,'nearest','extrap'); %DO NOT use spline interpolation!
end

%Interpolate acceleration at evenly spaced time points
if exist('tag_ac_data')
    a = zeros(T,3,n_tags);
    for i = 1:n_tags
        a(:,:,i) =  interp1(tag_ac_data{1,i}(:,8), tag_ac_data{1,i}(:,[3:5]),t,'nearest','extrap'); %DO NOT use spline interpolation!
    end
    a_abs = squeeze(vecnorm(a,2,2));    %Modulus
    a_flt = bandpass(a_abs,[7 9],100);  %Filtered at the wing-beat frequency
end

% % Uncomment to control interpolation quality
% ax(1) = subplot(311);  plot(tag_data{1,4}(:,8),tag_data{1,4}(:,3)','.',t,r(:,1,4),'-');
% ax(2) = subplot(312);  plot(tag_data{1,4}(:,8),tag_data{1,4}(:,4)','.',t,r(:,2,4),'-');
% ax(3) = subplot(313);  plot(tag_data{1,4}(:,8),tag_data{1,4}(:,5)','.',t,r(:,3,4),'-');
% linkaxes(ax,'x');
% ax(1) = subplot(311);  plot(tag_ac_data{1,1}(:,8),tag_ac_data{1,1}(:,3)','.',t,a(:,1,1),'-');
% ax(2) = subplot(312);  plot(tag_ac_data{1,1}(:,8),tag_ac_data{1,1}(:,4)','.',t,a(:,2,1),'-');
% ax(3) = subplot(313);  plot(tag_ac_data{1,1}(:,8),tag_ac_data{1,1}(:,5)','.',t,a(:,3,1),'-');
% linkaxes(ax,'x');

%Calculate velocity and 2D-direction of motion
v = diff(r,1,1)./diff(t); v=cat(1,zeros(1,3,n_tags),v);   v_abs = squeeze(vecnorm(v,2,2));
angle = squeeze(heading_angle(v(:,1,:),v(:,2,:)));

%Detect flight epochs based on wing-beat signal
if exist('tag_ac_data')
    for i = 1:n_tags
        [up,lo] = envelope(a_flt(:,i)+normrnd(0,1e-3,length(a_flt),1),10,'peak');   %Envelope of the acceleration signal (noise addedd to avoid problems with splining)
        env = normalize(up - lo,'range');           %Amplitude of the envelope
        env_th = otsuthresh(histcounts(env));       %Threshold (based on Otsu method). Can be set at 0.35
        wBeats(:,i) = movsum(env>env_th,2*Fs)>Fs/5; %Euristic criterion for flight detection
%                         ax(i) = subplot(n_tags,1,i);  area(t,wBeats(:,i)*1,'FaceAlpha',0.3,'LineStyle','none');  hold on;
%                         plot(t,normalize(v_abs(:,i),'range'));
%                         plot(t,r(:,1,i),t,r(:,2,i));
%                         plot(t,normalize(a_flt(:,i),'range',[-1 1]));
%                         plot(t,normalize(movsum(env>env_th,2*Fs),'range'));
%                         hold off;
    end
%               linkaxes(ax,'x');
end

%% Correct position by median filtering when the bat is not flapping its wings

if exist('tag_ac_data')
    stat_periods = repmat(~wBeats,1,1,3);    stat_periods = permute(stat_periods,[1 3 2]);
    r_stat = r;
    r_stat(~stat_periods) = nan;
    
    for i = 1:n_tags
        r_stat(:,:,i) =  smoothdata(r_stat(:,:,i),1,'movmedian',Fs*5,'omitnan');
    end
    r_stat(~stat_periods) = nan;
    
    % % Uncomment to control filtering quality
    % ax(1) = subplot(311);  plot(tag_data{1,4}(:,8),tag_data{1,4}(:,3)',t,r_stat(:,1,4));
    % ax(2) = subplot(312);  plot(tag_data{1,4}(:,8),tag_data{1,4}(:,4)',t,r_stat(:,2,4));
    % ax(3) = subplot(313);  plot(tag_data{1,4}(:,8),tag_data{1,4}(:,5)',t,r_stat(:,3,4));
    % linkaxes(ax,'x');
    
    %Substitute median filtered data when bat is not flapping its wings
    r_stat = fillmissing(r_stat,'constant',0,1);
    r_corr = r_stat.*stat_periods+r.*(~stat_periods);
    r_corr = smoothdata(r_corr,1,'loess',Fs*1);     %DO NOT USE lowess!!
    
    % ax(1) = subplot(311);  plot(tag_data{1,4}(:,8),tag_data{1,4}(:,3)',t,r_corr(:,1,4));
    % ax(2) = subplot(312);  plot(tag_data{1,4}(:,8),tag_data{1,4}(:,4)',t,r_corr(:,2,4));
    % ax(3) = subplot(313);  plot(tag_data{1,4}(:,8),tag_data{1,4}(:,5)',t,r_corr(:,3,4));
    % linkaxes(ax,'x');
    
else
    %Aggressive median filtering on original data
    tag_data_stat = tag_data;
    r_stat = zeros(T,3,n_tags);
    for i = 1:n_tags
        tag_data_stat{1,i}(:,[3:5]) = smoothdata(tag_data_filt{1,i}(:,[3:5]),1,'movmedian',Fs*3);
        r_stat(:,:,i) =  csaps(tag_data_stat{1,i}(:,8), tag_data_stat{1,i}(:,[3:5])', 1, t)';
    end
    %Calculate xy-velocity and smooth
    v_filt = smoothdata(squeeze( vecnorm( v(:,1:2,:),2,2)), 1, 'movmedian', Fs*3);
    %Substitute median filtered data when velocity is less than threshold
    stat_periods = repmat((v_filt < v_th),1,1,3);    stat_periods = permute(stat_periods,[1 3 2]);
    r_corr = r_stat.*stat_periods+r.*(~stat_periods);
    r_corr = smoothdata(r_corr,1,'loess',Fs*1);
    r_corr_1 = smoothdata(r_corr,1,'loess',Fs*1);
end

%% Plot Position data (raw and corrected)

for j = 1:3
    figure();   set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    tiledlayout(n_tags,1,'TileSpacing','tight');
    for i=1:n_tags
        ax(i) = nexttile;   %ax(i) = subplot(n_tags,1,i);
        plot(tag_data{1, i}(:,8), tag_data{1, i}(:,2+j),t,r_corr(:,j,i));
        sgtitle(['D-' num2str(j)]);   ylim(r_lim(j,:));   legend('raw','corrected');
        
    end
    linkaxes(ax,'x');    xlabel('Time(s)');
end

%% Correct position, velocity and angle

if options.use_r_corr
    r_old = r;
    r = r_corr; v_old = v_abs;
    v = diff(r,1,1)./diff(t); v=cat(1,zeros(1,3,n_tags),v);   v_abs = squeeze(vecnorm(v(:,1:3,:),2,2));
    angle = squeeze(heading_angle(v(:,1,:),v(:,2,:)));
end

%% Plot velocity (raw and corrected)

figure();   set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
tiledlayout(n_tags,1,'TileSpacing','tight');
for i=1:n_tags
    ax(i) = nexttile;   %ax(i) = subplot(n_tags,1,i);
    plot(t,v_old(:,i),'.');    hold on;     sgtitle('v');
    plot(t,v_abs(:,i),'.');    hold off;    legend('raw','corrected');
end
linkaxes(ax,'x');    xlabel('Time(s)');

%% Flight segmentation

%Detect flight starts and stops by using risetime and falltime
% !! Be particularly careful with the levels of risetime and falltime, these
%are fundamental in order to exclude stationary tails in the flight
bflying = zeros(T,n_tags); f_num = zeros(n_tags,1);   f_smp = cell(n_tags,2);
for i = 1:n_tags
    [bflying(:,i),f_num(i),f_smp{i,1},f_smp{i,2}] = FlightSegm_AF_v0(v_abs(:,i),v_th,Fs);
end

%Define staionary periods when the bat is not flying
stat_periods = repmat(~bflying,1,1,3);   stat_periods = permute(stat_periods,[1 3 2]);
r_qt = r;   r_qt(~stat_periods)=nan;     angle(squeeze(stat_periods(:,1,:))) = nan;

%% Flight clustering

if options.clusterFl
    for i = 1:n_tags
        f_clus(i) = FlightClus_AF_v1(squeeze(r(:,:,i)),bflying(:,i),'Alpha',alpha_clus,'Frechet',0,'Points',6);
    end
    for i = 1:n_tags
        f_clus(i).name = bat_nms(i,:);
    end
else
    f_clus=[];
end

%% Bat landing on a bat and social flights

for i = 1:n_tags
    if f_num(i)
        c = 0;
        f_land_table{i,1} = zeros(f_num(i),n_tags);
        for n =f_smp{i,2}'
            c = c+1;
            f_land_table{i,1}(c,:) = squeeze(vecnorm(r(n,:,i)-r(n,:,:),2,2)<0.25)';
            social_flag{i,1}(c,:) = sum(f_land_table{i,1}(c,:),2)>1;
        end
    else
        f_land_table{i,1} = [];
    end
end

%% Ego-centric vectors and interbat distances

r_ego = zeros(T,3,n_tags,n_tags);
for i = 1:n_tags
    for j =1:n_tags
        r_ego(:,:,i,j) = r(:,:,j)-r(:,:,i);
    end
end

bat_dist = zeros(T,length(bat_pairs));
for i = 1:length(bat_pairs)
    bat_dist(:,i) = vecnorm(r(:,:,bat_pairs(i,1))-r(:,:,bat_pairs(i,2)),2,2);
end

%% Probabilities and angles when flying together

p_fly1 = sum(bflying,1)./T;
heading_diff = nan(T,length(bat_pairs));
for i = 1:length(bat_pairs)
    %Probabilities
    p_fly2(i,1) = sum(bflying(:,bat_pairs(i,1)).*bflying(:,bat_pairs(i,2)),1)/T;
    p_fly2(i,2) = p_fly1(bat_pairs(i,1))*p_fly1(bat_pairs(i,2));
    
    %Angles
    cpled = find(bflying(:,bat_pairs(i,1)).*bflying(:,bat_pairs(i,2)));
    heading_diff(cpled,i) = rad2deg(angdiff(angle(cpled,bat_pairs(i,1)),angle(cpled,bat_pairs(i,2))));
end

%% Proximity index and social network graph

n_rep = 100;
d_th = 0.25;
p_idx = zeros(n_rep+1,length(bat_pairs));   r_sh = zeros(T,3);  p_val = zeros(length(bat_pairs),1);
A = zeros(n_tags);

for i = 1:length(bat_pairs)
    p_idx(1,i) = nnz(bat_dist(:,i)<d_th)/T;
    frames_to_shift = randi([60*Fs T-60*Fs],1,n_rep);
    for n = 1:n_rep
        r_sh = circshift(r(:,:,bat_pairs(i,2)),frames_to_shift(n),1);
        p_idx(n+1,i) = nnz(vecnorm(r(:,:,bat_pairs(i,1))-r_sh,2,2)<0.25)/T;
    end
    p_val(i) = nnz(p_idx(2:end,i)>p_idx(1,i))/n_rep;
    
    % Values for the Adjacency matrix
    if p_val(i)<1e-3
        A(bat_pairs(i,1),bat_pairs(i,2)) = 0.7;
    elseif p_val(i)<1e-2
        A(bat_pairs(i,1),bat_pairs(i,2)) = 0.5;
    elseif p_val(i)<5e-2
        A(bat_pairs(i,1),bat_pairs(i,2)) = 0.3;
    else
        A(bat_pairs(i,1),bat_pairs(i,2)) = 0.01;
    end
end

%Create Graph
G = graph(A,cellstr(bat_nms),'upper');

%% Make a few figures

if options.show_fig
    % FIG: Scatter plot all bats
    figure();       set(gcf, 'units','normalized','outerposition',[0 0.25 1 0.5]);
    for i=1:n_tags
        subplot(131);  plot3(r(:,1,i),r(:,2,i),r(:,3,i),'Color', bat_clr(i,:));  xlim(r_lim(1,:)); ylim(r_lim(2,:)); zlim(r_lim(3,:));  title('3D view');                 hold on;
        subplot(132);  plot3(r(:,1,i),r(:,2,i),r(:,3,i),'Color', bat_clr(i,:));  xlim(r_lim(1,:)); ylim(r_lim(2,:)); zlim(r_lim(3,:));  title('Top view');   view(0,90);  hold on;
        subplot(133);  plot3(r(:,1,i),r(:,2,i),r(:,3,i),'Color', bat_clr(i,:));  xlim(r_lim(1,:)); ylim(r_lim(2,:)); zlim(r_lim(3,:));  title('Door view');  view(-30,0); hold on;
    end
    hold off;
    
    % FIG: Trajectories individual bats
    figure();   set(gcf, 'units','normalized','outerposition',[0 0.25 1 0.35]);
    for i=1:n_tags
        subplot(1,n_tags,i);  plot3(r(:,1,i),r(:,2,i),r(:,3,i),'-','Color', bat_clr(i,:));  xlim(r_lim(1,:)); ylim(r_lim(2,:)); zlim(r_lim(3,:));  title(bat_nms(i,:));  view(90,90);
        axis equal;
    end
    
    % FIG: Density Histograms 3D
    figure();   set(gcf, 'units','normalized','outerposition',[0 0.25 1 0.35]);
    for i=1:n_tags
        subplot(1,n_tags,i); hist3(r(:,1:2,i),'edges',edges_d,'CdataMode','auto','FaceColor',bat_clr(i,:));  xlim(r_lim(1,:)); ylim(r_lim(2,:));  title(bat_nms(i,:));
    end
    
    % FIG: Density Histograms heat-map
    figure();   set(gcf, 'units','normalized','outerposition',[0 0.25 1 0.35]);
    for i=1:n_tags
        sbpt = subplot(1,n_tags,i);
        hist3(r(:,1:2,i),'edges',edges_d,'CdataMode','auto','edgecolor','none','FaceColor','interp');
        xlabel('x');
        xlim(r_lim(1,:)); ylim(r_lim(2,:));   title(bat_nms(i,:));  view(90,90);  colormap(sbpt,custom_map(:,:,i)); % Change color scheme
        axis square;
    end
    
    % FIG: Velocity and flight segmentation
    figure();       set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    tiledlayout(n_tags,1,'TileSpacing','tight');
    for i = 1:n_tags
        ax(i) = nexttile;   %ax(i) = subplot(n_tags,1,i);
        area(t,bflying(:,i)*5,'FaceAlpha',0.3,'LineStyle','none');  hold on;
        area(t,wBeats(:,i)*-1,'FaceAlpha',0.3,'LineStyle','none');
        plot(t,v_abs(:,i),'.','Color', bat_clr(i,:));     plot(t,r(:,1,i),'k--');  ylabel('Velocity (m/s)');     hold off;
        legend('Fly','Wing-B','Vel','x(m)');
        title([num2str(f_num(i)) ' flights']);
    end
    linkaxes(ax,'x');   xlabel('Time(s)');
    
    % FIG: Few example flights and direction of motion
    for i = 1:n_tags
        figure(); set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
        nb = 1; sgtitle(bat_nms(i,:));
        for n = randperm(min(f_num),min([32;f_num]))
            chunk = f_smp{i,1}(n,1):f_smp{i,2}(n,1);
            subplot(4,8,nb);     plot(r(chunk,1,i),r(chunk,2,i),'k','LineWidth',3);    hold on;
            quiver(downsample(r(chunk,1,i),20),...
                downsample(r(chunk,2,i),20),...
                downsample(v_abs(chunk,i),20).*downsample(cos(angle(chunk,i)),20),...
                downsample(v_abs(chunk,i),20).*downsample(sin(angle(chunk,i)),20),1.5,'k');
            hold off;
            xlim(r_lim(1,:)); ylim(r_lim(2,:)); view(90,90);
            nb = nb+1;
        end
    end
    
    % FIG: Proximity index histograms
    figure();   set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    tiledlayout(round(n_tags/2),n_tags-1,'TileSpacing','tight');
    for i = 1:length(bat_pairs)
        ax(i) = nexttile;   %subplot(round(n_tags/2),n_tags-1,i);
        histogram(p_idx(:,i),'edgecolor','none');
        xlabel('Proximity Index'); title([bat_nms(bat_pairs(i,1),:) '-' bat_nms(bat_pairs(i,2),:)]);
        yticks([]);  xlim([0 1]);
        y1=get(gca,'ylim');  hold on; plot([p_idx(1,i) p_idx(1,i)],y1); hold off;
        title([bat_nms(bat_pairs(i,1),:) '-' bat_nms(bat_pairs(i,2),:) ': ' num2str(p_val(i))]);
    end
    
    % FIG: Network-graph
    figure();
    G.Edges.LWidths = 20*G.Edges.Weight; %7*G.Edges.Weight/max(G.Edges.Weight);
    p_NG = plot(G);
    p_NG.LineWidth = G.Edges.LWidths;  p_NG.MarkerSize = 10;  p_NG.NodeLabelColor = bat_clr;
    p_NG.NodeFontSize = 15;    p_NG.NodeFontWeight = 'bold';
    
    % FIG: Trajectories individual bats on top of histogram of the others
    figure();   set(gcf, 'units','normalized','outerposition',[0 0.25 1 0.35]);
    for i=1:n_tags
        colormap('copper') % Change color scheme
        subplot(1,n_tags,i);
        hist3(reshape(permute(r(:,1:2,find(1:n_tags ~= i)),[1 3 2]),[],2),'edges',edges_d,'CdataMode','auto','edgecolor','none','FaceColor','interp');
        xlabel('x');    view(0,90);   hold on;
        plot3(r(:,1,i),r(:,2,i),1e9.*r(:,3,i),'-','Color', 'w');  xlim(r_lim(1,:)); ylim(r_lim(2,:));   title(bat_nms(i,:));  view(90,90);
        hold off;
        axis square;
    end
    
    % FIG: Probabilities of flying together
    figure();   set(gcf, 'units','normalized','outerposition',[0 0 0.5 1]);
    scatter(p_fly2(:,2),p_fly2(:,1),'filled');
    xlabel('P coupled fly (independent)');  ylabel('P coupled fly (real)');
    text(p_fly2(:,2),p_fly2(:,1),bat_pair_nms,'VerticalAlignment','bottom','HorizontalAlignment','right')
    axis square;    xlim(ylim);     hold on;    h = refline(1,0);   hold off;   h.Color = 'k';  h.LineStyle = '--';
    
    % FIG: Percentages of time flying together
    figure();
    histogram(sum(bflying,2),'Normalization','probability');
    ylabel('Probability');  xlabel('Bats simult. flying');
    
    % FIG: Landing in ego-centric coordinates, trajectories
    figure();       set(gcf, 'units','normalized','outerposition',[0 0.25 1 0.5]);
    sphere_r = 0.5;
    for i=1:n_tags
        subplot(1,n_tags,i);
        title(bat_nms(i,:));
        for j=1:n_tags
            for n = 1:f_num(j)
                land_ck = f_smp{j,2}(n,1)-2*Fs:f_smp{j,2}(n,1);
                plot3(r_ego(land_ck,1,i,j),r_ego(land_ck,2,i,j),r_ego(land_ck,3,i,j),'-','Color', bat_clr(j,:));
                xlim([-sphere_r sphere_r]); ylim([-sphere_r sphere_r]); zlim([-sphere_r sphere_r]);   view(0,90);
                hold on;
                axis square;
            end
            plot3(0,0,0,'o','MarkerFaceColor', bat_clr(i,:),'MarkerEdgeColor','none','MarkerSize',10);
        end
        hold off;
    end
    
    % FIG: Landing in ego-centric coordinates, histograms
    figure();       set(gcf, 'units','normalized','outerposition',[0 0 0.5 1]);
    tiledlayout(n_tags,n_tags-1,'TileSpacing','compact','Padding','compact');
    sphere_r = 1;
    edges_s = {-sphere_r:(2*sphere_r)/20:sphere_r -sphere_r:(2*sphere_r)/20:sphere_r};                 %Edges for density histogram
    for i=1:n_tags
        for j=1:n_tags
            if j ~= i
                sbpt = nexttile;
                l_point = [];
                for n = 1:f_num(j)
                    l_point = cat(1,l_point,squeeze(r_ego(f_smp{j,2}(n,1),:,i,j)));
                end
                hist3(l_point(:,1:2),'edges',edges_s,'CdataMode','auto','edgecolor','none','FaceColor','interp');
                xlabel('x');
                xlim([-sphere_r sphere_r]); ylim([-sphere_r sphere_r]);    title(['Land on ' bat_nms(i,:)]);  view(0,90);
                axis square;
                colormap(sbpt,custom_map(:,:,j));
                hold on;
                line([0,0],[1,-1]); line([-1,1],[0,0]);
                hold off;
            end
        end
    end
    
    % FIG: Landing in ego-centric coordinates, grand-histogram
    l_point = [];
    for i=1:n_tags
        for j=1:n_tags
            if j ~= i
                for n = 1:f_num(j)
                    l_point = cat(1,l_point,squeeze(r_ego(f_smp{j,2}(n,1),:,i,j)));
                end
                
            end
        end
    end
    figure();       set(gcf, 'units','normalized','outerposition',[0 0.25 1 0.5]);
    tiledlayout(1,3,'TileSpacing','tight','Padding','compact');
    ax_nms = ['x';'y';'z'];
    for i=1:3
        sgtitle('Landing Spot Distributions');
        ax_idx = setdiff(1:3,i);
        nexttile;   hist3(l_point(:,ax_idx),'edges',edges_s,'CdataMode','auto','edgecolor','none','FaceColor','interp');
        xlabel(ax_nms(ax_idx(1)));    ylabel(ax_nms(ax_idx(2)));    xlim([-sphere_r sphere_r]); ylim([-sphere_r sphere_r]);    axis square;    hold on;    line([0,0],[1,-1],'Color','white'); line([-1,1],[0,0]); hold off;    view(0,90);
        colormap('gray');
    end
end

%% Save figures and data
if options.save_data
    figHandles = findall(0,'Type','figure');
    for i = 1:numel(figHandles)
        saveas(figHandles(i),[analysis_directory, '/', batdate '_figure' num2str(numel(figHandles)+1-i) '.png']);
    end
    close all;
    save([analysis_directory,'/Analyzed_DATA_', batdate, '.mat'],...
        'a','a_abs','a_flt','alpha_clus','angle','bat_clr','bat_dist','bat_nms','bat_pair_nms','bat_pairs',...
        'batdate','bflying','CDPmtdata','d_th','edges_d','edges_s','env','extracted_CDPfile',...
        'f_clus','f_land_table','f_num','f_smp','Fs','heading_diff','n_rep','n_tags','options',...
        'p_fly1','p_fly2','p_idx','p_val','r','r_ego','r_lim','r_old','r_qt','social_flag','stat_periods',...
        't','T','v','v_abs','v_th','wBeats');
end

%% Save movie at 10Hz with bat trajectories
if options.savemovie
    close all;
    r_ds = downsample(r,10);
    figure();       set(gcf, 'units','pixels','position',[300,300,280,260]);
    n=1;
    tic;    disp('...Saving Session Movie...');
    for j=1:length(r_ds)
        for i = 1:n_tags
            %plot(r_ds(j,1,i),r_ds(j,2,i),'.','MarkerFaceColor', bat_clr(i,:),'MarkerSize',60);  hold on;
            %mov_tx(i) = text(r_ds(j,1,i),r_ds(j,2,i), bat_nms(i,:));    hold on;
            mov_tx(i) = text(r_ds(j,2,i),-r_ds(j,1,i), bat_nms(i,:));    hold on; %rotated view, with feeders on the bottom
        end
        hold off;
        xlim(r_lim(1,:)); ylim(r_lim(2,:));   xticks([]); yticks([]);
        drawnow();
        cdp_movie(n) = getframe(gcf) ;
        cdp_movie(n).cdata=cdp_movie(n).cdata(:,:,1);   %convert to grayscale, if needed
        delete(mov_tx);
        n=n+1;
    end
    
    % create the video writer with 10 fps
    writerObj = VideoWriter([analysis_directory, '/', batdate '_10Hz_movie.avi'],'Grayscale AVI');
    writerObj.FrameRate = 100;  %10x speed
    open(writerObj);
    for i=1:length(cdp_movie)
        writeVideo(writerObj, cdp_movie(i));
    end
    close(writerObj);
    toc
    close all;
end