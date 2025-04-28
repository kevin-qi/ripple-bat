% Script for the detection of echolocation calls

%===Parameters and data
n_mics = 1;
fs = 192000;
n_samples = round(0.002*fs);
%raw_audio = audio_cell{6, 1}(:,1:n_mics);
raw_audio = data;
t = [0:1:size(raw_audio,1)-1]./fs;

%% Filter design, based on Yovel et al. "Optimal localization by pointing off axis."(2010), fig. S2

%===Define filter based on echolocation spectrum
fr_pts_yovel = [15 20 25 30 35 40 45 50 55 60 65 70 75 80]*1e3;
db_pts_yovel = -[40 13 2 0 5 12 14 18 20 26 25 26 33 43];
fr_pts = linspace(0,fs/2,10);
db_pts = interp1(fr_pts_yovel,db_pts_yovel,fr_pts,'linear','extrap');
plot(fr_pts,db_pts,fr_pts_yovel,db_pts_yovel);

%===Design Yulewalk filter (or bandpass)
f = linspace(0,1,10);
m = db2mag(db_pts);
[b,a] = yulewalk(8,f,m);
%[b,a] = butter(10,[20e3*2/fs 60e3*2/fs],'bandpass');
[h,w] = freqz(b,a);
plot(w/pi,mag2db(abs(h)));

%===Visualize filter
plot(w/pi*fs/2,mag2db(abs(h)));     hold on;
plot(fr_pts_yovel,db_pts_yovel);    hold off;
yl = ylim;  xlabel('\omega/\pi');   ylabel('Magnitude');

%% Filter data 
tic
flt_audio = filter(b,a,raw_audio,[],1);
toc
% figure; set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
% tiledlayout(n_mics,1,'TileSpacing','none','Padding','compact');
% for i = 1:n_mics
% ax(i) = nexttile; plot(t,raw_audio(:,i),t,flt_audio(:,i)-0.5);
% end
% linkaxes(ax,'x');
% xlim([t(1) t(end)]);

%% Extract echolocation calls

%===Clicks are emitted in pairs, typically separated by ~18 ms 
%===see Holland et al.,"Echolocation signal structure in the Megachiropteran bat Rousettus aegyptiacus Geoffroy 1810."2004

clear loc pot_sampls pot_clicks clicks;

%===Threshold on the envelope!!!

%===Threshold (3 times standard deviation) and get peaks with 50 ms minimal distance
tic
trace = flt_audio(:,1);
th = 3*std(trace,[],1);
t_snip = [1:n_samples]./fs;
[~,pot_sampls] = findpeaks(trace,'MinPeakHeight',th,'MinPeakDistance',0.05*fs);
pot_clicks = zeros(n_samples,length(pot_sampls));
toc

%===Populate the matrix with suprathreshold events
for j = 1:length(pot_sampls)
    %pot_clicks(:,j) = normalize(trace(pot_sampls(j)+[1:n_samples]-round(n_samples/3)),'range',[-1 1]);
    pot_clicks(:,j) = trace(pot_sampls(j)+[1:n_samples]-round(n_samples/3));
end

%% Visualize some statistics

edges = 10.^(-4:0.1:0);
figure; set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
tiledlayout(1,3,'TileSpacing','tight','Padding','compact');
nexttile;   plot(t_snip,pot_clicks,'k');   hold on;    line(round(n_samples/3)*[1 1]./fs,ylim);    hold off;
nexttile;   plot(t_snip,median(pot_clicks,2));
nexttile;   histogram(diff(pot_sampls./fs),edges);     set(gca,'XScale','log');

%% Perform PCA and clustering on the extracted events
figure; set(gcf, 'units','normalized','outerposition',[0.3 0.3 0.1 0.3]);
[coeff,score,latent,tsquared,explained] = pca(pot_clicks');
scatter3(score(:,1),score(:,2),score(:,3),'filled');
xlabel('PCA1'); ylabel('PCA2'); zlabel('PCA3');

%===Cluster putative clicks (hierarchical or k-means)
X = score(:,1:3);
%Y = pdist(X,'euclidean');
%Z = linkage(Y,'single');
%dendrogram(Z,0);
%idx = cluster(Z,'Cutoff',0.03,'Criterion','distance');
idx = kmeans(X,2);
echo_clus = mode(idx);
col = [0.2 0.3 0.5].*ones(size(pot_clicks,2),1);      col(idx == echo_clus,:) = [0.1 0.5 0.2].*ones(size(col(idx == echo_clus,:)));
scatter3(score(:,1),score(:,2),score(:,3),10,col,'filled');

%===Keep the main cluster
clicks = pot_clicks(:,idx == echo_clus);
clk_lc = pot_sampls(idx == echo_clus);

%===Manually select cluster
if 1
figure; set(gcf, 'units','normalized','outerposition',[0.2 0.2 0.4 0.25]);
tiledlayout(1,4,'TileSpacing','tight','Padding','tight');
nexttile;   scatter(score(:,1),score(:,2),10,'filled');                       xlabel('PC1');  ylabel('PC2');  cluster_12 = drawpolygon;   tf_12 = inROI(cluster_12,score(:,1),score(:,2));   
nexttile;   scatter(score(:,1),score(:,3),10,tf_12*[0.1 0.5 0.2],'filled');   xlabel('PC1');  ylabel('PC3');  cluster_13 = drawpolygon;   tf_13 = inROI(cluster_13,score(:,1),score(:,3));   
nexttile;   scatter(score(:,2),score(:,3),10,tf_13*[0.1 0.5 0.2],'filled');   xlabel('PC2');  ylabel('PC3');  cluster_23 = drawpolygon;   tf_23 = inROI(cluster_23,score(:,2),score(:,3));
tf = tf_12 & tf_13 & tf_23;
nexttile;   scatter3(score(:,1),score(:,2),score(:,3),10,tf*[0.1 0.5 0.2],'filled')
clicks = pot_clicks(:,tf);
clk_lc = pot_sampls(tf);
end

%===Visualize some statistics
edges = 10.^(-4:0.1:0);
figure; set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
tiledlayout(1,3,'TileSpacing','tight','Padding','compact');
nexttile;   plot(t_snip,clicks,'k');   hold on;    line(round(n_samples/3)*[1 1]./fs,ylim);    hold off;
nexttile;   plot(t_snip,median(clicks,2));

click_times = clk_lc./fs; 

%% Visualize and hear trace and detected clicks

plot(t,normalize(raw_audio(:,i),'range'));  hold on;    stem(clk_lc./fs,ones(size(clk_lc)));    hold off;

audio_clicks = [];
for j = 1:size(clicks,2)
    audio_clicks = [audio_clicks; raw_audio(clk_lc(j)+[1:0.01*fs]-3840,1)];
end

audiowrite('clicks.wav',audio_clicks,fs);

%===Look at cross correlation
[c,lags] = xcorr(audio_clicks,0.5*fs,'normalized');
plot(lags/fs,c,'-');