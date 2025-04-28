function [click_samples,click_times,click_amplitude] = EchoDetect_AF_v1(audio_data,fs)
%% Usage example

%[click_samples,click_times,click_amplitude] = EchoDetect_AF_v1(data,fs);  save('clicks.mat','click_times','click_amplitude');


%% Load data and initialize matrices

%fs = 192e3;
%chunk = [fs*3:fs*1200];
disp('Loading data...');
chunk = (1:length(audio_data));
raw_audio = audio_data(chunk);
n_samples = round(0.001*fs);
t = [0:1:size(raw_audio,1)-1]./fs;
th = 3;                                 %number of standard deviations for thresholding

%% Filter design, based on Yovel et al. "Optimal localization by pointing off axis."(2010), fig. S2

%===Define filter based on echolocation spectrum
fr_pts_yovel = [15 20 25 30 35 40 45 50 55 60 65 70 75 80]*1e3;
db_pts_yovel = -[40 13 2 0 5 12 14 18 20 26 25 26 33 43];
fr_pts = linspace(0,fs/2,10);
db_pts = interp1(fr_pts_yovel,db_pts_yovel,fr_pts,'linear','extrap');
%plot(fr_pts,db_pts,fr_pts_yovel,db_pts_yovel);

%===Design Yulewalk filter (or bandpass)
f = linspace(0,1,10);
m = db2mag(db_pts);
%[b,a] = yulewalk(8,f,m);
[b,a] = butter(10,[10e3*2/fs 80e3*2/fs],'bandpass');
[h,w] = freqz(b,a);

%===Visualize filter
if 0
    figure();
    plot(w/pi*fs/2,mag2db(abs(h)));     hold on;
    plot(fr_pts_yovel,db_pts_yovel);    hold off;
    yl = ylim;  xlabel('\omega/\pi');   ylabel('Magnitude');
end

%===Filter data
disp('Filtering data...');
flt_audio = filter(b,a,raw_audio,[],1);

%% Extract echolocation calls by thresholding

%===Clicks are emitted in pairs, typically separated by ~20 ms 
%===see Holland et al.,"Echolocation signal structure in the Megachiropteran bat Rousettus aegyptiacus Geoffroy 1810."2004
%===see also Waters et al., "Echolocation performance and Call Structure in the Magachiropteran Fruit-Bat Rousettus aegyptiacus" 2003

disp('Thresholding data...');
t_snip = [1:n_samples]./fs;

%===Generate envelope of the audio signal (better doing it in chunks, 10)
%=== If 10 does not work, consider increasing this number

%env_audio = zscore(envelope(flt_audio));
tic;
env_audio = [];
L = round(length(flt_audio)/10);
flt_audio_cut = flt_audio;
while ~isempty(flt_audio_cut)
    env_audio = [env_audio; zscore(envelope(flt_audio_cut(1:min(L,length(flt_audio_cut)))))];
    flt_audio_cut = flt_audio_cut(L+1:end);
end
toc;

%===Threshold (1 standard deviation) and get peaks with 50 ms minimal (to minimize double counting clicks)
%distance (work on the downsampled envelope)
[pot_peaks,pot_sampls] = findpeaks(downsample(env_audio,10),'MinPeakHeight',th,'MinPeakDistance',0.05*fs/10);
pot_sampls = pot_sampls*10;
pot_clicks = zeros(n_samples,length(pot_sampls));

%===Populate the matrix with suprathreshold events
for j = 1:length(pot_sampls)
    pot_clicks(:,j) = flt_audio(pot_sampls(j)+[1:n_samples]-round(n_samples/2));
end     

%% Look at the spectral domain 

disp('Calculating FFT...');
%===Improve speed of the FFT by getting the next power of 2 from the signal length
n = 2^nextpow2(n_samples);

%===Calculate fft and power spectrum
Y = fft(pot_clicks,n);
f = fs*(0:(n/2))/n;
P = abs(Y/n).^2;

%% Perform PCA and clustering

%===Perform PCA
[~,score,~,~,~] = pca(P');
alim = [min(score(:,1:3))'  max(score(:,1:3))'];


%===Manually select cluster, working on the first 3 PCAs
figure; set(gcf, 'units','normalized','outerposition',[0.2 0.2 0.55 0.33]);
tiledlayout(1,5,'TileSpacing','tight','Padding','tight');
nexttile;   scatter3(score(:,1),score(:,2),score(:,3),5,'filled'); xlabel('PCA1'); ylabel('PCA2'); zlabel('PCA3');
nexttile;   scatter(score(:,1),score(:,2),10,'filled');                       xlabel('PC1');  ylabel('PC2'); xlim('padded'); ylim('padded');  cluster_12 = drawpolygon;   tf_12 = inROI(cluster_12,score(:,1),score(:,2));   
nexttile;   scatter(score(:,1),score(:,3),10,tf_12*[0.1 0.5 0.2],'filled');   xlabel('PC1');  ylabel('PC3'); xlim('padded'); ylim('padded');  cluster_13 = drawpolygon;   tf_13 = inROI(cluster_13,score(:,1),score(:,3)) & tf_12;   
nexttile;   scatter(score(:,2),score(:,3),10,tf_13*[0.1 0.5 0.2],'filled');   xlabel('PC2');  ylabel('PC3'); xlim('padded'); ylim('padded');  cluster_23 = drawpolygon;   tf_23 = inROI(cluster_23,score(:,2),score(:,3)) & tf_12 & tf_13;
tf = tf_12 & tf_13 & tf_23;
nexttile;   scatter3(score(:,1),score(:,2),score(:,3),10,tf*[0.1 0.5 0.2],'filled');

%===Cluster putative clicks automatically via k-means
if 0
X = score(:,1:3);
idx = kmeans(X,3);
echo_clus = mode(idx);
col = [0.2 0.3 0.5].*ones(size(pot_clicks,2),1);      col(idx == echo_clus,:) = [0.1 0.5 0.2].*ones(size(col(idx == echo_clus,:)));
figure();
scatter3(score(:,1),score(:,2),score(:,3),10,col,'filled');
tf = idx == echo_clus;
end

%===Assign clicks to clustered data
clicks = pot_clicks(:,tf);
clk_lc = pot_sampls(tf);
clk_am = pot_peaks(tf);

%% Visualize features of the clicks

edges = 10.^(-4:0.1:0);
figure; set(gcf, 'units','normalized','outerposition',[0.2 0.1 0.6 0.5]);
tiledlayout(2,5,'TileSpacing','tight','Padding','compact');

nexttile;   plot(t_snip,pot_clicks,'k');   hold on;    line(round(n_samples/3)*[1 1]./fs,ylim);    hold off;
nexttile;   plot(t_snip,median(pot_clicks,2));
nexttile;   histogram(diff(pot_sampls./fs),edges);     set(gca,'XScale','log');
nexttile;   histogram(pot_peaks);
nexttile;   plot(f,P(1:n/2+1,:),'k');  xlabel('F (Hz)'); ylabel('Power');

nexttile;   plot(t_snip,clicks,'k');   hold on;    line(round(n_samples/3)*[1 1]./fs,ylim);    hold off;
nexttile;   plot(t_snip,median(clicks,2));
nexttile;   histogram(diff(clk_lc./fs),edges);     set(gca,'XScale','log');
nexttile;   histogram(pot_peaks(tf));   xlabel('Amplitude');    ylabel('Counts');
nexttile;   
plot(f,P(1:n/2+1,tf),'k');  hold on;    plot(f,P(1:n/2+1,~tf),'r'); hold off;   xlabel('F (Hz)'); ylabel('Power');
axes('Position',[.85 .2 .1 .2]); box on;
plot(f,mag2db(mean(P(1:n/2+1,tf),2)),'k');  %hold on;    plot(f,mean(P(1:n/2+1,~tf),2),'r'); hold off;

%% Plot a few examples
figure; set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
tiledlayout(5,10,'TileSpacing','none','Padding','none');

for i = randi(size(clicks,2),[1 5*5])
    y = flt_audio(clk_lc(i)+round([-fs*0.01:fs*0.01]));
    nexttile;   plot([1:length(y)]./fs,y);  xlim([0 0.02]);
    nexttile;   spectrogram(y,50,25,[],fs,'yaxis');
    colormap turbo;
    colorbar off;
end

%% Assign output
click_samples = clk_lc; 
click_times = clk_lc./fs;
click_amplitude = clk_am;

end

