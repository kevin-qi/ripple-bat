function [click_time,click_shape,click_power,fs] = EchoDetect_AF_v3(filename)
%% Brief script to detect echolocation clicks form audio recordings in the flight room
% AF, September 2022

%=== INPUTS
% data:     Audio trace
% fs:       Sampling frequency

%=== OUTPUTS
% click_time:   Times of the clicks,
% click_shape:  Matrix with the waveform for all the detected clicks (1 ms)
% click_power:  Matrix with the power spectrum for all the detected clicks
% f_FFT:        Frequencies for the power spectrum
% fs:           Sampling frequency

%=== USAGE EXAMPLE
% [click_time,click_shape,click_power,fs] = EchoDetect_AF_v2(filename);

%% LOAD DATA AND PARAMETERS

%=== Load data
[data,fs] = audioread(filename);
debug = 0;
disp('Loading data...');
if debug,chunk = [fs*1:fs*3000];
else, chunk = (1:length(data));   end
raw_audio = data(chunk);

%=== Definitions and parameters
file_name_parts = split(filename,'_');
batdate = file_name_parts{3,1};                                             % Get identifier of the session
mic_ref = file_name_parts{2,1};                                             % Get Mic identifier
threshold = 5;                                                              % Threshold (sd) for detection of events (default: 10)
kmeansclus = 0;                                                             % If using kmeans or hierarchical
filter_freq = [10e3 min(40e3,fs/2-1)];                                      % Bandpass Filter frequencies
[b,a] = butter(8,[filter_freq(1)*2/fs filter_freq(2)*2/fs],'bandpass');     % Filter definition
min_time = 10e-3;                                                           % Minimum time between two events
iei_edges = 10.^[-2:0.02:0];                                                % Edges for the inter-event-interval histogram
n_samples = round(0.001*fs);                                                % Number of samples (1 ms) to look at the waveform
time_snippet = [0:n_samples-1]/fs-n_samples/(2*fs);                         % Time around the click
n_optimal = 2^nextpow2(n_samples);                                          % Optimal number of samples for FFT
f_FFT = fs*(0:(n_optimal/2))/n_optimal;                                     % FFT Frequencies
manual_cluster = 0;                                                         % If using automatic clustering

%% PROCESS AUDIO DATA

%=== Initial processing
disp('Filtering data...');      flt_audio = filter(b,a,raw_audio,[],1);
disp('Z-scoring data...');      zscore_audio = zscore(flt_audio);
disp('Finding peaks... ');      [pot_ampl,pot_samples] = findpeaks(zscore_audio,'MinPeakHeight',threshold,'MinPeakDistance',fs*min_time);

%=== Visualize inter-event-intervals
figure; set(gcf, 'units','normalized','outerposition',[0.2 0.2 0.3 0.4]);
nexttile;   histogram(diff(pot_samples)/fs,iei_edges,'edgecolor','none');   set(gca, 'XScale', 'log');
hold on;    plot(0.02*[1 1],ylim,'r--');    plot(0.1*[1 1],ylim,'r--'); hold off;   xlabel('Inter event time (s)');

%=== Populate the matrix of potential clicks with suprathreshold events
pot_clicks = zeros(n_samples,numel(pot_samples));
for j = 1:length(pot_samples)
    pot_clicks(:,j) = flt_audio(pot_samples(j)+[1:n_samples]-round(n_samples/2));
end

%===Calculate fft and power spectrum, apply PCA
disp('Calculating FFT...');     Y = fft(pot_clicks,n_optimal); P = abs(Y/n_optimal).^2;    P = P(1:n_optimal/2+1,:);
disp('Applying PCA...');        [~,score] = pca(P');

%% DENOISE OR CLUSTER THE CLICKS

if manual_cluster
    
    %=== Manually select cluster, working on the first 3 PCAs
    figure('units','normalized','outerposition',[0 0.2 1 0.4]);
    tiledlayout(1,5,'TileSpacing','tight','Padding','tight');
    nexttile;   scatter3(score(:,1),score(:,2),score(:,3),5,'filled'); xlabel('PCA1'); ylabel('PCA2'); zlabel('PCA3');
    nexttile;   scatter(score(:,1),score(:,2),5,'filled');                       xlabel('PC1');  ylabel('PC2');   cluster_12 = drawpolygon;   tf_12 = inROI(cluster_12,score(:,1),score(:,2));
    nexttile;   scatter(score(:,1),score(:,3),5,tf_12*[0.1 0.5 0.2],'filled');   xlabel('PC1');  ylabel('PC3');   cluster_13 = drawpolygon;   tf_13 = inROI(cluster_13,score(:,1),score(:,3)) & tf_12;
    nexttile;   scatter(score(:,2),score(:,3),5,tf_13*[0.1 0.5 0.2],'filled');   xlabel('PC2');  ylabel('PC3');   cluster_23 = drawpolygon;   tf_23 = inROI(cluster_23,score(:,2),score(:,3)) & tf_12 & tf_13;
    tf = tf_23;
    nexttile;   scatter3(score(:,1),score(:,2),score(:,3),10,tf*[0.1 0.5 0.2],'filled');
    
else
    
    if kmeansclus
        
        %=== Automatically cluster with k means and keep heaviest cluster
        idx = kmeans(score(:,1:3),4);
        figure('units','normalized','outerposition',[0 0.2 0.3 0.5]);
        gscatter(score(:,1),score(:,2),idx);
        [~,heavy_clusters] = sort(groupcounts(idx),'descend');
        tf = (idx == heavy_clusters(1));
        title(['Cluster ',num2str(heavy_clusters(1)),' likely composed of clicks']);
        
    else
        
        %=== Use hierarchical clustering
        distance = 10e-5;
        X = score(:,1:3);                                                                            
        Y = pdist(X,'euclidean'); 
        Z = linkage(Y,'average');
        %Z = linkage(Y,'single');
        id_temp = cluster(Z,'Cutoff',distance,'Criterion','distance');  % perform clustering
        [Ns,~,b] = histcounts(id_temp,[1:max(id_temp)+1]-0.5);          % calculate number of points/cluster
        id_temp(Ns(b)<10) = 0;                                          % Force clusters with less than 10 points to have id 0
        [unique_labels, ~, idx_tmp] = unique(id_temp);                  % Find unique cluster labels 
        counts = accumarray(idx_tmp, 1);                                % Count points in them
        [~, sort_order] = sort(counts, 'descend');                      % Sort counts in descending order and get the sort order
        new_labels = zeros(size(unique_labels));                        % Create a mapping from original labels to new labels
        new_labels(sort_order) = 1:length(unique_labels);
        idx= new_labels(idx_tmp);                                       % Reassign cluster identities using the new labels
        
        figure('units','normalized','outerposition',[0 0.2 0.3 0.5]);
        gscatter(score(:,1),score(:,2),idx);
        tf = (idx == 1);
        title('Cluster 1 likely composed of clicks');
        
        % Add circles for optimal clustering eval
        hold on;
        center = mean(score(tf,1:2));   theta = linspace(0, 2*pi, 100);
        x1 = center(1) + 0.5*10e-5 * cos(theta);     y1 = center(2) + 0.5*10e-5 * sin(theta); 
        x2 = center(1) + 1.0*10e-5 * cos(theta);     y2 = center(2) + 1.0*10e-5 * sin(theta);
        x3 = center(1) + 3.0*10e-5 * cos(theta);     y3 = center(2) + 3.0*10e-5 * sin(theta);
        plot(x1,y1,'b-',x2,y2,'r-',x3,y3,'g-'); % plot the circle in blue color
        axis equal;
        
    end
    
end

%===Assign clicks to clustered data
click_shape = pot_clicks(:,tf);
click_sampl = pot_samples(tf);
click_power = P(:,tf);
click_time = click_sampl/fs;

%% Visualize features of the clicks

figure('units','normalized','outerposition',[0.2 0.1 0.7 0.3]);
tiledlayout(1,5,'TileSpacing','tight','Padding','compact');
nexttile;   plot(time_snippet,click_shape,'k');     title('All detected events');   xlabel('Time (s)');
nexttile;   plot(time_snippet,mean(click_shape,2)); title('Average waveform');      xlabel('Time (s)');
nexttile;   histogram(max(click_shape,[],1),'edgecolor','none');    title('Amplitude distribution');    xlabel('Mic level');
nexttile;   histogram(diff(click_time),iei_edges,'edgecolor','none');   hold on;    plot(0.02*[1 1],ylim,0.1*[1 1],ylim,'r--'); hold off;
xlabel('Inter-event-time (s)'); set(gca,'XScale','log');    title('Inter-event time');
nexttile;   plot(f_FFT,mag2db(click_power),'k');  xlabel('F (Hz)'); ylabel('dB');    title('Power spectrum');
hold on;    plot(f_FFT,mag2db(mean(click_power,2)),'r','LineWidth',2);  xlabel('F (Hz)'); ylabel('dB');

%% Save figures and data
figHandles = findall(0,'Type','figure');
for i = 1:numel(figHandles)
    saveas(figHandles(i),['Detected_Clicks', num2str(numel(figHandles)+1-i),'_',mic_ref, '.png']);
    saveas(figHandles(i),['Detected_Clicks', num2str(numel(figHandles)+1-i),'_',mic_ref, '.fig']);
end
close all;

%=== Save structure
Detected_Clicks.times = click_time;
Detected_Clicks.shape = click_shape';
Detected_Clicks.amp = max(click_shape)';
Detected_Clicks.power = click_power';
Detected_Clicks.fs = fs;
save(['Detected_Clicks_',batdate,'_',mic_ref,'.mat'],'Detected_Clicks');

%% Plot a few examples
% figure; set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
% tiledlayout(5,10,'TileSpacing','none','Padding','none');
% 
% for i = randi(numel(click_sampl),[1 5*5])
%     interval = round([-fs*0.05:fs*0.05]);
%     y = flt_audio(click_sampl(i)+interval);
%     nexttile;   plot(interval/fs,y);
%     nexttile;   spectrogram(y,50,25,[],fs,'yaxis');
%     colormap turbo;
%     colorbar off;
% end

end

