function audio_concat_AF_v0()

% Brief script to concatenate audio data recorded in the flight room with 'record_audio_tobias'
% (modified from T.A.S, 2020)

%Recorded audio chunks are 12s long, saved every 10s (or 12s) and with a
%superposition of about 2s between the end and the start of adjacent chunks
%Each chunk is sampled at 192 kHz and contains n channels (usually 5 mics
%and 1 synchronization channel recording a 50ms TTL pulse happening every 3s)
% Mics are ordered 3,1,2,4 from the left to the right of the feeders

%===Get Files, sorted by creation date
filePath = cd;
fileList = dir([filePath '\*audio_trial_*']);
name_parts = strsplit(fileList(1).name,'_');
batdate = char(name_parts(1,2));
[~,idx] = sort(datenum({fileList.date},'dd-mm-yyyy HH:MM:SS'));
fileList = fileList(idx);
fileFirst = fileList(1);
n_chunks = size(fileList,1);
n_mics = 6;
n_ds = 1; 
debug = 0;
check_c3d = 0;

%===Check when the chunks were created
D_cum = [0 0 0 0 0 0];          %[y m d h m s];
for i = 1:size(fileList,1)-1
    D_cum = [D_cum; datevec(fileList(i+1).datenum-fileList(i).datenum)];
end
D_cum = D_cum(2:end,:);
disp('Chunks saved every (s):')
disp(unique(D_cum(:,6)));

%% ALIGNMENT 
% Alignment is done by cutting repeated segments
% It is crucial - for speed purposes - to allocate data in a cell array
warning('off','all');

%===Initialize the cell containing the audio and 
%===detect the first chunk containing a TTL pulse
audio_cell = cell(n_chunks,1);
i=1; no_TTL = 1;
while no_TTL
    load(fileList(i).name);
    signal = downsample(recbuf,n_ds);
    if max(signal(:,end))>0.2
        audio_cell{i} = signal;
        no_TTL = 0;
    end
    i = i+1;
end

first_chunk = i-1;  

%===Align the chunks by cutting repeated data
for i = first_chunk:n_chunks-1
    disp(['Alignment of chunk ', num2str(i,3)]);    
    %--Load 1st chunk
    load(fileList(i).name);
    signal_1 = downsample(recbuf,n_ds);
    %--Check if TTL is present
    if max(signal_1(:,end))>0.2
        %--Get TTL samples
        [~,~,UT,~,~] = risetime(signal_1(:,end));
        %--Get samples between the last TTL and the end of the chunk
        smpl2end = 12*round(fs/n_ds)-UT(end);
    else
        smpl2end = nan;
    end
    
    %--Load 2nd chunk
    load(fileList(i+1).name);
    signal_2 = downsample(recbuf,n_ds);
    %--Check if TTL is present
    if max(signal_2(:,end))>0.2
        %--Get TTL samples
        [~,~,UT,~,~] = risetime(signal_2(:,end));
        %--Get samples between the start of the chunk and the 1st TTL
        smpl2str = UT(1);
    else
        smpl2str = nan;
    end
    
    %--Proceed only with files containing TTLs
    if ~isnan(smpl2end) && ~isnan(smpl2str)
        if smpl2end+smpl2str<3*round(fs/n_ds)
            next_chunk = signal_2(smpl2str+smpl2end:end,:);
        else
            next_chunk = signal_2(smpl2end+smpl2str-3*round(fs/n_ds):end,:);
        end
        audio_cell{i+1} = next_chunk;
    end
end
warning('on','all');

%===Delete empty cells
audio_cell = audio_cell(~cellfun('isempty',audio_cell));

%===Check alignment quality
audio_TTL = [];
for i = 1:size(audio_cell,1)
    audio_TTL = [audio_TTL; downsample(audio_cell{i}(:,end),100)];
end
[~,~,UT,~,~] = risetime(audio_TTL,round(fs/100));
figure; set(gcf, 'units','normalized','outerposition',[0.2 0.3 .2 0.6]);
tiledlayout(2,1,'TileSpacing','tight','Padding','compact');
nexttile;   plot(diff(UT),'.'); ylabel('Time difference between TTLs (s)')
nexttile;   histogram(diff(UT));ylabel('Counts'); xlabel('Time difference between TTLs (s)');
title(num2str(unique(round(diff(UT)*1000))));

%% CUT AUDIO FILES BETWEEN FIRST AND LAST TTLs AND SAVE

%=== Save Data
for m = 1:n_mics  
    %===Concatenate chunks from a given microphone
    conc_signal = [];
    for i = 1:size(audio_cell,1)
        conc_signal = [conc_signal; downsample(audio_cell{i}(:,m),1)];
    end
    
    %===Cut file between first and last TTLs
    conc_signal = conc_signal(UT(1)*round(fs/n_ds):UT(end)*round(fs/n_ds),:);
    
    %===Write audio file
    disp(['Saving data from Mic ' num2str(m)]);
    audiowrite(['Flight_Mic' num2str(m) '_' batdate '_.wav'],conc_signal,round(fs/n_ds));
end

%=== Save figures
figHandles = findall(0,'Type','figure');
for i = 1:numel(figHandles)
    saveas(figHandles(i),[batdate, '_micfigure', num2str(numel(figHandles)+1-i), '.png']);
    saveas(figHandles(i),[batdate, '_micfigure', num2str(numel(figHandles)+1-i), '.fig']);
end
close all;

%% Check if the same number of TTLs were detected by cortex and the mic system

if check_c3d
    c3d_clus_file = dir(fullfile(cd, '*_c3d_*-Bat_*.mat'));
    clus_data = load(c3d_clus_file.name);
    c3d_t = [0: 1/clus_data.AnalogFrameRate : (length(clus_data.AnalogSignals)-1)/clus_data.AnalogFrameRate]';
    c3d_TTL.evn_signal = [0; normalize(diff(clus_data.AnalogSignals(:,2)),'range',[-1 1])];
    [~, c3d_TTL.evn_times] = findpeaks(c3d_TTL.evn_signal,c3d_t,'MinPeakHeight',0.5,'MinPeakDistance',2);
    r = risetime(data);
    if length(c3d_TTL.evn_times)~=length(r)
        error('Something is wrong with TTLs');
    else
        disp('All TTLs detected during fly session!!');
    end
end



%% ALIGNMENT (Debug mode)

if debug
    conc_signal = cell(n_chunks,1);
    tic;
    %warning('off','all');
    n_ds = 10;
    %===Detect the first chunk containing a TTL pulse
    i=1; no_TTL = 1;
    while no_TTL
        load(fileList(i).name);
        signal = downsample(recbuf(:,end),n_ds);
        if max(signal(:,end))>0.2
            conc_signal = signal;
            %conc_signal{i} = signal;
            no_TTL = 0;
        end
        i = i+1;
    end
    first_chunk = i-1;
    %===Align the chunks by cutting repeated data
    for i = first_chunk:n_chunks-1
        disp(i);
        %--Load 1st chunk
        load(fileList(i).name);
        signal_1 = downsample(recbuf(:,end),n_ds);
        %--Check if TTL is present
        if max(signal_1(:,end))>0.2
            %--Get TTL samples
            [~,~,UT,~,~] = risetime(signal_1(:,end));
            %--Get samples between the last TTL and the end of the chunk
            smpl2end = 12*round(fs/n_ds)-UT(end);
        else
            smpl2end = nan;
        end
        %--Load 2nd chunk
        load(fileList(i+1).name);
        signal_2 = downsample(recbuf(:,end),n_ds);
        %--Check if TTL is present
        if max(signal_2(:,end))>0.2
            %--Get TTL samples
            [~,~,UT,~,~] = risetime(signal_2(:,end));
            %--Get samples between the start of the chunk and the 1st TTL
            smpl2str = UT(1);
        else
            smpl2str = nan;
        end
        %--Proceed only with files containing TTLs
        if ~isnan(smpl2end) && ~isnan(smpl2str)
            if smpl2end+smpl2str<3*round(fs/n_ds)
                start = uint32(smpl2str+smpl2end);
                next_chunk = signal_2(start:end,:);
            else
                start = uint32(smpl2end+smpl2str-3*round(fs/n_ds));
                next_chunk = signal_2(start:end,:);
            end
            conc_signal = [conc_signal; next_chunk];
            %conc_signal{i} = next_chunk;
        end
    end
    toc;
    %warning('on','all');
end

end