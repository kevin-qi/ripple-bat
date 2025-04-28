function audio_concat_AF_v1()
%% Brief script to concatenate audio data recorded in the flight room with 'record_audio_tobias'
% (modified from T.A.S, 2020)
% Last update by AF on 220914

%Recorded audio chunks are 12s long, saved every 10s (or 12s) and with a
%superposition of about 2s between the end and the start of adjacent chunks
%Each chunk is sampled at 192 kHz and contains n channels (usually 5 mics
%and 1 synchronization channel recording a 50ms TTL pulse happening every 3s)
% !!!------------BE CAREFUL that the mic to channel mapping may change between datasets,
% double check that the TTL is in the last recorded channel and adjust the
% number of mics accordingly ----------- !!!
% Mics are ordered 3,1,2,4 from the left to the right of the feeders. 
% Mics 1 and 2 are high quality, mics 3 and 4 are those on the walls.

%===Get Files, sorted by creation date
filePath = cd;
fileList = dir([filePath '\*audio_trial_*']);
name_parts = strsplit(fileList(1).name,'_');
batdate = char(name_parts(1,2));
[~,idx] = sort(datenum({fileList.date},'dd-mm-yyyy HH:MM:SS'));
fileList = fileList(idx);
n_chunks = size(fileList,1);        % number of chunks
n_mics = 5;                         % number of microphones
mics_2_save = [1:4];                % ids of the microphones to save
n_ds = 2;                           % downsampling factor
check_c3d = 0;                      % if checking c3d TTLs

%===Check when the chunks were created
D_cum = [0 0 0 0 0 0];          %[y m d h m s];
for i = 1:size(fileList,1)-1
    D_cum = [D_cum; datevec(fileList(i+1).datenum-fileList(i).datenum)];
end
%D_cum = D_cum(2:end,:); disp('Chunks saved every (s):') disp(unique(D_cum(:,6)));
disp('Concatenating Audio Files...');

%% ALIGNMENT 

% Alignment is done by cutting repeated segments
% It is crucial - for speed purposes - to allocate data in a cell array
warning('off','all');

%===Initialize the cell containing the audio and 
%===detect the first chunk containing a TTL pulse
audio_cell = cell(n_chunks,1);
i=1; no_TTL = 1;    chunk_dur = [];
while no_TTL
    load(fileList(i).name);
    chunk_dur = [chunk_dur; size(recbuf,1)/fs];
    signal = downsample(recbuf,n_ds);
    if max(signal(:,end))>0.2
        for j=1:n_mics
        audio_cell{i,j} = signal(:,j);
        end
        no_TTL = 0;
    end
    i = i+1;
end
first_chunk = i-1;  

%===Align the chunks by cutting repeated data
for i = first_chunk:n_chunks-1
    
    %=== Display processing
    if ~mod(i,50)
    disp(['Alignment of chunks after ', num2str(i,3)]); 
    end
    chunk_dur = [chunk_dur; size(recbuf,1)/fs];
    
    %--Load 1st chunk
    load(fileList(i).name,'recbuf');
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
    load(fileList(i+1).name,'recbuf');
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
        for j=1:n_mics
            audio_cell{i+1,j} = next_chunk(:,j);
        end
    end
end
warning('on','all');

%=== Delete empty cells
audio_cell = audio_cell(~any(cellfun('isempty',audio_cell),2),:);

%=== Check distance beween TTLs
disp('Detecting TTLs...');
audio_TTL = cell2mat(cellfun(@(x) downsample(x,1),audio_cell(:,end),'UniformOutput',false));
[~,~,UT,~,~] = risetime(audio_TTL,round((fs/n_ds)/1));

%=== Display some diagnostic
disp('---------------------------------------------');
disp(['Unique Chunk Durations (s):',num2str(unique(chunk_dur))]);
disp('---------------------------------------------');
figure; set(gcf, 'units','normalized','outerposition',[0.2 0.3 .4 0.3]);
tiledlayout(1,3,'TileSpacing','tight','Padding','compact');
nexttile;   plot(signal(:,end));    hold on;    plot(xlim,[.2 .2],'r--');   hold off;   title('First chunk with TTL');
nexttile;   plot(diff(UT),'.'); ylabel('Time difference between TTLs (s)')
nexttile;   histogram(diff(UT));ylabel('Counts'); xlabel('Time difference between TTLs (s)');
sgtitle([num2str(numel(UT)),' Detected TTLs']);

%% CUT AUDIO FILES BETWEEN FIRST AND LAST TTLs AND SAVE

%=== Save Data
if ~exist('Concatenated_Audio','dir');mkdir('Concatenated_Audio');end
for m = 1:n_mics  
    if any(m == mics_2_save)
        warning('off');
        conc_signal = cell2mat(audio_cell(:,m));                                    %===Concatenate chunks from a given microphone
        conc_signal = conc_signal(round(UT(1)*fs/n_ds):round(UT(end)*fs/n_ds),:);   %===Cut file between first and last TTLs
        disp(['Saving data from Mic ' num2str(m)]);                                 %===Write audio file
        audiowrite(['Concatenated_Audio/Flight_Mic' num2str(m) '_' batdate '_.wav'],conc_signal,round(fs/n_ds));
        warning('on');
    end
end

%=== Save figures
figHandles = findall(0,'Type','figure');
for i = 1:numel(figHandles)
    saveas(figHandles(i),['Concatenated_Audio/', batdate, '_micfigure', num2str(numel(figHandles)+1-i), '.png']);
    saveas(figHandles(i),['Concatenated_Audio/', batdate, '_micfigure', num2str(numel(figHandles)+1-i), '.fig']);
end
close all;

%% Check if the same number of TTLs were detected by cortex and the mic system

if check_c3d
    grandparent = fileparts(fileparts(fileparts(fileList(1).folder)));
    for i=1:2
        if i==1,c3d_clus_file = dir(fullfile(grandparent,'tracking','**','*1-Bat_Cluster.mat'));
        else,   c3d_clus_file = dir(fullfile(grandparent,'tracking','**','*2-Bat_Cluster.mat'));end
        if ~isempty(c3d_clus_file)
            clus_data = load([c3d_clus_file.folder '/' c3d_clus_file.name]);
            c3d_t = [0: 1/clus_data.AnalogFrameRate : (length(clus_data.AnalogSignals)-1)/clus_data.AnalogFrameRate]';
            c3d_TTL.evn_signal = [0; normalize(diff(clus_data.AnalogSignals(:,2)),'range',[-1 1])];
            [~, c3d_TTL.evn_times] = findpeaks(c3d_TTL.evn_signal,c3d_t,'MinPeakHeight',0.5,'MinPeakDistance',2);
            if length(c3d_TTL.evn_times)~=numel(UT)
                disp(['Cortex File ',num2str(i),': Something is wrong with TTLs']);
            else
                disp(['Cortex File ',num2str(i),': All TTLs detected during fly session!!']);
            end
        end
    end
end

end