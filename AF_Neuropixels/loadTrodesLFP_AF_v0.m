function [out] = loadTrodesLFP_AF_v0
%% Modified from loadTrodesLFP, originally from KQ. Feb 2024
% To be used in the .rec folder

%=== Get file names and paths
path_to_recording_dir = cd;
[~, dirname, ext] = fileparts(path_to_recording_dir);
assert(strcmp(ext, ".rec"), "Trodes recording directory must end in .rec");
mergedLFP_dirname = dirname + "_merged.LFP";
mergedTime_filename = dirname + "_merged.timestamps.dat";

%% Channel Map
xPos = repmat([8 -24 24 -8].', [240 1]);                    % Horizontal pitch is 16 um
yPos = reshape( repmat( [0:960/2], 2,1 ), 1, [] )*20+100;   % Vertical pitch is 20 um

%% Load LFP sample timestamps
timestamps = readTrodesExtractedDataFile(fullfile(path_to_recording_dir, ...
    mergedLFP_dirname, ...
    mergedTime_filename));
local_sample_timestamps_usec = 1e6 * double(timestamps.fields.data) / double(timestamps.clockrate);
timestamp_at_creation_usec = 1e6 * double(timestamps.timestamp_at_creation) / double(timestamps.clockrate);
first_timestamp_usec = 1e6 * double(timestamps.first_timestamp) / double(timestamps.clockrate);

dio = loadTrodesDigital(path_to_recording_dir);
ttl_timestamps_usec = dio{1}.ttl_timestamp_usec;
first_sample_timestamp_usec = dio{1}.first_timestamp_usec;

%% Synchronize lfp sample timestamps with TTLs
global_sample_timestamps = local2GlobalTime(ttl_timestamps_usec, local_sample_timestamps_usec);
%global_sample_timestamps = local_sample_timestamps_usec;

datFiles = dir(fullfile(path_to_recording_dir, mergedLFP_dirname, '*_*.LFP_nt*ch*.dat'));
nChannels = length(datFiles);

fnames = {datFiles.name};
parsedTokens = regexp(fnames, '(\d\d\d\d)(\d\d)(\d\d)_(\d\d\d\d\d\d).*\.LFP_nt(\d)(\d\d\d)ch\d.dat', 'tokens');

% Restructure parsed LFP files to match channel numbers with probe number
chNums = [];
probeNums = [];
for i = 1:nChannels
    probeNums = [probeNums str2num(parsedTokens{i}{1}{end-1})];
    chNums = [chNums str2num(parsedTokens{i}{1}{end})];
end

nProbes = length(unique(probeNums));

channelMap = {};
numChannels = nan([1, nProbes]);
lfpData = {};
channelIDs = {};
channelPositions = {};
for probeIdx = 1:nProbes
    % Get lfp files for probeIdx
    datFiles = dir(fullfile(path_to_recording_dir, mergedLFP_dirname, sprintf('*_*.LFP_nt%d*ch*.dat', probeIdx)));
    parsedTokens = regexp({datFiles.name}, '(\d\d\d\d)(\d\d)(\d\d)_(\d\d\d\d\d\d).*\.LFP_nt(\d)(\d\d\d)ch\d.dat', 'tokens');
    numChannels(probeIdx) = length(parsedTokens);
    channelIDs{probeIdx} = nan([1, numChannels(probeIdx)]);
    channelMap{probeIdx} = nan([1, numChannels(probeIdx)]);
    for i = 1:numChannels(probeIdx)
        % Parse Trodes channel number from file name. Smaller ID is towards
        % the tip of the probe
        channelIDs{probeIdx}(i) = str2num(parsedTokens{i}{1}{end});
        
    end
    [a, sortIdx] = sort(channelIDs{probeIdx}); % Sort in case filesystem is not in order
    sortedDatFiles = datFiles(sortIdx);
    
    % Load LFP data
    chPositions_x = xPos(chNums(probeNums == probeIdx));
    chPositions_y = yPos(chNums(probeNums == probeIdx));
    channelPositions{probeIdx} = [chPositions_x.'; chPositions_y].';
    lfp = zeros([numChannels(probeIdx),length(local_sample_timestamps_usec)],"int16");
    for idx = 1:numChannels(probeIdx)
        file = sortedDatFiles(idx);
        tokens = regexp(file.name, '(\d\d\d\d)(\d\d)(\d\d)_(\d\d\d\d\d\d).*\.LFP_nt(\d)(\d\d\d)ch(\d).dat', 'tokens');
        tokens = tokens{1};
        yyyy = tokens{1};
        mm = tokens{2};
        dd = tokens{3};
        date = [yyyy mm dd];
        hhmmss = tokens{4};
        chNum = str2num(tokens{6});
        probeNum = str2num(tokens{7});
        disp(sprintf("Loading ch %d ...", chNum));
        
        chData = readTrodesExtractedDataFile(fullfile(file.folder, file.name));
        
        lfp(idx, :) = chData.fields.data;
        channelMap{probeIdx}(idx) = chNum;
    end
    lfpData{probeIdx} = lfp;
end

out = struct();
out.timestamps = timestamps;
out.lfp = lfpData;
out.channelMap = channelMap;
out.channelID = channelIDs;
out.channelPositions = channelPositions;
out.voltage_scaling = chData.voltage_scaling;
out.local_sample_timestamps_usec = local_sample_timestamps_usec;
out.timestamp_at_creation_usec = timestamp_at_creation_usec;
out.first_timestamp_usec = first_timestamp_usec;
out.global_sample_timestamps_usec = global_sample_timestamps;
out.ttl_timestamps_usec = ttl_timestamps_usec;
out.first_sample_timestamps_usec = first_sample_timestamp_usec;

%% Show summary output and save LFP data for each probe

%=== Make a figure with the channel map for each probe
figure('units','normalized','outerposition',[.3 .1 .15 .7]);
tiledlayout(1,nProbes,'TileSpacing','compact');
for probeIdx = 1:nProbes
    nexttile;   scatter(xPos,yPos(1:end-2)',10,'k','filled');   hold on;
    scatter(out.channelPositions{1,probeIdx}(:,1),out.channelPositions{1,probeIdx}(:,2),15,'r','filled');
    xlim([-100 100]);   axis off;   title(['Probe ',num2str(probeIdx)]);
end

%=== Clean-up duplicate samples and convert global_sample_timestamps_usec into seconds
[unique_t,unique_smp] = unique(out.global_sample_timestamps_usec);
unique_t = unique_t/1e6;

%=== Generate downsampled time (2x)
ds_factor = 2;
t_ds = downsample(unique_t,ds_factor)';
smp_ds = downsample(unique_smp,ds_factor)';

%=== Downsample LFP data and save one file for each probe
red_out = struct();
for probeIdx = 1:nProbes
    
    %=== Create output structure
    red_out.xPos = xPos;                                    % All x positions along the probe
    red_out.yPos = yPos';                                   % All y positions along the probe
    red_out.lfp = out.lfp{1,probeIdx}(:,smp_ds)';           % LFP sampled at the reduced time points
    red_out.channelMap = channelMap{1,probeIdx};            % Channel map (This one should match Trodes channels!)
    red_out.channelID = channelIDs{1,probeIdx};             % Channel IDs
    red_out.channelPositions = channelPositions{1,probeIdx};% Channel positions along the probe
    red_out.voltage_scaling = chData.voltage_scaling;       % Voltage Scaling Factor
    red_out.t_ds = t_ds;                                    % Time relative to first TTL
    red_out.sampling_freq = 1/mean(diff(t_ds));             % Average sampling frequency
    
    %=== Display size of the LFP data
    info = whos('red_out');
    size_in_GB = info.bytes / (1024^3);
    disp(['LFP Data Size: ',num2str(size_in_GB), 'GB']);
    
    %=== Create analysis folder for storing the results
    analysis_directory=fullfile(path_to_recording_dir,['Sorted_units_AF']);
    if ~exist(analysis_directory,'dir');mkdir(analysis_directory);end
    disp(['Saving LFP from probe ',num2str(probeIdx)]);
    save([analysis_directory,'\LFP_probe',num2str(probeIdx),'.mat'],'red_out','-v7.3');
    
end

end

