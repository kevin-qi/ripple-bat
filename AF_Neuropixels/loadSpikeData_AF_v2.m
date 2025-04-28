function [out] = loadSpikeData_AF_v2
%% Modified from loadSpikeData, originally from KQ. Feb 2024.
% To be used in the kilosort_outdir_probe folder. The function considers 2 alternatives:
%   1) Load direct output from Kilosort (using KSLabel from cluster_KSLabel.tsv)
%   2) Load 'units' labeled as good after manual curation, those reported in the phy.group column of the cluster_info.tsv file

%=== Params and init
folderPath = cd;        % Current folder
useKSLabel = 1;         % Default, if using KS Label (the alternative is using manual curation)
[path_to_recording_dir,kilosort_out_folder_name] = fileparts(folderPath);

%% Load spike sorted data

%=== Check if Kilosort output exists
if(exist(fullfile(path_to_recording_dir, kilosort_out_folder_name, 'cluster_KSLabel.tsv')))
    if useKSLabel
        phyFile = dir(fullfile(path_to_recording_dir,kilosort_out_folder_name, 'cluster_KSLabel.tsv'));
        disp('Using KSLabel...');
    else
        phyFile = dir(fullfile(path_to_recording_dir,kilosort_out_folder_name, 'cluster_info.tsv'));
        disp('Using Manual Curation...');
    end
else
    disp('No Kilosort outputs found');
    return
end

%=== Proceed with extraction
disp(['Extracting units from:', phyFile.folder]);
sp_unit  = double(readNPY(fullfile(phyFile.folder, 'spike_clusters.npy')));     % Cluster assigned to each spike (modified after each run of Phy)
sp_times = double(readNPY(fullfile(phyFile.folder, 'spike_times.npy')));        % Spike times in raw sample numbers
sp_amp   = double(readNPY(fullfile(phyFile.folder, 'amplitudes.npy')));         % Spike amplitude (template scaled)
sp_pos   = double(readNPY(fullfile(phyFile.folder, 'spike_positions.npy')));    % Spike positions on the probe
templates   = double(readNPY(fullfile(phyFile.folder, 'templates.npy')));       % Templates for each cluster (n_cluster, samples, channels)

%=== Get the IDs of good and mua clusters
phy = readtable(fullfile(phyFile.folder, phyFile.name), 'FileType', 'text', 'VariableNamingRule', 'preserve'); 
if useKSLabel  
    good_unit_idx = strcmp(phy.KSLabel, 'good');
    mua_units_idx = strcmp(phy.KSLabel, 'mua');
else
    good_unit_idx = strcmp(phy.group, 'good');
    mua_units_idx = strcmp(phy.group, 'mua');
end
good_units = phy(good_unit_idx, {'cluster_id'});
mua_units =  phy(mua_units_idx, {'cluster_id'});
num_good_units = size(good_units ,1);
num_mua_units = size(mua_units ,1);

disp([num2str(num_good_units),' good units found']);
disp([num2str(num_mua_units),' MUA units found']);


%% Synchronization

%=== Load TTL data
dio = loadTrodesDigital(path_to_recording_dir);
local_ttl_timestamps_usec = dio{1}.ttl_timestamp_usec;

%=== Convert raw sample numbers by dividing by clockrate
timeData = loadTrodesTimestamps(path_to_recording_dir);
sp_times = sp_times(sp_times <= length(timeData.sample_timestamps_usec));
sp_unit = sp_unit(sp_times <= length(timeData.sample_timestamps_usec));
local_spike_times_usec = timeData.sample_timestamps_usec(sp_times);

%=== Synchronize timestamps using TTLs
global_spike_times_usec = local2GlobalTime(local_ttl_timestamps_usec, local_spike_times_usec);

%% Populate data structure with good units and mua units

%=== Start with good units
out.good_units = good_units;
spikeTimes_usec = {};
localSpikeTimes_usec = {};
spikeTimeSampleIdx = {};
spikePos_um = {};
spikeTemplate = {};
for i = 1:num_good_units
    spikeTimes_usec{i} = global_spike_times_usec(sp_unit == good_units.cluster_id(i));
    localSpikeTimes_usec{i} = local_spike_times_usec(sp_unit == good_units.cluster_id(i));
    spikeTimeSampleIdx{i} = sp_times(sp_unit == good_units.cluster_id(i));
    spikePos_um{i} = sp_pos(sp_unit == good_units.cluster_id(i),:);
    spikeTemplate{i} = squeeze(templates(good_units.cluster_id == good_units.cluster_id(i),:,:));
end
out.good_units.spikeTimes_usec = spikeTimes_usec.';             % Timestamps of the spikes relative to first TTL
out.good_units.localSpikeTimes_usec = localSpikeTimes_usec.';   % Original Timestamps of the spikes
out.good_units.spikePos_um = spikePos_um.';                     % x,y positions of the spike on the probe
out.good_units.template = spikeTemplate.';                      % Template for each cluster

%=== Add mua units
out.mua_units = mua_units;
spikeTimes_usec = {};
localSpikeTimes_usec = {};
spikeTimeSampleIdx = {};
spikePos_um = {};
spikeTemplate = {};
for i = 1:num_mua_units
    spikeTimes_usec{i} = global_spike_times_usec(sp_unit == mua_units.cluster_id(i));
    localSpikeTimes_usec{i} = local_spike_times_usec(sp_unit == mua_units.cluster_id(i));
    spikeTimeSampleIdx{i} = sp_times(sp_unit == mua_units.cluster_id(i));
    spikePos_um{i} = sp_pos (sp_unit == mua_units.cluster_id(i),:);
    spikeTemplate{i} = squeeze(templates(mua_units.cluster_id == mua_units.cluster_id(i),:,:));
end
out.mua_units.spikeTimes_usec = spikeTimes_usec.';             % Timestamps of the spikes relative to first TTL
out.mua_units.localSpikeTimes_usec = localSpikeTimes_usec.';   % Original Timestamps of the spikes
out.mua_units.spikePos_um = spikePos_um.';                     % x,y positions of the spike on the probe
out.mua_units.template = spikeTemplate.';                      % Template for each cluster

%=== Add a few more variables
out.allLocalSpikeTimes_usec = local_spike_times_usec;           % Original time stamps
out.allGlobalSpikeTimes_usec = global_spike_times_usec;         % Timestamps after alignment
out.numGoodUnits = height(good_units);                          % Number of good units
out.numMUAUnits = height(mua_units);                            % Number of Mua units
out.local_ttl_timestamps_usec = local_ttl_timestamps_usec;      % Timestamps of the 3s TTLs
out.probe_id = kilosort_out_folder_name;                        % Name of the folder where results were loaded
out.curation_date = datestr(now, 'yymmdd_HHMM');                % Time this out was generated 
out.useKSLabel = useKSLabel;                                    % Flag for using KSLabel or output of manual curation

%% Save Data

%=== Create analysis folder for storing the results
analysis_directory=fullfile(path_to_recording_dir,['Sorted_units_AF']);
if ~exist(analysis_directory,'dir');mkdir(analysis_directory);end
save([analysis_directory,'\SU_',kilosort_out_folder_name,'.mat'],'out');

end
