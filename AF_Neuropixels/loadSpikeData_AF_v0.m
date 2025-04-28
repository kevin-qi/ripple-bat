function [out] = loadSpikeData_AF_v0
%% Modified from loadSpikeData, originally from KQ. Feb 2024.
% To be used in the kilosort_outdir_probe folder, after manual curation

% FileList = dir(fullfile(cd,'*kilosort_outdir_probe*'));
% for i=1:3
%     cd(FileList(i).name);
%     loadSpikeData_AF_v0;
%     cd .. 
% end

folderPath = cd;
useKSLabel = 1;
%folderPath = uigetdir('','Select the folder to be analyzed');

[path_to_recording_dir,kilosort_out_folder_name] = fileparts(folderPath);

%% Load spike sorted data
if(exist(fullfile(path_to_recording_dir, kilosort_out_folder_name, 'cluster_info.tsv')))
    disp("Kilosort outputs found");
    phyFile = dir(fullfile(path_to_recording_dir,kilosort_out_folder_name, 'cluster_info.tsv'));
    if(length(phyFile)>0)
        disp("Manual curation found");
        disp('Extracting units, aligning and saving...');
        phy = readtable(fullfile(phyFile.folder, phyFile.name), 'FileType', 'text', 'VariableNamingRule', 'preserve'); % import curated units
        
        % Clusters labeled 'good' after manual curation (or just KSLabel
        if useKSLabel
            good_unit_idx = strcmp(phy.KSLabel, 'good');
            mua_unit_idx = logical(~strcmp(phy.KSLabel, 'good') .* ~strcmp(phy.KSLabel, 'noise'));
        else
            good_unit_idx = strcmp(phy.group, 'good');
            mua_unit_idx = logical(~strcmp(phy.group, 'good') .* ~strcmp(phy.group, 'noise'));
        end

        good_units = phy(good_unit_idx, {'cluster_id', 'Amplitude', 'ch', 'depth', 'fr', 'n_spikes', 'KSLabel', 'group'});
        mua_units = phy(mua_unit_idx, {'cluster_id', 'Amplitude', 'ch', 'depth', 'fr', 'n_spikes', 'KSLabel', 'group'});
        
        sp_unit = double(readNPY(fullfile(phyFile.folder, 'spike_clusters.npy'))); % Unit ID
        sp_times = double(readNPY(fullfile(phyFile.folder, 'spike_times.npy'))); % Spike times in raw sample numbers
        sp_templates = double(readNPY(fullfile(phyFile.folder, 'amplitudes.npy')));
        
        % Load TTL data
        dio = loadTrodesDigital(path_to_recording_dir);
        isRising = dio{1}.state == 1;
        local_ttl_timestamps_usec = dio{1}.ttl_timestamp_usec;
        first_sample_timestamp_usec = dio{1}.first_timestamp_usec;
        
        
        % Convert raw sample numbers by dividing by clockrate
        timeData = loadTrodesTimestamps(path_to_recording_dir);
        sp_times = sp_times(sp_times <= length(timeData.sample_timestamps_usec));
        sp_unit = sp_unit(sp_times <= length(timeData.sample_timestamps_usec));
        local_spike_times_usec = timeData.sample_timestamps_usec(sp_times);

        % Synchronize timestamps using TTLs
        global_spike_times_usec = local2GlobalTime(local_ttl_timestamps_usec, local_spike_times_usec);
        
        
        num_good_units = size(good_units ,1);
        num_mua_units = size(mua_units, 1);
 
        %% Populate data structure with good units
        out.good_units = good_units;
        spikeTimes_usec = {};
        localSpikeTimes_usec = {};
        spikeTimeSampleIdx = {};
        for i = 1:num_good_units
            spikeTimes_usec{i} = global_spike_times_usec(sp_unit == good_units.cluster_id(i));
            localSpikeTimes_usec{i} = local_spike_times_usec(sp_unit == good_units.cluster_id(i));
            spikeTimeSampleIdx{i} = sp_times(sp_unit == good_units.cluster_id(i));
        end
        out.good_units.spikeTimes_usec = spikeTimes_usec.';
        out.good_units.localSpikeTimes_usec = localSpikeTimes_usec.';
        out.good_units.spikeTimeSampleIdx = spikeTimeSampleIdx.';
        
        %% Populate data structure with mua units
        out.mua_units = mua_units;
        spikeTimes_usec = {};
        localSpikeTimes_usec = {};
        for i = 1:num_mua_units
            spikeTimes_usec{i} = global_spike_times_usec(sp_unit == mua_units.cluster_id(i));
            localSpikeTimes_usec{i} = local_spike_times_usec(sp_unit == mua_units.cluster_id(i));
        end
        out.mua_units.spikeTimes_usec = spikeTimes_usec.';
        out.mua_units.localSpikeTimes_usec = localSpikeTimes_usec.';
        
        %% Populate data fields
        out.allLocalSpikeTimes_usec = local_spike_times_usec;
        out.allGlobalSpikeTimes_usec = global_spike_times_usec;
        out.numGoodUnits = height(good_units);
        out.numMuaUnits = height(mua_units);
        out.local_ttl_timestamps_usec = local_ttl_timestamps_usec;
        out.probe_id = kilosort_out_folder_name;
        out.curation_date = datestr(now, 'yymmdd_HHMM');
    else
        disp("Manual Curation Missing. Please manually curate kilosort results first!")
        out.good_units = table();
        out.mua_units = table();
        out.numGoodUnits = 0;
        out.numMuaUnits = 0;
        out.local_ttl_timestamps_usec = [];
    end
else
    disp("No Kilosort outputs found")
    out.good_units = table();
    out.mua_units = table();
    out.numGoodUnits = 0;
    out.numMuaUnits = 0;
    out.local_ttl_timestamps_usec = [];
end

%% Save Data

%=== Create analysis folder for storing the results
analysis_directory=fullfile(path_to_recording_dir,['Sorted_units_AF']);
if ~exist(analysis_directory,'dir');mkdir(analysis_directory);end
save([analysis_directory,'\SU_',kilosort_out_folder_name,'.mat'],'out');

end

