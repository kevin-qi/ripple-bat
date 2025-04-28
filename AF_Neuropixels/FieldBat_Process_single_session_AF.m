function FieldBat_Process_single_session_AF()

%% Brief script for Exporting and Kilosorting NP data

%=== Run this script in the folder containing the 3 standard items:
% yyyymmdd_hhmmss.rec                       Folder with the locally recorded Trodes data merged with SD card data
% batid_day_env_workspace.trodesconf        File for Env Workspace
% batid_day_neural_workspace.trodesconf     File for Neural Workspace

%=== The first section [Export Trodes Data] extracts the relevant data from the merged files
%=== The second section [Run Kilosort] runs Kilosort 3.5 on each of the 3 probes

%% Export Trodes Data
trodesPath = 'C:\Users\Angelo\Desktop\Software_AF\Trodes_2-5-0_Windows64\Trodes_2-5-0_Windows64';    %=== Adjust this one depending on the workstation you are using
tic;
disp("Exporting data...")
exportAllTrodesRecordings(cd,trodesPath,...
    'extractLFP', 1, ...
    'extractKilosort', 1, ...
    'extractSpikes', 1, ...
    'extractSpikeBand', 0, ...
    'extractTime', 1, ...
    'extractDIO', 1, ...
    'extractAnalog', 1);
toc;

%% Run Kilosort

% %=== Get the name of the trodesconf file, to count the number of channels from each probe
% trodes_configuration_file = dir('*neural*.trodesconf');
% trodes_configuration = readTrodesFileConfig(trodes_configuration_file.name);
% 
% %=== Get the name of the .rec folder and go inside
% dot_rec_folder = dir('*.rec');
% rec_folder_path = fullfile(pwd,dot_rec_folder.name);
% cd(rec_folder_path);
% 
% %=== Get the channels for each probe
% numch = zeros(1,3);
% all_channels = cellfun(@(x) str2double(x(2:end)),{trodes_configuration.nTrodes.id},'UniformOutput',false);
% all_probes = cellfun(@(x) str2double(x(1)),{trodes_configuration.nTrodes.id},'UniformOutput',false);
% for i=1:3
%     numch(i) = sum([all_probes{:}] == i);
% end
% 
% %=== Run Kilosort on each probe
% tot_time_kilo = tic;
% path_to_exp_dir = cd;
% disp("Kilosorting Data...")
% runKilosort3(path_to_exp_dir, 'probeNum', 1, 'numChannels', numch(1), 'nBlocks', 1);
% runKilosort3(path_to_exp_dir, 'probeNum', 2, 'numChannels', numch(2), 'nBlocks', 1);
% runKilosort3(path_to_exp_dir, 'probeNum', 3, 'numChannels', numch(3), 'nBlocks', 1);
% toc(tot_time_kilo);

end
