% NP_replay_detect_all_sessions

%=== Look for all folders
recFolders = dir('Dataset*');  % Look for .rec filenames
recFolders =  recFolders([recFolders.isdir]);    % Keep only folders
%=== Begin Processing...
for file_number = 1:length(recFolders)
    cd(recFolders(file_number).folder);
    cd(recFolders(file_number).name);
    disp(['Processing ', recFolders(file_number).name]);
    try
        NP_continuous_replay_detection_final()
    catch ME
        disp(['There was an error with session ' recFolders(file_number).name]);
        fprintf(2,'The message was:\n%s',ME.message);
        fprintf('\n');
    end
end

