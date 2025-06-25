%% BRIEF DESCRIPTION
tic;
for hide=1
end
%%=== LOAD DATA
for hide=1
    disp('Loading data...');
    bhvFile = dir('Extracted_Behavior_*');  load([bhvFile.folder '/' bhvFile.name]);
    imuFile = dir('IMU_data.mat');          load([imuFile.folder '/' imuFile.name]);
    nrnFile = dir('SU_kilosort*');
    NP_unit = table();
    MUA_unit = table();
    n_probes = size(nrnFile,1);
    for i=1:n_probes
        load([nrnFile(i).folder '/' nrnFile(i).name]);
        %out.good_units.group = NaN(size(out.good_units.group));                        % This is the structure containing the single units for each probe
        NP_unit = [NP_unit; out.good_units];                                            % Single units are assembled into a table
        MUA_unit = [MUA_unit; out.mua_units];
    end
    clear out;
    NP_unit.fr = cellfun(@(x) 1/mean(diff(x/1e6)),NP_unit.spikeTimes_usec);
    MUA_unit.fr = cellfun(@(x) 1/mean(diff(x/1e6)),MUA_unit.spikeTimes_usec);
    %=== Load LFP from Probe1
    load('LFP_probe1.mat');                                                             % This file contains LFP data from probe 1
    unique_ID = options.unique_ID;      
    
    %=== Create analysis folder for storing the results
    options.savefigures = 1;                                   % If creating folder for saving the data or not
    if options.savefigures == 1
        analysis_directory=fullfile(pwd);
        if ~exist(analysis_directory,'dir')
            mkdir(analysis_directory);
        end
    end

    
end 
%% loading RPL time vector 
FileList = dir(fullfile(cd,"RPL*"));
curr_RPL =load(FileList(1).name);
RPL_t_max = curr_RPL.RPL_out.time_vector(end); %lowerbound of the max
RPL_t_min = curr_RPL.RPL_out.time_vector(1);

%% loading spike time 
spk = load("spike_data.mat").spk;
spike_t_max = max(spk(:,1));
spike_t_min = min(spk(:,1));

%% behavior time 
t_max = max(t);
t_min = min(t);

%% pooling together
session_st = max([0,t_min,spike_t_min,RPL_t_min]);
session_en = min([t_max,spike_t_max,RPL_t_max]);

session_info.curr_session_id = [unique_ID{1} '_' unique_ID{2} '_' unique_ID{3}];
session_info.t = t;
session_info.r = r;
session_info.session_st = session_st;
session_info.session_en = session_en;
save([analysis_directory,'/session_info_', unique_ID{1}, '_', unique_ID{2}, '_' unique_ID{3},'.mat'],'session_info')
toc;
