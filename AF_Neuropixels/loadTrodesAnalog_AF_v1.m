function [sensors] = loadTrodesAnalog_AF_v1
%% Modified from loadTrodesAnalog, originally from KQ. Feb 2024.
% To be used in the .rec folder

%=== Get file names and paths
path_to_recording_dir = cd;
[~, dirname, ext] = fileparts(path_to_recording_dir);
assert(strcmp(ext, ".rec"), "Trodes recording directory must end in .rec");
mergedAnalog_dirname = dirname + "_merged.analog";
mergedAnalog_filename = dirname + "_merged.timestamps.dat";

%% Analog data sample timestamps
disp(['Extracting analog signals from:', mergedAnalog_dirname]);
timestampData = readTrodesExtractedDataFile(fullfile(path_to_recording_dir, mergedAnalog_dirname, mergedAnalog_filename));
first_sample_timestamp_usec = double(timestampData.first_timestamp)*1e6/timestampData.clockrate;
sample_timestamps = timestampData.fields.data;
sample_timestamps_usec = 1e6 * double(sample_timestamps) / double(timestampData.clockrate);

%% TTL data
dio = loadTrodesDigital(path_to_recording_dir);
isRising = dio{1}.state == 1;
ttl_timestamps_usec = dio{1}.ttl_timestamp_usec;
global_sample_timestamps_usec = local2GlobalTime(ttl_timestamps_usec, sample_timestamps_usec);
close all;

%% Read and organize Acceleromer and Gyroscope data structure

% 3-axis accelerometer – The accelerometer range is +/- 2g, encoded in signed 16-bit integers: each step is 2/32767g, 0.000061g.
% 3-axis gyro – The gyro range is +/- 2000 degrees/second at 16 bits: 0.061deg/sec per step.

sensorNames = ["Accel", "Gyro"];
euclideanAxes = ["X", "Y", "Z"];
sensors = struct();
for i = 1:length(sensorNames)
    for j = 1:length(euclideanAxes)
        fname = join([mergedAnalog_dirname, "_Headstage_", sensorNames(i), euclideanAxes(j), ".dat"],"");
        data = readTrodesExtractedDataFile(fullfile(path_to_recording_dir, mergedAnalog_dirname, fname));
        sensorVoltage = data.fields.data;
        sensors.(lower(sensorNames(i)))(j,:) = sensorVoltage;
        sensors.raw.(lower(sensorNames(i)))(j,:)  = data;
    end
end
sensors.global_sample_timestamps_usec = global_sample_timestamps_usec;
sensors.local_sample_timestamps_usec = sample_timestamps_usec;

%% Save a downsampled version of the data

%=== Clean-up duplicate samples and convert into seconds
[unique_t,unique_smp] = unique(sensors.global_sample_timestamps_usec);
unique_t = unique_t/1e6;

%=== Generate downsampled time at 500 Hz (which corresponds to the original sampling frequency, according to SpikeGadgets Manual)
tot_time = (unique_t(end)-unique_t(1));
t_ds = linspace(unique_t(1),unique_t(end),round(tot_time*500));

%=== Extract acceleration and gyroscope data and interpolate
raw_accel = double(sensors.accel(:,unique_smp))'.*0.000061;
raw_gyros = double(sensors.gyro(:,unique_smp))'.*0.061;
acc_ds = interp1(unique_t,raw_accel,t_ds);
gyr_ds = interp1(unique_t,raw_gyros,t_ds);

%===Create structure
NP_imu.acc = acc_ds;                            % Accelerometer
NP_imu.gyr = gyr_ds;                            % Gyroscope
NP_imu.t = t_ds;                                % Downsampled time vector
NP_imu.Fs = 1/mean(diff(t_ds));                 % Sampling Frequency

%=== Create analysis folder for storing the results
analysis_directory=fullfile(path_to_recording_dir,['Sorted_units_AF']);
if ~exist(analysis_directory,'dir');mkdir(analysis_directory);end
save([analysis_directory,'\IMU_data.mat'],'NP_imu');



