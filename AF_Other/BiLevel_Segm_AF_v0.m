function [state_vector,state_num,state_str,state_stp,width] = BiLevel_Segm_AF_v0(signal,threshold,Fs,db_time,min_dur)
%% Theshold input signal and debounce

% ---- INPUT -----
% signal:       input signal
% threshold:    threshold for HI/LOW determination
% Fs:           sampling Frequency
% db_time:      time for joining adjacent intervals
% min_dur:      minimal duration of an event

% ---- OUTPUT ----
% state_vector: logical vector of 0s (LOWs) and 1s (HIs)  
% state_num:    number of HI states
% state_str:    start sample of HI states
% state_stp:    stop sample of HI states
% width:        duration of the HIs in seconds

%=== Make sure signal is a column vector
if isrow(signal), signal = signal'; end

% Initialize and define vectors 
state_vector = false(size(signal));                     % Vector with the final state value (LO/HI)
raw_state = signal > threshold;                         % Signal above threshold
raw_max_r = movmax(raw_state,[0 round(Fs*db_time)]);    % Left movmax
raw_max_l = movmax(raw_state,[round(Fs*db_time) 0]);    % Right movmax
jnd_state = raw_max_r & raw_max_l;                      % Debounced signal

% Find HI states and keep only longer than min_dur
[width,t_str,t_stp] = pulsewidth(double(jnd_state),Fs); % Find HIs
t_str = t_str(width>min_dur); state_str = round(t_str*Fs);  % Keep only HIs longer than min_dur
t_stp = t_stp(width>min_dur); state_stp = round(t_stp*Fs);  % Keep only HIs longer than min_dur
width = width(width>min_dur);

% Fill-up the state_vector
state_num = numel(state_str);                               % Number of HIs
for f = 1:state_num
    state_vector(state_str(f,1):state_stp(f,1),1) = 1;
end

% %=== Uncomment to debug
% t = linspace(0,numel(signal)/Fs,numel(signal))';
% figure('units','normalized','outerposition',[0 .2 1 .5]);
% tiledlayout(4,1,'TileSpacing','none');
% ax(1) = nexttile;   plot(t,signal);   hold on;    refline(0,threshold);
% ax(2) = nexttile;   area(t,raw_state);   
% ax(3) = nexttile;   area(t,jnd_state);
% ax(4) = nexttile;   area(t,state_vector);
% linkaxes(ax,'x');

end
