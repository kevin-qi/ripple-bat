function [new_st, new_en] = horiz_segment_v1(event)
    %=== setting up event for testing 
    test = 0;
    if test 
       event = RP.post{5726} 
    end
    nS = size(event,2);
    nT = size(event,1);
    uniform_prob = 1/nS;
    horiz_perc = 0.05;                                                      % if position moved less than this percent in the selected timebins, then erase
   
    %     min_dur = ceil(nT*0.25);  
    min_dur = 0;
% segment only horizontal if it's horizontal by at least this long
    %=== setting each time vector as horizontal or not (<3*uniform
    %distribution, before/after nan)
    event(find(isnan(event)))=1/size(event,2);
    [position_prob, decoded_position] = max(event');
    rescale_pos =  rescale([1:size(event,2)]);
    predicted_pos = rescale(decoded_position);
    predicted_pos_diff = diff(predicted_pos);
    predicted_pos_diff = [predicted_pos_diff predicted_pos_diff(end)];

    
    horiz_segment = position_prob<3*uniform_prob;                              % setting low uniform position to be horizontal
%     horiz_segment = conv2(horiz_segment, [1 1 1], "same");                        % setting before and after nan positions to be horizontal
%     horiz_segment (abs(predicted_pos_diff)<= horiz_perc) = 1;                  % setting adjacent ones to horizontal
    event_time_vector = [1:numel(decoded_position)];    
        

    horiz_segment=horiz_segment>0;

    % calculate percent of time bins that are horizontal
%     horiz_sum = sum(horiz_segment);
%     horizontal_percent = horiz_sum/size(event,1);
    
    %
    %=== filter out events that are too short
    event_st = strfind(~horiz_segment, [0 1]);
    event_en = strfind(~horiz_segment, [1 0]);
    if ~horiz_segment(1)==1
        event_st = [1 event_st];
    end
    if ~horiz_segment(end)==1
        event_en = [event_en numel(horiz_segment)];
    end
    event_dur = event_en - event_st; 
    min_dur = 5; %25 ms
    event_erase_idx = find(event_dur<min_dur); 
    for i = event_erase_idx
        horiz_segment([event_st(i): event_en(i)])=1;
    end


    
    % calculate start and end of cut 
    horiz_binary = horiz_segment>0;
    if horiz_binary(1) == 0 
       new_st = 1;
    end
    if horiz_binary(end)==0
       new_en = nT;
    end
    if horiz_binary(1) == 1
       new_st = strfind(horiz_binary, [1 0]);
       if numel(new_st)>1
       	new_st = new_st(1);
       end
    end
    if horiz_binary(end)==1
       new_en = strfind(horiz_binary, [0 1]);
       if numel(new_en)>1
        new_en = new_en(end);
       end
    end
    if test 
        figure('units','normalized','outerposition',[.5 0.5 .3 .3]);
       subplot(1,2,1);
        imagesc(event');
        axis  xy;
%         hold on;scatter(event_time_vector,decoded_position); 
        subplot(1,2,2);
        imagesc(event([new_st:new_en],:)')
                axis  xy;
        colormap(hot);
    end
 end 