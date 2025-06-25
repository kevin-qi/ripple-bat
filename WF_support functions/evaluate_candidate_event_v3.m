function [weighted_corr, max_jump_dist, avg_jump_distance,num_peaks, replay_score, slope, fitted_y, posterior_spread, com_jump] = evaluate_candidate_event_v3(event)
    try
    event(find(isnan(event)))=1/size(event,2);
    
    [position_prob, decoded_position] = max(event');
    event_time_vector = [1:numel(decoded_position)];
    weighted_mean_time = sum(position_prob.*event_time_vector)/sum(position_prob);
    weighted_mean_position = sum(position_prob.*decoded_position)/sum(position_prob);
    weighted_covariance_x_y = sum(position_prob.*(event_time_vector-weighted_mean_time).*(decoded_position-weighted_mean_position))/sum(position_prob);
    weighted_covariance_x = sum(position_prob.*(event_time_vector-weighted_mean_time).^2)/sum(position_prob);
    weighted_covariance_y = sum(position_prob.*(decoded_position-weighted_mean_position).^2)/sum(position_prob);
    weighted_corr = weighted_covariance_x_y/sqrt(weighted_covariance_x*weighted_covariance_y); 
    
    max_jump_dist = max(abs(diff(decoded_position)));
    y_hist = hist(decoded_position,10);
    [pks_y,pks_x] = findpeaks(y_hist);
    num_peaks = numel(pks_y);
    
    uniform_prob = 1/size(event,2);
    %=== rescaling y 
    event_position_vector = [1:size(event,2)];
    event_pos_vector_rescaled = rescale (event_position_vector);
    event_pos_scaled = event_pos_vector_rescaled(decoded_position);

    %=== computing average jump distance 
    avg_jump_distance = mean(abs(diff(event_pos_scaled)));
    
    pos_potential_events = decoded_position;
    pos_potential_events(find(decoded_position==0))= NaN;
    max_pos_diff_potential_events = nan(size(pos_potential_events));
    max_pos_diff_potential_events(find(isnan(pos_potential_events)==0))=[diff(pos_potential_events(find(isnan(pos_potential_events)==0))),NaN];
    %=== computing average posterior_spread and average COM jump
    post_spread=[];
    post_COM =[];
    for i =  1:size(event,1) % loop through every time point 
        current_post = event(i,:);
        [~,predicted_pos] = max(current_post'); 
        predicted_pos = event_pos_vector_rescaled(predicted_pos);
        post_spread(i) = sqrt(sum(((event_pos_vector_rescaled- predicted_pos).^2).*current_post));
        post_COM(i) = sum(current_post.*event_pos_vector_rescaled);
    end 
    posterior_spread = sum(post_spread)/size(event,1);
    com_jump = sum(abs(diff(post_COM)))/size(event,1);

    %=== calculating replay score 
    % fitting a line
    line_exclude = 0;
    max_posterior_time_pos =[];
    max_posterior_time_pos(:,1) = event_time_vector;
    max_posterior_time_pos(:,2) = event_pos_scaled;
    max_posterior_time_pos_filtered = max_posterior_time_pos(find( position_prob' > (3*uniform_prob) & max_posterior_time_pos(:,2) > line_exclude+ min(max_posterior_time_pos(:,2)) & max_posterior_time_pos(:,2) <  -line_exclude + max(max_posterior_time_pos(:,2))),:);
%     max_posterior_time_pos_filtered = max_posterior_time_pos;
    
    b = robustfit(max_posterior_time_pos_filtered(:,1),max_posterior_time_pos_filtered(:,2));
%     b = robustfit(max_posterior_time_pos(:,1),max_posterior_time_pos(:,2));

    line_slope_intercept = [b(2) b(1)]; 
    
    fitted_y = polyval(line_slope_intercept,event_time_vector);
    
    % calculating the average decoded probability within range d% of fitted line
    d=0.1;
    avg_prob_within_d = [];
    for i = 1:numel(fitted_y)
        if (event_pos_scaled(i)==0) % introduced to specifically penalize the time bins with uniform prob
            avg_prob_within_d(i) = 0;
%         elseif (fitted_y(i) < min(event_pos_vector_rescaled)+d |  fitted_y(i) > max(event_pos_vector_rescaled)-d)
%             avg_prob_within_d (i) = mean(event(i,:));

        elseif  (fitted_y(i) < min(event_pos_vector_rescaled)+d)
            range_low =  min(event_pos_vector_rescaled);
            range_high = d;
            current_prob = event(i,:);
            range_prob = current_prob (find (event_pos_vector_rescaled > range_low & event_pos_vector_rescaled < range_high));
%             avg_prob_within_d (i) = mean(range_prob);    
              avg_prob_within_d (i) = sum(range_prob);    
        elseif  (fitted_y(i) > max(event_pos_vector_rescaled)-d)
            range_low =  max(event_pos_vector_rescaled)-d;
            range_high = max(event_pos_vector_rescaled);
            current_prob = event(i,:);
            range_prob = current_prob (find (event_pos_vector_rescaled > range_low & event_pos_vector_rescaled < range_high));
%             avg_prob_within_d (i) = mean(range_prob); 
            avg_prob_within_d (i) = sum(range_prob);    
        else 
            range_low = fitted_y(i) -d;
            range_high = fitted_y(i) + d;
            current_prob = event(i,:);
            range_prob = current_prob (find (event_pos_vector_rescaled > range_low & event_pos_vector_rescaled < range_high));
%             avg_prob_within_d (i) = mean(range_prob);
            avg_prob_within_d (i) = sum(range_prob);    
        end
    end 
      % visualization of max position and the line of best fit, and the average posterior
      % within +-d, at each time point
%     figure; imagesc(event_time_vector, event_pos_vector_rescaled,event');
%     axis xy;
%     hold on; scatter (max_posterior_time_pos(:,1), max_posterior_time_pos(:,2),"filled",'MarkerFaceColor',[1 1 1])
%     hold on; scatter (max_posterior_time_pos_filtered(:,1), max_posterior_time_pos_filtered(:,2),"filled",'MarkerFaceColor',[1 1 1])
%     hold on; plot(event_time_vector,fitted_y,'-','Color',[1 1 1])
%     plot(avg_prob_within_d); % visualize average prob on top of plot 
    replay_score = mean(avg_prob_within_d);
    slope = line_slope_intercept(1);
    catch
        weighted_corr = NaN; 
        max_jump_dist = NaN; 
        avg_jump_distance = NaN;
        num_peaks = NaN;
        replay_score  = NaN;
        slope = NaN; 
        fitted_y = NaN;
        posterior_spread = NaN;
        com_jump = NaN;
    end
    
    % calculating replay score 
end

% 
% % %     % visualization of max position and the line of best fit
%     figure; imagesc(event_time_vector, event_pos_vector_rescaled,event');
%     axis xy;
%     hold on; scatter (max_posterior_time_pos(:,1), max_posterior_time_pos(:,2),"filled",'MarkerFaceColor',[1 1 1])
%     hold on; scatter (max_posterior_time_pos_filtered(:,1), max_posterior_time_pos_filtered(:,2),"filled",'MarkerFaceColor',[1 1 0])
%     hold on; scatter (event_time_vector, post_COM,"filled",'MarkerFaceColor',[1 0 1])
%     hold on; plot(event_time_vector,fitted_y,'-','Color',[1 1 1])
%     hold on; plot(event_time_vector,[0 abs(diff(post_COM))],'-','Color',[1 1 0])
%     
% %     hold on; plot(event_time_vector,[0 abs(diff( max_posterior_time_pos(:,2)))'],'-','Color',[1 1 0])
%     
%     plot(avg_prob_within_d); % visualize average prob on top of plot 
%     colormap(hot)