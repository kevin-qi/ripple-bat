function [weighted_corr, max_jump_dist, avg_jump_distance,num_peaks] = evaluate_candidate_event(event)
    [position_prob, decoded_position] = max(event');
    event_time_vector = [1:numel(decoded_position)];
    weighted_mean_time = sum(position_prob.*event_time_vector)/sum(position_prob);
    weighted_mean_position = sum(position_prob.*decoded_position)/sum(position_prob);
    weighted_covariance_x_y = sum(position_prob.*(event_time_vector-weighted_mean_time).*(decoded_position-weighted_mean_position))/sum(position_prob);
    weighted_covariance_x = sum(position_prob.*(event_time_vector-weighted_mean_time).^2)/sum(position_prob);
    weighted_covariance_y = sum(position_prob.*(decoded_position-weighted_mean_position).^2)/sum(position_prob);
    weighted_corr = weighted_covariance_x_y/sqrt(weighted_covariance_x*weighted_covariance_y); 
    max_jump_dist = max(abs(diff(decoded_position)));
    avg_jump_distance = mean(abs(diff(decoded_position)));
    y_hist = hist(decoded_position,[1:5:size(event,2)]);
    [pks_y,pks_x] = findpeaks(y_hist);
    num_peaks = numel(pks_y);
end