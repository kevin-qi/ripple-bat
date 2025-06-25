function [weighted_corr] = calc_weighted_corr(event)
    event(find(isnan(event)))=1/size(event,2);
    [position_prob, decoded_position] = max(event');
    event_time_vector = [1:numel(decoded_position)];
    weighted_mean_time = sum(position_prob.*event_time_vector)/sum(position_prob);
    weighted_mean_position = sum(position_prob.*decoded_position)/sum(position_prob);
    weighted_covariance_x_y = sum(position_prob.*(event_time_vector-weighted_mean_time).*(decoded_position-weighted_mean_position))/sum(position_prob);
    weighted_covariance_x = sum(position_prob.*(event_time_vector-weighted_mean_time).^2)/sum(position_prob);
    weighted_covariance_y = sum(position_prob.*(decoded_position-weighted_mean_position).^2)/sum(position_prob);
    weighted_corr = weighted_covariance_x_y/sqrt(weighted_covariance_x*weighted_covariance_y); 
end