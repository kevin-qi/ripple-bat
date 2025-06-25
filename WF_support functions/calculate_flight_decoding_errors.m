function [errors] = calculate_flight_decoding_errors (time_intervals,field_smoothed, Time_Spike_Position,position_distribution_smoothed, n_cells,position_bin_size)
    errors = [];
    for current_flight = 1: size(time_intervals,1)
        flight_st = time_intervals(current_flight,1);
        flight_en = time_intervals(current_flight,2);

        tau = 0.03; % each time bin is 30 ms 
        step = tau; % there is no overlap between time bins 
        bin_st = [flight_st:step:flight_en-tau];
        field_smoothed(find(field_smoothed==0)) = 0.0001;

        posterior = [];
        posterior_uniform_prior = [];
        for i = 1: numel(bin_st)
            time_spike_position_flight = Time_Spike_Position(find(Time_Spike_Position >= bin_st(i) & Time_Spike_Position <= bin_st(i)+tau),:);
            n_i = hist(time_spike_position_flight(:,2), [1:n_cells]);%number of spikes per cell
        %     figure;plot(n_i);xlabel("cell id"); ylabel("number of firing during flight");
            prod_term = prod(field_smoothed.^(n_i'));
            exp_term = exp(-tau*(sum(field_smoothed)));
             posterior(i,:) = exp_term.*prod_term.*position_distribution_smoothed;
             posterior_uniform_prior (i,:) = exp_term.*prod_term;
        end
        posterior_norm = normalize(posterior');
        posterior_uniform_prior_norm = normalize(posterior_uniform_prior');

        % calculating the error from actual position
        [~,max_posteriors] = max(posterior_norm, [],1);

        % plotting actual position 
        actual_positions = [];
        for i = 1: numel(bin_st)
            time_spike_position_flight = Time_Spike_Position(find(Time_Spike_Position >= bin_st(i) & Time_Spike_Position <= bin_st(i)+tau),:);
            actual_positions = [actual_positions mean(time_spike_position_flight(:,3))];
        end
        % scale actual position and position bins 
        max_posteriors_scaled = max_posteriors*position_bin_size;
        % calculate root mean error 
        errors(current_flight) = sqrt(sum((max_posteriors_scaled-actual_positions).^2)/numel(max_posteriors_scaled));
    end
    errors=mean(errors);
end
