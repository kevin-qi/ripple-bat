function [new_st, new_en] = segment_replay_v0(event)
    event(find(isnan(event)))=1/size(event,2);
    [position_prob, decoded_position] = max(event');
    event_time_vector = [1:numel(decoded_position)];
    % figure;
    % imagesc(event');
    % hold on;scatter(event_time_vector,decoded_position);

    for not_used = []

        segment_replay_v1
    % weighted_corr_test_L = [];
    % weighted_corr_test_L_x = [];
    % === cut off on the left
    % for i= 1:5:numel(decoded_position)
    %     weighted_corr_test_L_x = [weighted_corr_test_L_x, i]
    %     weighted_corr_test_L=[weighted_corr_test_L,calc_weighted_corr(event(i:end,:))];
    % end
    % hold on; plot(weighted_corr_test_L_x,abs(weighted_corr_test_L)*size(event,2))
    % === cut off on the right
    % weighted_corr_test_R = [];
    % weighted_corr_test_R_x = [];
    % for i= 1:5:numel(decoded_position)
    %     weighted_corr_test_R_x = [weighted_corr_test_R_x, i]
    %     weighted_corr_test_R=[weighted_corr_test_R,calc_weighted_corr(event(1:i,:))];
    % end
    % hold on; plot(weighted_corr_test_R_x,abs(weighted_corr_test_R)*size(event,2))
    end

    %=== paired random int 
    min_bin_interval = 10;
    if numel(event_time_vector) <= min_bin_interval
        new_st = 0;
        new_en = 0;
        return;
    end

    cut_st = [];
    cut_en = [];
    cut_wc = [];
    for n = 1:1000
        checking = 1;
        while checking 
            r1 = randi(size(event,1),1);
            r2 = randi(size(event,1),1);
            if r2-r1>=min_bin_interval
                break
            end
        end
        cut_st = [cut_st, r1];
        cut_en = [cut_en, r2];
        cut_wc = [cut_wc, calc_weighted_corr(event(r1:r2,:))];
    end 
    % figure; hist(cut_wc);
    top_prctile = prctile(abs(cut_wc),95);
    top_idx = find(abs(cut_wc)>top_prctile);
    % top_idx = find(abs(cut_wc)>0.85);
    
    cut_dur = cut_en (top_idx) - cut_st(top_idx);

    % using the segment with the maximum duration
    % [~,chosen_idx] = max(cut_dur);
    if isempty(top_idx)
        new_st = 0;
        new_en = 0; 
        return;
    end

    % using the segment parameters randomly 
    chosen_idx = randi(numel(top_idx));

    new_st = cut_st(top_idx(chosen_idx));
    new_en = cut_en(top_idx(chosen_idx));


    %=== plotting, to be commented out while use
%     figure('units','normalized','outerposition',[.2 .3 .6 .55]);
%     subplot(1,2,1); imagesc(event')
%     subplot(1,2,2); imagesc(event(new_st:new_en,:)')
%     colormap(hot);
    % figure; hist(cut_en(top_idx));
end 