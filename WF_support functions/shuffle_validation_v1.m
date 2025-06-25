function [weigthted_corr_p, avg_jump_p, replay_score_p, slope_p] = shuffle_validation_v1(event)
    %=== data cleaning 
    event(isnan(event)) =0;
    [weighted_corr_unshuffled, max_jump_distance_unshuffled, avg_jump_distance_unshuffled,num_peaks_unshuffled, replay_score_unshuffled, slope_unshuffled, fitted_y_unshuffled, posterior_spread_unshuffled, com_jump_unshuffled] = evaluate_candidate_event_v2(event);

    
    %=== shuffling by randperm 
    weighted_corr_test_rp = [];
    max_jump_distance_test_rp = [];
    avg_jump_distance_test_rp = [];
    num_peaks_test_rp = [];
    replay_score_test_rp =  [];
    slope_test_rp = [];
    fitted_y_test_rp = {};
    posterior_spread_test_rp = [];
    com_jump_test_rp = [];
    for i = 1:100
        prob_shuffled_bins = event(randperm(size(event,1)),:);
       [weighted_corr_test_rp(i), max_jump_distance_test_rp(i), avg_jump_distance_test_rp(i),num_peaks_test_rp(i), replay_score_test_rp(i), slope_test_rp(i), fitted_y_test_rp{i}, posterior_spread_test_rp(i), com_jump_test_rp(i)] = evaluate_candidate_event_v2(prob_shuffled_bins);
    end
    [weighted_corr_p_randperm] = p_value_WF(weighted_corr_unshuffled, weighted_corr_test_rp);
    [avg_jump_distance_p_randperm] = p_value_WF(avg_jump_distance_unshuffled, avg_jump_distance_test_rp);
    [replay_score_p_randperm] = p_value_WF(replay_score_unshuffled, replay_score_test_rp);
    [post_spread_p_randperm] = p_value_WF(posterior_spread_unshuffled, posterior_spread_test_rp);
    [slope_p_randperm] = p_value_WF(slope_unshuffled, slope_test_rp);
    [com_jump_p_randperm] = p_value_WF(com_jump_unshuffled, com_jump_test_rp);
    
    


    %=== shuffling by circshift
    weighted_corr_test_cs = [];
    max_jump_distance_test_cs = [];
    avg_jump_distance_test_cs = [];
    num_peaks_test_cs = [];
    replay_score_test_cs =  [];
    slope_test_cs = [];
    fitted_y_test_cs = {};
    posterior_spread_test_cs = [];
    com_jump_test_cs = [];
    for i = 1:100
        prob_shuffled_circshift=[];
        for ii = 1:size(event,1)
            prob_shuffled_circshift(ii,:) = circshift(event(ii,:),randi([5 size(event,2)-5]));
        end 
       [weighted_corr_test_cs(i), max_jump_distance_test_cs(i), avg_jump_distance_test_cs(i),num_peaks_test_cs(i), replay_score_test_cs(i), slope_test_cs(i), fitted_y_test_cs{i}, posterior_spread_test_cs(i), com_jump_test_cs(i)] = evaluate_candidate_event_v2(prob_shuffled_circshift);
    end 
    [weighted_corr_p_circshift] = p_value_WF(weighted_corr_unshuffled, weighted_corr_test_cs);
    [avg_jump_distance_p_circshift] = p_value_WF(avg_jump_distance_unshuffled, avg_jump_distance_test_cs);
    [replay_score_p_circshift] = p_value_WF(replay_score_unshuffled, replay_score_test_cs);
    [com_jump_p_circshift] = p_value_WF(com_jump_unshuffled, com_jump_test_cs);
    [slope_p_circshift] = p_value_WF(slope_unshuffled, slope_test_cs);
    
    
    %=== saving variables [metric_unshuffled, p_val_randperm, p_value_circshift]
    weigthted_corr_p = [weighted_corr_p_randperm weighted_corr_p_circshift];
    avg_jump_p = [avg_jump_distance_p_randperm avg_jump_distance_p_circshift];
    replay_score_p = [replay_score_p_randperm replay_score_p_circshift];
    slope_p = [slope_p_randperm slope_p_circshift];
    
        %=== visualization of randperm shuffling results 
%     fig = figure; 
%     set(gcf,'position',[10,300,1250,260]);
%     set(gca,'YDir','normal');
%     colormap(hot);
%     subplot(1,6,1); imagesc (event',prctile(event,[1 99],"all")'); title("Data before shuffling"); xlabel("Time bins"); ylabel("Position bins");
%     axis xy
%     subplot(1,6,2);imagesc(prob_shuffled_bins',prctile(event,[1 99],"all")'); title("Time bin shuffle example"); xlabel("Time bins"); ylabel("Position bins");
%     axis xy
%     subplot(1,6,3);
%     plot_p_value_WF(weighted_corr_unshuffled, weighted_corr_test_rp,weighted_corr_p_randperm)
%     xlabel("Weighted correlation"); ylabel("count");
%     title("Weighted correlation");
%     subplot(1,6,4);
%     plot_p_value_WF(avg_jump_distance_unshuffled, avg_jump_distance_test_rp,avg_jump_distance_p_randperm)
%     xlabel("Average jump distance");  ylabel("count");
%     title("Average jump distance");
%     subplot(1,6,5);
%     plot_p_value_WF(replay_score_unshuffled, replay_score_test_rp, replay_score_p_randperm);
%     xlabel("Replay score");  ylabel("count");
%     title("Replay score");
%     subplot(1,6,6);
%     plot_p_value_WF(slope_unshuffled, slope_test_rp, slope_p_randperm);
%     xlabel("slope");  ylabel("count");
%     title("slope");
%     
% %     subplot(1,6,6);
% %     plot_p_value_WF(com_jump_unshuffled, com_jump_test_rp, com_jump_p_randperm);
% %     xlabel("COM jump");  ylabel("count");
% %     title("COM jump");
% 
%     %=== visualization of circshift shuffling results 
%     fig = figure; 
%     set(gcf,'position',[10,10,1250,260]);
%     colormap(hot);
%     subplot(1,6,1); imagesc (event',prctile(event,[1 99],"all")'); title("Data before shuffling"); xlabel("Time bins"); ylabel("Position bins");
%     axis xy
%     subplot(1,6,2);imagesc(prob_shuffled_circshift',prctile(event,[1 99],"all")'); 
%     title("Data after shuffling position bins"); xlabel("Time bins"); ylabel("Position bins");
%     set(gca,'YDir','reverse');
%     axis xy
%     subplot(1,6,3);
%     plot_p_value_WF(weighted_corr_unshuffled, weighted_corr_test_cs, weighted_corr_p_circshift);
%     xlabel("Weighted correlation"); ylabel("count");
%     title("Weighted correlation");
%     subplot(1,6,4);
%     plot_p_value_WF(avg_jump_distance_unshuffled, avg_jump_distance_test_cs,avg_jump_distance_p_circshift);
%     xlabel("Average jump distance");  ylabel("count");
%     title("Average jump distance");
%     subplot(1,6,5);
%     plot_p_value_WF(replay_score_unshuffled, replay_score_test_cs, replay_score_p_circshift);
%     xlabel("Replay score");  ylabel("count");
%     title("Replay score");
%     subplot(1,6,6);
%     plot_p_value_WF(slope_unshuffled, slope_test_cs, slope_p_circshift);
%     xlabel("slope");  ylabel("count");
%     title("slope");

%     if we want to use com as a parameter 
%     subplot(1,6,6);
%     plot_p_value_WF(com_jump_unshuffled, com_jump_test_cs, com_jump_p_circshift);
%     xlabel("COM jump");  ylabel("count");
%     title("COM jump");
    
    

end