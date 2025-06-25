function [] = shuffle_validation(event, name)

    %=== shuffling times vs weighted corr and max jump distance  (to see that
    %value of metric go down with each iteration of shuffling)
    % shuffle_times = 10;
    % event_shuffled = event; 
    % weighted_corr_test = [];
    % % max_jump_distance_test = [];
    % avg_jump_distance_test = [];
    % % num_peaks_test = [];
    % for i = 1:shuffle_times 
    %     first_col = randi([1,30]);
    %     second_col = randi([1,30]);
    %     placehold_col = event_shuffled(first_col,:);
    %     event_shuffled(first_col,:) = event_shuffled(second_col,:);
    %     event_shuffled(second_col,:) = placehold_col;
    %     [weighted_corr_test(i),max_jump_distance_test(i), avg_jump_distance_test(i),num_peaks_test(i)] = evaluate_candidate_event(event_shuffled);
    % end 
    % figure; 
    % subplot(2,2,1); imagesc(event_shuffled'); xlabel("time bins"); ylabel("position");
    % subplot(2,2,2); plot(weighted_corr_test); xlabel("number of shuffles"); ylabel("weighted correlation");
    % subplot(2,2,3); plot(max_jump_distance_test);  xlabel("number of shuffles"); ylabel("max jump distance");
    % subplot(2,2,4); plot(avg_jump_distance_test);  xlabel("number of shuffles"); ylabel("avg jump distance");
    %=== data cleaning 
    event(isnan(event)) =0;
    [weighted_corr_unshuffled,max_jump_distance_unshuffled, avg_jump_distance_unshuffled, num_peaks_unshuffled] = evaluate_candidate_event(event);

    %=== shuffling by randperm 
    weighted_corr_test_rp = [];
    max_jump_distance_test_rp = [];
    avg_jump_distance_test_rp = [];
    num_peaks_test_rp = [];
    for i = 1:1000 
        prob_shuffled_bins = event(randperm(size(event,1)),:);
        % figure; imagesc(prob_shuffled_bins');
        [weighted_corr_test_rp(i),max_jump_distance_test_rp(i), avg_jump_distance_test_rp(i),num_peaks_test_rp(i)] = evaluate_candidate_event(prob_shuffled_bins);
    end 
%     fig = figure;
%     subplot(1,2,1); hist(weighted_corr_test); xlabel("weighted correlation from shuffles"); ylabel("count");
%     subplot(1,2,2); hist(avg_jump_distance_test);  xlabel("avg jump distance");  ylabel("count");
%     ax=axes(fig,'visible','off'); ax.Title.Visible='on'; 
%     title (ax,"generated data shuffled 1000x (bin shuffle)"); 

    fig = figure; 
    set(gcf,'position',[10,300,1250,260]);
    set(gca,'YDir','normal');
    colormap(hot);
    subplot(1,4,1); imagesc (event',prctile(event,[1 99],"all")'); title("Data before shuffling"); xlabel("Time bins"); ylabel("Position bins");
    axis xy
    subplot(1,4,2);imagesc(prob_shuffled_bins',prctile(event,[1 99],"all")'); title("Time bin shuffle example"); xlabel("Time bins"); ylabel("Position bins");
    axis xy
    subplot(1,4,3);[weighted_corr_p_randperm_l,weighted_corr_p_randperm_r] = plot_shuffled([weighted_corr_unshuffled weighted_corr_test_rp]);
    xlabel("Weighted correlation"); ylabel("count");
    title("Weighted correlation(shuffled 1000x)");
    subplot(1,4,4);
    [avg_jump_distance_p_randperm_1,avg_jump_distance_p_randperm_r] = plot_shuffled([avg_jump_distance_unshuffled avg_jump_distance_test_rp]);
    xlabel("Average jump distance");  ylabel("count");
    title("Average jump distance (shuffled 1000x)");
     saveas(gcf, ['validation_fig/validate_synthetic_randperm' name '.png']);
%     ax=axes(fig,'visible','off'); ax.Title.Visible='on'; 
%     title (ax,"Time bins shuffled 1000x"); 
    
    %=== shuffling by circshift
    weighted_corr_test = [];
    max_jump_distance_test = [];
    avg_jump_distance_test = [];
    num_peaks_test = [];
    prob_shuffled_circshift = [];
    for i = 1:1000
        for ii = 1:size(event,1)
            prob_shuffled_circshift(ii,:) = circshift(event(ii,:),randi([5 size(event,2)-5]));
        end 
        % figure; imagesc(prob_shuffled_circshift');
        [weighted_corr_test(i),max_jump_distance_test(i), avg_jump_distance_test(i),num_peaks_test(i)] = evaluate_candidate_event(prob_shuffled_circshift);
    end 

%     fig = figure; 
%     subplot(1,2,1); hist(weighted_corr_test); xlabel("weighted correlation from shuffles"); ylabel("count");
%     subplot(1,2,2); hist(avg_jump_distance_test);  xlabel("avg jump distance");  ylabel("count");
%     ax=axes(fig,'visible','off'); ax.Title.Visible='on'; 
%     title (ax,"generated data shuffled 1000x (peak shuffle)"); 

    fig = figure; 
    set(gcf,'position',[10,10,1250,260]);
    colormap(hot);
    subplot(1,4,1); imagesc (event',prctile(event,[1 99],"all")'); title("Data before shuffling"); xlabel("Time bins"); ylabel("Position bins");
    axis xy
    subplot(1,4,2);imagesc(prob_shuffled_circshift',prctile(event,[1 99],"all")'); 
    title("Data after shuffling position bins"); xlabel("Time bins"); ylabel("Position bins");
    set(gca,'YDir','reverse');
    axis xy
    subplot(1,4,3);[weighted_corr_p_randperm_l,weighted_corr_p_randperm_r] = plot_shuffled([weighted_corr_unshuffled weighted_corr_test]);
    xlabel("Weighted correlation"); ylabel("count");
    title("Weighted correlation shuffeled 1000x");
    subplot(1,4,4);
    [avg_jump_distance_p_randperm_1,avg_jump_distance_p_randperm_r] = plot_shuffled([avg_jump_distance_unshuffled avg_jump_distance_test]);
    xlabel("Average jump distance");  ylabel("count");
    title("Average jump distance (shuffled 1000x)");
    saveas(gcf, ['validation_fig/validate_synthetic_circ_shift' name '.png']);
%     ax=axes(fig,'visible','off'); ax.Title.Visible='on'; 
%     title (ax,"Position bins shuffled 1000x"); 
end