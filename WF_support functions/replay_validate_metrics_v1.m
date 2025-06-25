% replay_validate_metrics_v1.m
% this code generates synthetic "replay" and then validate different metrics for
% evaluating replay using the generated data. 

% update from v0:
% - aims to include "replay score" as a metric 

%% synthetic dataset generation - generate a synthetic position and then apply a gaussian around that max posibility
test_time_vector = [1:38];
num_positions = 38;
test_position = test_time_vector*num_positions/numel(test_time_vector);
test_prob = zeros(numel(test_time_vector), numel(test_position));
for i = 1:(numel(test_time_vector))
    test_prob(i,round(test_position(i))) = 1;
end
test_prob_gau = conv2(test_prob,normpdf([-10:1:10],0,2),"same");
figure; imagesc (test_prob_gau'); title("Synthetic data before shuffling"); xlabel("Time bins"); ylabel("Position bins");
%% evaluate how the value of different replay metrics changes with each iteration of shuffling 
% this can probably be improved by averaging each data point over many tries 
shuffle_times = 50; % this determines how many shuffles 
test_prob_gau_shuffled = test_prob_gau; 
weighted_corr_test = [];
max_jump_distance_test = [];
avg_jump_distance_test = [];
num_peaks_test = [];
replay_score_test = [];
for i = 1:shuffle_times 
    first_col = randi([1,38]);
    second_col = randi([1,38]);
    placehold_col = test_prob_gau_shuffled(first_col,:);
    test_prob_gau_shuffled(first_col,:) = test_prob_gau_shuffled(second_col,:);
    test_prob_gau_shuffled(second_col,:) = placehold_col;
    [weighted_corr_test(i),max_jump_distance_test(i), avg_jump_distance_test(i),num_peaks_test(i),replay_score_test(i)] = evaluate_candidate_event_v1(test_prob_gau_shuffled);
end 
figure; 
subplot(2,3,1); imagesc(test_prob_gau_shuffled'); xlabel("time bins"); ylabel("position");
subplot(2,3,2); plot(weighted_corr_test); xlabel("number of shuffles"); ylabel("weighted correlation");
subplot(2,3,3); plot(max_jump_distance_test);  xlabel("number of shuffles"); ylabel("max jump distance");
subplot(2,3,4); plot(avg_jump_distance_test);  xlabel("number of shuffles"); ylabel("avg jump distance");
subplot(2,3,5); plot(replay_score_test);  xlabel("number of shuffles"); ylabel("replay score");

%% generate high noise and low noise data and compare metrics 
% high noise: random, simulates "no replay"
% low noise: potential replay with some noise
shuffle_times_noise_high = 50;
shuffle_times_noise_low = 7;
total_shuffles = 1000; % this is the amount of repeat times 
weighted_corr_test_low_noise = [];
max_jump_distance_low_noise = [];
avg_jump_distance_low_noise = [];
replay_score_low_noise = [];

weighted_corr_test_high_noise = [];
max_jump_distance_high_noise = [];
avg_jump_distance_high_noise = [];
replay_score_high_noise = [];


for i = 1:total_shuffles
    test_prob_gau_shuffled = test_prob_gau;
    for ii = 1:shuffle_times_noise_high 
        first_col = randi([1,38]);
        second_col = randi([1,38]);
        placehold_col = test_prob_gau_shuffled(first_col,:);
        test_prob_gau_shuffled(first_col,:) = test_prob_gau_shuffled(second_col,:);
        test_prob_gau_shuffled(second_col,:) = placehold_col;
    end 
    [weighted_corr_test_high_noise(i),max_jump_distance_test_high_noise(i), avg_jump_distance_test_high_noise(i),num_peaks_test_high_noise(i), replay_score_high_noise(i)] = evaluate_candidate_event_v1(test_prob_gau_shuffled);
end
test_prob_gau_high_noise = test_prob_gau_shuffled;

for i = 1:total_shuffles
    test_prob_gau_shuffled = test_prob_gau;
    for ii = 1:shuffle_times_noise_low 
        first_col = randi([1,38]);
        second_col = randi([1,38]);
        placehold_col = test_prob_gau_shuffled(first_col,:);
        test_prob_gau_shuffled(first_col,:) = test_prob_gau_shuffled(second_col,:);
        test_prob_gau_shuffled(second_col,:) = placehold_col;
    end 
    [weighted_corr_test_low_noise(i),max_jump_distance_test_low_noise(i), avg_jump_distance_test_low_noise(i),num_peaks_test_low_noise(i), replay_score_low_noise(i)] = evaluate_candidate_event_v1(test_prob_gau_shuffled);
end

% plotting some examples of high noise and low noise results 
fig = figure;
subplot(1,2,1); imagesc(test_prob_gau_high_noise'); xlabel("time bins"); ylabel("position"); title("Example of synthetic data with high noise");
axis xy;
subplot(1,2,2); imagesc(test_prob_gau_shuffled'); xlabel("time bins"); ylabel("position"); title("Example of synthetic data with low noise");
axis xy;
colormap(hot);
% plotting the histogram distribution of shuffled results
fig = figure;
subplot(1,4,1);hist([max_jump_distance_test_high_noise' max_jump_distance_test_low_noise']); xlabel("Maximum jump distance"); ylabel("count");
legend('high noise','low noise','Location','northwest')
subplot(1,4,2);hist([weighted_corr_test_high_noise' weighted_corr_test_low_noise']); xlabel("Weighted correlation"); ylabel("count");
subplot(1,4,3); hist([avg_jump_distance_test_high_noise' avg_jump_distance_test_low_noise']);  xlabel("Average jump distance");  ylabel("count");
subplot(1,4,4); hist([replay_score_high_noise' replay_score_low_noise']);  xlabel("Replay score");  ylabel("count");

ax=axes(fig,'visible','off'); ax.Title.Visible='on'; 
title (ax,"generated data shuffled 1000x (bin shuffle)"); 

