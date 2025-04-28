function [rankcorr,p_value2,shuffled_rankcorr,a,cell_pairs,pw_timing] = spikeSeq_analysis_AF_v1(s,n_shuffle)

% === INPUTS:
% s is a cell array with the spike times of each unit of the sequence
% When the array is empty, it means there where no spikes from that unit
% n_shuffle is the number of shuffles

%=== Basic preprocessing
if isrow(s), s=s'; end      % Make sure the inputs are in the right format
N = size(s,1);              % Total number of cells
tgt_seq = [1:N]';           % Target sequence of activation
s_original = s;             % Keep a copy of the original sequence

%=== Left and right extremes
tL = min(vertcat(s{:}));
tR = max(vertcat(s{:}));

%=== Keep the first spike in the sequence and set to Nan the cells that did not fire
s(cellfun(@isempty,s))={0};
s1 = cellfun(@(x) sort(x(1)),s);
s1(s1==0)=NaN;

%=== Calculate the rank of the cells in the sequence
cdt_seq = tiedrank(s1);

%=== Calculate the rank-correlation
[rankcorr,p_value1] = corr(s1,tgt_seq,'Type','Kendall','Rows','Complete');

%=== Shuffling procedure
shuffled_rankcorr = NaN(n_shuffle,1);
for i=1:n_shuffle
    shuffled_rankcorr(i) = corr(s1(randperm(N)),tgt_seq,'Type','Kendall','Rows','Complete');
end

%=== P value calculation
if rankcorr>0,    p_value2 = sum(shuffled_rankcorr>rankcorr)/n_shuffle;
else,             p_value2 = sum(shuffled_rankcorr<rankcorr)/n_shuffle;
end

%=== Warning if NaN Correlation
if isnan(rankcorr)
    warning('NaN Correlation Value');
end

%=== Calculate Replay Speed in terms of cells/second
x = s1(~isnan(s1));
y = tgt_seq(~isnan(s1));
p = polyfit(x, y, 1);
a = p(1);   b = p(2);
y_fit = a*x+b;

%=== Calculate pairwise intervals between first spikes
cell_pairs = nchoosek(1:N,2);
pw_timing = cell_pairs(:,1)*0;
for i=1:size(cell_pairs,1)
    pw_timing(i) = s1(cell_pairs(i,2))-s1(cell_pairs(i,1));
end

cell_pairs = {cell_pairs};
pw_timing = {pw_timing};

% %=== Plotting for debugging
% figure('units','normalized','outerposition',[.2 .3 .07 .3]);
% for nc= 1:N
%     plot(s_original{nc}, nc*ones(size(s_original{nc})), 'r|','MarkerSize', 5);  hold on;             % Raster for each cell
% end
% ylim([0 N+1]); xticks([tL, tR]); xlim([tL, tR]); yticks([]); xticklabels({[],[num2str((tR-tL)*1000,3), ' ms']});
% plot(x,y,'b.',x,y_fit,'k');

end

