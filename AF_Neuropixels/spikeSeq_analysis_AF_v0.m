function [rankcorr,p_value2,shuffled_rankcorr] = spikeSeq_analysis_AF_v0(s,n_shuffle)

% === INPUTS:
% s is a cell array with the spike times of each unit of the sequence
% When the array is empty, it means there where no spikes from that unit
% n_shuffle is the number of shuffles

%=== Basic preprocessing
if isrow(s), s=s'; end      % Make sure the inputs are in the right format
N = size(s,1);              % Total number of cells
tgt_seq = [1:N]';           % Target sequence of activation

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

end

