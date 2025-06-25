function [p_val] = p_value_WF(real, shuffled)
%% Calculate p value and plot distribution and real value for a shuffled dataset

%=== Number of repetitions
non_nan = nnz(~isnan(shuffled));

%=== Calculate p_value
if isnan(real)
    p_val_r =  NaN;
    p_val_l =  NaN;
    p_val= NaN;
else
    p_val_r =  nnz(shuffled>real)/non_nan;
    p_val_l =  nnz(shuffled<real)/non_nan;
    p_val = 1-(max(p_val_r,p_val_l));
end
end