function [p_val_l,p_val_r] = plot_shuffled(data)
%% Calculate p value and plot distribution and real value for a shuffled dataset
% data(1,:) is the real datum
% data(2:end,:) contains the shuffled data

%=== Make sure data is a column vector
if isrow(data); data = data';   end

%=== Number of repetitions
non_nan = nnz(~isnan(data(2:end)));

%=== Calculate p_value
if isnan(data(1))
    p_val_l = NaN;
    p_val_r = NaN;
else
    p_val_r =  nnz(data(2:end)>data(1))/non_nan;
    p_val_l =  nnz(data(2:end)<data(1))/non_nan;
end

%=== Plot histogram and real value
histogram(data(2:end),'FaceColor','k','edgecolor','none');   hold on;    plot(data(1)*[1 1],ylim,'LineWidth',3,'Color','r'); ylabel('Counts');   
% title(['p = [', num2str(p_val_l,3) '    ' num2str(p_val_r,3) ']']);
text = ['p = [' num2str(p_val_l,3) '    ' num2str(p_val_r,3) ']'];
hold on;textscatter(max(xlim)*0.5, max(ylim)*0.9,convertCharsToStrings(text),'FontSize',12,'ColorData',[1 0 0]);

end