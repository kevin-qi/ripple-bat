function [] = plot_p_value_WF(real, shuffled, p_val)
%=== Plot histogram and real value
histogram(shuffled,'FaceColor','k','edgecolor','none');   hold on;    plot(real*[1 1],ylim,'LineWidth',3,'Color','r'); ylabel('Counts');   
text = ['p = ' num2str(p_val)];
hold on;textscatter(max(xlim)*0.5, max(ylim)*0.9,convertCharsToStrings(text),'FontSize',12,'ColorData',[1 0 0]);
end