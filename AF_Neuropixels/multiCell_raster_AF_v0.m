function multiCell_raster_AF_v0(s,extra_t)

%=== Make sure the inputs are in the right format
if isrow(s), s=s'; end

%=== Extract relevant quantities
n_cells = size(s,1);

%=== Define y coordinates for the cells
y_cells = num2cell([1:n_cells]');

%=== Associate y coordinate to each spike
s_cell_ids = cellfun(@(x,y) y*ones(size(x)),s,y_cells,'UniformOutput',false);

%=== Define x and y coordinates for the raster
x = vertcat(s{:})';
y = vertcat(s_cell_ids{:})'-0.5;

%=== Define line endpoints
x_start = repmat(x, 2, 1);
x_end = repmat(x, 2, 1);
y_start = [y; y+1];
y_end = [y+1; y+1];

%=== Plot vertical lines
plot([x_start(:), x_end(:)]', [y_start(:), y_end(:)]','r-','LineWidth',1); 

%=== Set axis limits and revert axis
xlim([min(x)-extra_t, max(x)+extra_t]);
ylim([.5,n_cells + .5]); 
set(gca, 'YDir','reverse');

end

