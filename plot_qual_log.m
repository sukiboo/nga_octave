# this function generates a log-scale plot of the quality values given in the matrix M, 
# where each row represents an algorithm
# the first half of the rows is the average error, the second half is the maximal error
# algorithm names and error type will be specified in the legend
function plot_qual_log(M, step_size=1, desc='')

### set up the figure
fig = figure('Position', [0,0,1200,900]);
colors = {[1.,.2,.1],[.0,.5,1.],[1.,.5,.0],[.0,.8,.4]};
markers = 'odv^';

### construct grid based on step size
[K,N] = size(M);
x_grid = step_size * [1 : floor(N/step_size)];

### adapt data for log-scale
M_log = M(:,x_grid);
M_max = max(M_log(:));
M_min = min(M_log(M_log ~= 0));
a = M_max / M_min;
M_log = log(M_log/M_min) / log(a);

### set up axes
# x-axis
xlim([step_size-1 N+1]);
set(gca, 'xtick', x_grid);
#xlabel('cardinality of reduced basis');
set(gca, 'xgrid', 'off');
# y-axis
ylim([-.05 1.05]);
set(gca, 'ytick', 0 : .1 : 1);
set(gca, 'ygrid', 'on', 'gridlinestyle', '-', 'gridalpha', .5, 'linewidth', 1);
yticks = get(gca, 'ytick'); 
set(gca,'yaxislocation','right');
ylabels = arrayfun(@(y) sprintf('%.2e', M_min + (M_max-M_min)*(a^y-1)/(a-1)), yticks, 'UniformOutput', false); 
set(gca, 'yticklabel', ylabels, 'fontsize', 14);
set(gca, 'box', 'on');
# title
#title(sprintf('%s quality of reduced bases', desc), 'FontWeight', 'Normal');

### plot each process
hold on
for k = 1 : K
	if nnz(M(k,:)) > 0
		plot(x_grid, M_log(k,:), 'color', colors{k}, 'linewidth', 3, ...
			'marker', markers(k), 'markersize', 15, 'markerfacecolor', colors{k});
	endif
endfor
legend({'    OGA', '    NGA', '    EIM', '    POD'});
hold off

### print the graph
set(gca, 'fontname', 'Times', 'fontsize', 28);
set(legend, 'fontname', 'Times', 'fontsize', 32);
set(gca, 'position', [.025 .075 .85 .875])
set(legend, 'position', [.025 .75 .175 .2])
set(get(legend, 'children')(1:2:7), 'markersize', 12);
set(get(legend, 'children')(2:2:8), 'linewidth', 3);
set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0 0 16 9])
name = strftime("%Y-%m-%d_%H-%M-%S", localtime(time()));
print(fig, sprintf('pictures/%s_%s.svg', name, desc), '-dsvg');

endfunction