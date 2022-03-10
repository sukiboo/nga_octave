# this function generates a log-scale plot of the approximation errors given in the matrix M, 
# where each row represents an algorithm
# the first half of the rows is the average error, the second half is the maximal error
# algorithm names and error type will be specified in the legend
function plot_appr_log(M, step_size=1)

### set up the figure
fig = figure('Position', [0,0,1200,900]);
colors = {[1.,.2,.1],[.0,.5,1.],[1.,.5,.0],[.0,.8,.4],[1.,.2,.1],[.0,.5,1.],[1.,.5,.0],[.0,.8,.4]};
markers = 'odv^odv^';
styles = {'-','--'};

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
xlabel('cardinality of reduced basis');
set(gca, 'xgrid', 'off');
# y-axis
ylim([-.05 1.05]);
set(gca, 'ytick', 0 : .1 : 1);
set(gca, 'ygrid', 'on', 'gridlinestyle', '-', 'gridalpha', .5, 'linewidth', 1.1);
yticks = get(gca, 'ytick');
set(gca,'yaxislocation','right');
ylabels = arrayfun(@(y) sprintf('%.2e', M_min + (M_max-M_min)*(a^y-1)/(a-1)), yticks, 'UniformOutput', false); 
set(gca, 'yticklabel', ylabels, 'fontsize', 14);
set(gca, 'box', 'on');
# title
title('approximation error on the training set', 'FontWeight', 'Normal');

### plot each process
hold on
for k = 1 : K
	if nnz(M(k,:)) > 0
		plot(x_grid, M_log(k,:), 'color', colors{k}, styles{ceil(k/4)}, 'linewidth', 2, ...
			'marker', markers(k), 'markersize', 8, 'markerfacecolor', colors{k});
	endif
endfor

### add legend
plot(	1,-1, 'o;  GA ;', 'color', colors{1}, 'markersize', 10, 'markerfacecolor', colors{1}, ...
		1,-1, 'd;  NGA;', 'color', colors{2}, 'markersize', 10, 'markerfacecolor', colors{2}, ...
		1,-1, 'v;  EIM;', 'color', colors{3}, 'markersize', 10, 'markerfacecolor', colors{3}, ...
		1,-1, '^;  POD;', 'color', colors{4}, 'markersize', 10, 'markerfacecolor', colors{4}, ...
		1,-1, '-k;  Average error;', 'linewidth', 2.5, ...
		1,-1, '--k;  Maximal error;', 'linewidth', 2.5);
legend();
hold off

### print the graph
set(gca, 'fontname', 'Times', 'fontsize', 32);
set(legend, 'fontname', 'Times', 'fontsize', 32);
set(gca, 'position', [.025 .125 .85 .8])
set(legend, 'position', [.025 .125 .25 .3])
set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0 0 15 10])
name = strftime("%Y-%m-%d_%H-%M-%S", localtime(time()));
print(fig, sprintf('pictures/%s_appr.eps', name), '-depsc2');

endfunction