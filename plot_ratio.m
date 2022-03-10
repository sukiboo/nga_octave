# this function generates a plot of the values given in the matrix M, 
# where each row represents a process
# process_name and process_time will be specified in the legend
# p and N_tr will be specified in the title
function plot_ratio(ratio_avg, ratio_max, step_size=1)

### set up the figure
fig = figure('Position', [0,0,1200,600]);
colors = {[1.,.2,.1],[.0,.5,1.]};
markers = 'ds';

### title
%title('Ratio of approximation errors of OGA and NGA', 'fontsize', 14);
set(gca, 'box', 'on');

### construct grid based on step size
N = size(ratio_avg,2);
x_grid = step_size * [1 : floor(N/step_size)];

### set up axes
grid minor
# x-axis
xlim([0 N+1]);
set(gca, 'xtick', x_grid);
set(gca, 'xgrid', 'off', 'xminorgrid', 'off');
%xlabel('Cardinality of reduced basis');
# y-axis
ylim([.8 1.1]);
set(gca, 'ytick', .8 : .05 : 1.1);
%ylabel('OGA/NGA approximation error ratio');
set(gca, 'ygrid', 'on', 'gridlinestyle', '--', 'gridalpha', .5, 'minorgridlinestyle', '--', 'minorgridalpha', .5);
yticks = get(gca, 'ytick'); 
set(gca,'yaxislocation','right');

### plot each process
hold on
plot([1:N], ratio_avg(2,:), '-dr', 'linewidth', 2, 'markersize', 9, 'markerfacecolor', 'red');
plot([1:N], ratio_max(2,:), '-sb', 'linewidth', 2, 'markersize', 6, 'markerfacecolor', 'blue');
patch([1:N,N:-1:1], [ratio_avg(1,:), fliplr(ratio_avg(3,:))], 'red', 'FaceAlpha', .3, 'EdgeAlpha', .0)
patch([1:N,N:-1:1], [ratio_max(1,:), fliplr(ratio_max(3,:))], 'blue', 'FaceAlpha', .2, 'EdgeAlpha', .0)
hold off

### add legend
%legend({'   Ratio of average errors', '   Ratio of maximal errors'}, 'location', 'north');
%set(legend, 'fontsize', 20, 'orientation', 'horizontal');

### print the graph
set(gca, 'fontname', 'Times', 'fontsize', 28);
%set(legend, 'fontname', 'Times', 'fontsize', 32, 'orientation', 'horizontal', 'box', 'off');
%set(get(legend, 'children')(1), 'markersize', 12);
%set(get(legend, 'children')(3), 'markersize', 16);
%set(get(legend, 'children')(2:2:4), 'linewidth', 3);
set(gca, 'position', [.005 .075 .925 .875])
%set(legend, 'position', [.006 .85 .923 .123])
set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0 0 16 8])
name = strftime("%Y-%m-%d_%H-%M-%S", localtime(time()));
#print(fig, sprintf('pictures/%s_ratio.eps', name), '-depsc2');
print(fig, sprintf('pictures/%s_ratio.svg', name), '-dsvg');

endfunction