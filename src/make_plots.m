function make_plots(opts,time,varargin)

% Extract options information
symbols = opts.symbols;
title = opts.title;
legends = opts.legends;
saveFigs = opts.saveFigs;
filename = opts.filename;

% Define sizes
n = numel(symbols);
m = numel(varargin);

% Create subplots
figure
for i = 1:n
    ax(i) = subplot(n,1,i); %#ok<AGROW>
    for j = 1:m
        X = varargin{j};
        plot(ax(i),time,X(i,:),'-', 'LineWidth', 1.5);
        hold(ax(i),'on');
    end
    legend(ax(i),legends,'Location','best');
    hold(ax(i),'off');
    grid(ax(i),'on');grid(ax(i),'minor')
    ylabel(ax(i),symbols{i},'Interpreter','latex', 'FontSize', 20)
end
xlabel('Time (s)');
linkaxes(ax,'x')
sgtitle(title);

if saveFigs
    saveas(gcf,filename);
end