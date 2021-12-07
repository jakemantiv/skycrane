function make_plots(opts,time,varargin)

% Extract options information
symbols = opts.symbols;
title = opts.title;
legends = opts.legends;
saveFigs = opts.saveFigs;
filename = opts.filename;
limfrac = opts.limfrac;
linespecs = opts.linespecs;

% Define sizes
n = numel(symbols);
m = numel(varargin);

% Limit y-axis based on last limfrac of data
N = numel(time);
mask = round(N*(1-limfrac),0):N;

% Create subplots
f = figure('Units','normalized','Position',[0,0,.8,.8]);
for i = 1:n
    ylims = [inf,-inf];
    ax(i) = subplot(n,1,i); %#ok<AGROW>
    for j = 1:m
        X = varargin{j};
        plot(ax(i),time,X(i,:),linespecs{i,j}, 'LineWidth', 1.5);
        hold(ax(i),'on');
        ylims = [min([ylims(1),X(i,mask)]), max([ylims(2),X(i,mask)])];
    end
    legend(ax(i),legends,'Location','EastOutside');
    hold(ax(i),'off');
    grid(ax(i),'on');grid(ax(i),'minor')
    ylabel(ax(i),symbols{i},'Interpreter','latex', 'FontSize', 20)
    ylim(ax(i),ylims + [-.1,.1].*diff(ylims));
end
xlabel('Time (s)');
linkaxes(ax,'x')
sgtitle(title);

if saveFigs
    saveas(f,filename);
    exportgraphics(f,[filename, '.png']);
end