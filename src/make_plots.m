function make_plots(time,X,symbols,title,saveFigs,filename)
n = numel(symbols);

% Create subplots
figure
for i = 1:n
    ax(i) = subplot(n,1,i); %#ok<AGROW>
    plot(ax(i),time,X(i,:),'-', 'LineWidth', 1.5);
    grid(ax(i),'on');grid(ax(i),'minor')
    ylabel(ax(i),symbols{i},'Interpreter','latex', 'FontSize', 20)
end
xlabel('Time (s)');
linkaxes(ax,'x')
sgtitle(title);

if saveFigs
    saveas(gcf,filename);
end