function plot_run(spike,roi_e,roi_i,cfg,summary,destination)
%PLOT_RUN Show every raw ROI spike; E red, I blue, original anatomical order.
f=figure('Visible','off','Position',[100 100 1000 450],'Color','w');
cleanup=onCleanup(@() close(f));
ax=axes(f); hold(ax,'on');
ids=[roi_e roi_i];
for k=1:numel(ids)
    id=ids(k); t=spike(2:spike(1,id)+1,id);
    if k<=numel(roi_e),color=[0.7 0.15 0.12];else,color=[0.1 0.3 0.75];end
    scatter(ax,t,repmat(k,size(t)),2,color,'filled');
end
xlim(ax,[0 cfg.duration_s]);ylim(ax,[0 numel(ids)+1]);
xlabel(ax,'Time (s)');ylabel(ax,'ROI index (E then I)');
title(ax,sprintf('%s | %s | E %.2f Hz, I %.2f Hz in (%.2f, %.2f] s', ...
    strrep(cfg.preset,'_',' '),cfg.condition,summary.mean_hz(1),summary.mean_hz(end),summary.window_s));
exportgraphics(f,destination,'Resolution',150);
end
