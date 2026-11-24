function report=reproduce_historical_onset()
%REPRODUCE_HISTORICAL_ONSET Reanalyse frozen complex-cell histograms.
% This checks an archived analysis, not a fresh stochastic simulation.
root=fileparts(fileparts(mfilename('fullpath')));
s=load(fullfile(root,'data','historical_onset_reference.mat'));
[T50,best_ms]=historical_t50(s.counts,s.legacy_rate_factor);
assert(isequal(T50,s.T50) && isequal(best_ms,s.t50_best));
report.passed=true;
report.cell_phase_comparisons=numel(T50);
report.maximum_error_ms=max(abs(T50-double(s.T50)),[],'all');
report.plot_offset_ms=double(s.plot_offset_ms);
groups={s.Cat(:,2)>9,s.Cat(:,2)>=7 & s.Cat(:,2)<=9,s.Cat(:,2)>=2 & s.Cat(:,2)<=5};
report.group_order={'Cx high','Cx medium','Cx low'};
report.group_count=cellfun(@nnz,groups);
report.median_with_historical_plot_offset_ms=cellfun(@(g) median(best_ms(g)+double(s.plot_offset_ms)),groups);
fprintf('Historical complex onset: %d / %d exact matches.\n',nnz(T50==s.T50),numel(T50));
fprintf('High / medium / low medians, with original +4 ms plot offset: %.1f / %.1f / %.1f ms\n', ...
    report.median_with_historical_plot_offset_ms);
end
