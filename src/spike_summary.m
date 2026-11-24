function summary=spike_summary(spike,Cat,roi_e,roi_i,window_s,preset)
%SPIKE_SUMMARY Per-cell rates and within-cell ISIs in (start, end].
if nargin<6,preset='high_cx';end
cat=Cat(roi_e,:);
if strcmp(preset,'low_lgn')
    groups={roi_e,roi_e(cat(:,1)~=2),roi_e(cat(:,1)==2), ...
        roi_e(cat(:,1)==2 & cat(:,2)==0),roi_e(cat(:,1)==2 & cat(:,2)==1), ...
        roi_e(cat(:,1)==2 & cat(:,2)==2),roi_i};
    summary.labels={'All E','Simple','Complex','Cx 0 LGN','Cx 1 LGN','Cx 2 LGN','I'};
else
groups={roi_e,roi_e(cat(:,1)~=2),roi_e(cat(:,1)==2 & cat(:,2)>=2 & cat(:,2)<=5), ...
    roi_e(cat(:,1)==2 & cat(:,2)>=7 & cat(:,2)<=9),roi_e(cat(:,1)==2 & cat(:,2)>9),roi_i};
summary.labels={'All E','Simple','Cx low','Cx medium','Cx high','I'};
end
summary.window_s=window_s;
summary.cell_count=cellfun(@numel,groups);
summary.mean_hz=nan(1,numel(groups));
summary.spike_count=zeros(1,numel(groups));
summary.isi_count=zeros(1,numel(groups));
summary.median_isi_ms=nan(1,numel(groups));
for g=1:numel(groups)
    intervals=[];
    for id=groups{g}
        t=spike(2:spike(1,id)+1,id);
        t=t(t>window_s(1) & t<=window_s(2));
        summary.spike_count(g)=summary.spike_count(g)+numel(t);
        intervals=[intervals;diff(t)*1000]; %#ok<AGROW>
    end
    if ~isempty(groups{g})
        summary.mean_hz(g)=summary.spike_count(g)/(numel(groups{g})*diff(window_s));
    end
    summary.isi_count(g)=numel(intervals);
    if ~isempty(intervals),summary.median_isi_ms(g)=median(intervals);end
end
end
