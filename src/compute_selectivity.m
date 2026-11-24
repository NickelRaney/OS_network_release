function sel = compute_selectivity(fr_fwd, fr_rev, fr_ortho, Cat, roi_e, preset)
%COMPUTE_SELECTIVITY Directional selectivity (DS) and orthogonal/preferred (OP) ratio.
%   sel = compute_selectivity(fr_fwd, fr_rev, fr_ortho, Cat, roi_e, preset)
%
% Inputs:
%   fr_fwd   - Firing rates (Hz) for optimal orientation, forward drift (e.g. [1 0])
%   fr_rev   - Firing rates (Hz) for optimal orientation, opposite drift (e.g. [-1 0])
%   fr_ortho - Firing rates (Hz) for orthogonal orientation (e.g. [0 1])
%   Cat      - Excitatory cell category matrix (13689 x 3)
%   roi_e    - Global excitatory indices in the central ROI
%   preset   - 'disordered_lgn' (default; alias 'high_cx') or 'low_lgn'
%
% Definitions:
%   DS = max(fr_fwd, fr_rev) / min(fr_fwd, fr_rev)
%   OP = fr_ortho / max(fr_fwd, fr_rev)

if nargin < 6 || isempty(preset), preset = 'disordered_lgn'; end
if strcmp(preset, 'high_cx'), preset = 'disordered_lgn'; end

n_e = size(Cat, 1);
fr_fwd = reshape(fr_fwd, 1, []);
fr_rev = reshape(fr_rev, 1, []);
fr_ortho = reshape(fr_ortho, 1, []);
if numel(fr_fwd) > n_e, fr_fwd = fr_fwd(1:n_e); end
if numel(fr_rev) > n_e, fr_rev = fr_rev(1:n_e); end
if numel(fr_ortho) > n_e, fr_ortho = fr_ortho(1:n_e); end

fr_pref = max(fr_fwd, fr_rev);
fr_opp = min(fr_fwd, fr_rev);

ds_raw = fr_pref ./ fr_opp;
ds_capped = ds_raw;
ds_capped(isnan(ds_capped) | ds_capped > 8) = 9;

op_ratio = fr_ortho ./ fr_pref;
op_ratio(fr_pref == 0) = nan;

sel.fr_pref = fr_pref;
sel.fr_opp = fr_opp;
sel.fr_ortho = fr_ortho;
sel.preferred_is_fwd = fr_fwd >= fr_rev;
sel.ds_raw = ds_raw;
sel.ds_capped = ds_capped;
sel.op_ratio = op_ratio;

cat_roi = Cat(roi_e, :);
if strcmp(preset, 'low_lgn')
    groups = {roi_e, roi_e(cat_roi(:,1)~=2), roi_e(cat_roi(:,1)==2), ...
        roi_e(cat_roi(:,1)==2 & cat_roi(:,2)==0), ...
        roi_e(cat_roi(:,1)==2 & cat_roi(:,2)==1), ...
        roi_e(cat_roi(:,1)==2 & cat_roi(:,2)==2)};
    sel.group_labels = {'All E','Simple','Complex','Cx 0 LGN','Cx 1 LGN','Cx 2 LGN'};
else
    groups = {roi_e, roi_e(cat_roi(:,1)~=2), roi_e(cat_roi(:,1)==2), ...
        roi_e(cat_roi(:,1)==2 & cat_roi(:,2)>=2 & cat_roi(:,2)<=5), ...
        roi_e(cat_roi(:,1)==2 & cat_roi(:,2)>=7 & cat_roi(:,2)<=9), ...
        roi_e(cat_roi(:,1)==2 & cat_roi(:,2)>9)};
    sel.group_labels = {'All E','Simple','Complex','Cx low','Cx medium','Cx high'};
end

n_g = numel(groups);
sel.group_count = cellfun(@numel, groups);
sel.median_ds = nan(1, n_g);
sel.median_op = nan(1, n_g);
for g = 1:n_g
    idx = groups{g};
    if ~isempty(idx)
        sel.median_ds(g) = median(ds_capped(idx), 'omitnan');
        sel.median_op(g) = median(op_ratio(idx), 'omitnan');
    end
end
end
