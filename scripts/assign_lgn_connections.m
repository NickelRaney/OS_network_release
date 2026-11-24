function [LVwire, Cat, LVwire_i, Cat_i] = assign_lgn_connections(x, y, template, model_type)
%ASSIGN_LGN_CONNECTIONS Assign presynaptic LGN afferents to L4 E and I cells.
%   [LVwire, Cat, LVwire_i, Cat_i] = assign_lgn_connections(x, y, template, model_type)
%
% Inputs:
%   x          - Ne x 2 excitatory coordinates (from build_lgn_grid)
%   y          - Ni x 2 inhibitory coordinates (from build_lgn_grid)
%   template   - Simple-cell template struct (from build_lgn_templates)
%   model_type - 'disordered_lgn' (default) or 'low_lgn'
%
% Outputs:
%   LVwire     - Ne x 14 matrix of presynaptic LGN indices (1..110, zero-padded)
%   Cat        - Ne x 3 [class, nLGN, polarity] (class==2 denotes complex)
%   LVwire_i   - Ni x 7 matrix of presynaptic LGN indices for I cells
%   Cat_i      - Ni x 1 presynaptic LGN counts for I cells

if nargin < 4 || isempty(model_type), model_type = 'disordered_lgn'; end
if strcmp(model_type, 'high_cx'), model_type = 'disordered_lgn'; end
assert(ismember(model_type, {'disordered_lgn','low_lgn'}), ...
    'OSNetwork:ModelType', 'model_type must be ''disordered_lgn'' or ''low_lgn''.');

ne = size(x, 1);
LVwire = zeros(ne, 14);
Cat = zeros(ne, 3);

for i = 1:ne
    pt = x(i, :);
    ind3 = 6 - floor(pt(2)*12 + 0.5); % nearest LGN row
    ind4 = floor(pt(1)*12) + 6;       % nearest LGN column
    r = rand(1);
    orient = orientation_domain(pt);

    if r < 0.67
        [wire, polar] = pick_simple(ind3, ind4, template, orient);
        n_lgn = numel(wire);
        LVwire(i, 1:n_lgn) = wire;
        Cat(i, :) = [n_lgn, n_lgn, polar];
    elseif r < 0.70
        [wire, polar] = pick_simple_wide(ind3, ind4, template, orient);
        n_lgn = numel(wire);
        LVwire(i, 1:n_lgn) = wire;
        Cat(i, :) = [n_lgn + 3, n_lgn, polar];
    else
        wire = pick_complex(ind3, ind4, model_type);
        if isempty(wire) || (isscalar(wire) && wire == 0)
            Cat(i, :) = [2, 0, 0];
        else
            n_lgn = numel(wire);
            LVwire(i, 1:n_lgn) = wire;
            Cat(i, :) = [2, n_lgn, 0];
        end
    end
end

ni = size(y, 1);
LVwire_i = zeros(ni, 7);
Cat_i = zeros(ni, 1);
for i = 1:ni
    pt = y(i, :);
    ind3 = 6 - floor(pt(2)*12 + 0.5);
    ind4 = floor(pt(1)*12) + 6;
    wire = pick_inhibitory(ind3, ind4);
    n_lgn = numel(wire);
    LVwire_i(i, 1:n_lgn) = wire;
    Cat_i(i) = n_lgn;
end
end

function orient = orientation_domain(pt)
ind1 = floor(pt*[1; 1]*4 + 0.5);
res1 = pt*[1; 1]*4 + 0.5 - ind1;
b1 = (res1 < 0.05) || (res1 > 0.95);
ind2 = floor(pt*[-1; 1]*4 + 0.5);
res2 = pt*[-1; 1]*4 + 0.5 - ind2;
b2 = (res2 < 0.05) || (res2 > 0.95);

if b1 && b2
    orient = ceil(rand(1)*4);
    return;
end
if mod(ind1, 2) == 0
    if mod(ind2, 2) == 0 % 0 deg
        orient = 1;
        if b1 && rand(1) < 0.5, orient = 2;
        elseif b2 && rand(1) < 0.5, orient = 4; end
    else % 135 deg
        orient = 4;
        if b1 && rand(1) < 0.5, orient = 3;
        elseif b2 && rand(1) < 0.5, orient = 1; end
    end
else
    if mod(ind2, 2) == 0 % 45 deg
        orient = 2;
        if b1 && rand(1) < 0.5, orient = 1;
        elseif b2 && rand(1) < 0.5, orient = 3; end
    else % 90 deg
        orient = 3;
        if b1 && rand(1) < 0.5, orient = 4;
        elseif b2 && rand(1) < 0.5, orient = 2; end
    end
end
end

function [wire, polar] = pick_simple(ind3, ind4, template, orient)
polar = ceil(rand(1)*2);
r = rand(1);
if r < 0.25
    n = 2; % 5 LGN cells (1:3 ratio vs 6 LGN cells)
else
    n = 3; % 6 LGN cells
end
pool = template.P{n, mod(orient-1, 2) + 1};
idx = pool(ceil(rand(1)*numel(pool)));
raw = template.temp{n, orient, polar}(idx, :);
wire = raw(1:2:end) - 3 + ind3 + (raw(2:2:end) - 4 + ind4)*11;
end

function [wire, polar] = pick_simple_wide(ind3, ind4, template, orient)
while true
    polar = ceil(rand(1)*2);
    n = ceil(rand(1)*3);
    pool = template.P_w{n, mod(orient-1, 2) + 1};
    idx = pool(ceil(rand(1)*numel(pool)));
    raw = template.temp_w{n, orient, polar}(idx, :);
    wire = raw(1:2:end) - 3 + ind3 + (raw(2:2:end) - 4 + ind4)*11;
    if max(wire) <= 110 && min(wire) > 0
        break;
    end
end
end

function wire = pick_complex(ind3, ind4, model_type)
if strcmp(model_type, 'low_lgn')
    n = floor(rand(1)*3); % 0, 1, or 2 with equal probability 1/3
    if n == 0
        wire = 0;
        return;
    end
    a = zeros(2, 2);
    a(randperm(4, n)) = 1;
    [row, col] = find(a == 1);
    row = row + floor(rand(1)*4);
    col = col + floor(rand(1)*4);
    wire = (row - 3 + ind3 + (col - 4 + ind4)*11)';
    wire = wire(wire >= 1 & wire <= 110);
    return;
end

% Disordered-LGN complex distribution:
% 40% high (10-11), 40% mid (7-8), 20% low (2-5)
r = rand(1);
if r < 0.40
    n = 10 + (rand(1) >= 0.5);
    a = zeros(4, 4); a(randperm(16, n)) = 1;
    [row, col] = find(a == 1);
    wire = (row - 3 + ind3 + (col - 4 + ind4)*11)';
elseif r < 0.80
    n = 7 + (rand(1) >= 0.5);
    a = zeros(4, 4); a(randperm(16, n)) = 1;
    [row, col] = find(a == 1);
    wire = (row - 3 + ind3 + (col - 4 + ind4)*11)';
else
    n = 2 + floor(rand(1)*4); % 2, 3, 4, or 5 uniformly
    a = zeros(4, 4); a(randperm(16, n)) = 1;
    [row, col] = find(a == 1);
    wire = (row - 3 + ind3 + (col - 4 + ind4)*11)';
end
wire = wire(wire >= 1 & wire <= 110);
end

function wire = pick_inhibitory(ind3, ind4)
% Assign 5, 6, or 7 LGN cells uniformly in a 3x3 subgrid within 5x5
n = 5 + floor(rand(1)*3);
while true
    a = zeros(3, 3);
    a(randperm(9, n)) = 1;
    if any(a(:, 2)) && any(a(2, :))
        [row, col] = find(a == 1);
        row = row + floor(rand(1)*2);
        col = col + floor(rand(1)*2);
        wire = (row - 3 + ind3 + (col - 4 + ind4)*11)';
        if max(wire) <= 110 && min(wire) >= 1
            break;
        end
    end
end
end
