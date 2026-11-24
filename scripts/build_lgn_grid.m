function [grid, x, y] = build_lgn_grid()
%BUILD_LGN_GRID Generate the 10x11 LGN lattice and cortical 4HC E/I grids.
%   [grid, x, y] = build_lgn_grid()
%
% Outputs:
%   grid - 2 x 110 LGN receptive-field centers (degrees); odd indices are
%          ON cells and even indices are OFF cells.
%   x    - 13689 x 2 evenly spaced E-cell coordinates on a 117 x 117 grid
%          spanning [-0.25, 0.25] x [-0.25, 0.25] deg (4 hypercolumns).
%   y    - 2401 x 2 evenly spaced I-cell coordinates on a 49 x 49 grid
%          spanning [-0.25, 0.25] x [-0.25, 0.25] deg.

grid = zeros(2, 110);
grid(1, :) = kron([-3/2, -7/6, -5/6, -1/2, -1/6, 1/6, 1/2, 5/6, 7/6, 3/2] * 0.25, ones(1, 11));
grid(2, :) = kron(ones(1, 10), [5/3, 4/3, 1, 2/3, 1/3, 0, -1/3, -2/3, -1, -4/3, -5/3] * 0.25);

ne_side = 117;
dx = 0.5 / ne_side;
coords_e = (-0.25 + dx/2) : dx : (0.25 - dx/2);
x = [kron(coords_e, ones(1, ne_side))', kron(ones(1, ne_side), coords_e)'];

ni_side = 49;
dy = 0.5 / ni_side;
coords_i = (-0.25 + dy/2) : dy : (0.25 - dy/2);
y = [kron(coords_i, ones(1, ni_side))', kron(ones(1, ni_side), coords_i)'];
end
