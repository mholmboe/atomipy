function [pairList, distList] = cell_list_test(atom, Box_dim,varargin)
% COMPUTETRICLINICCELLLISTDISTANCES  Computes pairwise distances between
% atoms in a triclinic box (with periodic boundary conditions),
% but only for pairs within a specified cutoff. Uses a cell-list approach
% for efficiency and memory savings.
%
% INPUT:
%   atom   = struct array of size N, with fields:
%               atom(i).x, atom(i).y, atom(i).z  (Cartesian coordinates)
%   a, b, c     = unit cell lengths
%   alpha, beta, gamma = angles (in degrees) of the triclinic cell
%   cutoff = distance cutoff for computing distances
%
% OUTPUT:
%   pairList  = Mx2 array of atom indices (i, j) for pairs within cutoff
%   distList  = Mx1 array of distances corresponding to each row of pairList
%
% NOTE:
%   1. We assume the input coordinates (x, y, z) are already in Cartesian
%      coordinates consistent with the triclinic cell. If they are in
%      fractional coordinates, you must first convert them to Cartesian.
%   2. This example stores all pairs i < j within the cutoff. If you
%      have extremely large systems, consider streaming the output or
%      storing only what you need (e.g., adjacency).

if nargin>2
    cutoff=varargin{1};
else
    cutoff=3;
end

if numel(Box_dim)==9
    Cell=Box_dim2Cell(Box_dim);
elseif numel(Box_dim)==6
    Cell=Box_dim;
elseif numel(Box_dim)==3
    Cell=Box_dim;
end

a=Cell(1);
b=Cell(2);
c=Cell(3);

if numel(Cell)==3
    % Convert angles from degrees to radians
    alpha_rad = deg2rad(90);
    beta_rad = deg2rad(90);
    gamma_rad = deg2rad(90);
else
    % Convert angles from degrees to radians
    alpha_rad = deg2rad(Cell(4));
    beta_rad = deg2rad(Cell(5));
    gamma_rad = deg2rad(Cell(6));
end

%% 0. Constants & Setup
% Make sure angles are in radians internally

% Number of atoms
N = size(atom,2);

% Extract coordinates into a convenient Nx3 matrix
coords = [[atom.x]' [atom.y]' [atom.z]'];

%% 1. Construct the triclinic box matrix H
% The standard convention for a triclinic lattice vectors (in Cartesian):
%   H = [ ax, bx, cx
%         0,  by, cy
%         0,   0, cz ]
%
% where:
%   ax = a
%   bx = b*cos(gamma)
%   cx = c*cos(beta)
%   by = b*sin(gamma)
%   cy = c * (cos(alpha) - cos(beta)*cos(gamma)) / sin(gamma)
%   cz = sqrt( c^2 - cx^2 - cy^2 )
%
% NOTE: This is one standard way. Some codes store H differently.

ax = a;
bx = b*cos(gamma_rad);
by = b*sin(gamma_rad);
cx = c*cos(beta_rad);
cy = c*(cos(alpha_rad) - cos(beta_rad)*cos(gamma_rad)) / sin(gamma_rad);
cz = sqrt(c^2 - cx^2 - cy^2);

H = [ax bx cx;
    0  by cy;
    0   0  cz];

% For minimum-image calculations, we need the inverse of H
Hinv = inv(H);

%% 2. Build the cell list
% Heuristic: Use a cell size that is >= cutoff so that no two atoms
% that could be within 'cutoff' end up in cells that are not neighbors.
cellSize = cutoff;  % you can tweak this (e.g. 0.5*cutoff, etc.)

% We need to figure out how many cells along each direction.
% We'll do this by dividing each cell vector length by cellSize.
%   cellVectorLengths = [a, b, c].
% But strictly for triclinic, we can approximate or take the bounding
% box. Below is a simpler bounding-box approach (not the tightest):
boxCart = H * [1; 1; 1];  % one corner to the opposite corner
% bounding box extents in x,y,z:
boundingBoxSize = [max(abs([ax bx cx])), ...
    max(abs([0  by cy])), ...
    max(abs([0   0  cz]))];

nCells = max(floor(boundingBoxSize / cellSize), 1);

% For each atom, we find which cell it belongs to in a fractional sense.
% We'll do a "fractional" approach: fractionalCoord = Hinv * (cartCoord)
fracCoords = (Hinv * coords')';  % Nx3 in fractional
% Ensure that each fractional coordinate is in [0,1) by mod
fracCoords = fracCoords - floor(fracCoords);

% The cell index for each dimension is floor(fractionalCoord * nCells)
% so that it runs from 0 to nCells(d)-1
% cellIndex = floor(fracCoords .* nCells);
% cellIndex(cellIndex == nCells) = nCells(cellIndex == nCells) - 1;
% cellIndex is Nx3. Convert it to a single linear index or store it as Nx3.
cellIndex = floor(fracCoords .* nCells);

% For dimension x
outOfBound = (cellIndex(:,1) >= nCells(1));
cellIndex(outOfBound,1) = nCells(1) - 1;

% For dimension y
outOfBound = (cellIndex(:,2) >= nCells(2));
cellIndex(outOfBound,2) = nCells(2) - 1;

% For dimension z
outOfBound = (cellIndex(:,3) >= nCells(3));
cellIndex(outOfBound,3) = nCells(3) - 1;


% Let's store each (ix, iy, iz) in a single linear index for convenience
% (Alternatively, you could keep a 3D array of lists.)
cellLinIdx = sub2ind(nCells, ...
    cellIndex(:,1)+1, ...
    cellIndex(:,2)+1, ...
    cellIndex(:,3)+1);

% Create a cell array that for each cell linear index,
% holds a list of the atoms in that cell
numCells = prod(nCells);
cellList = cell(numCells,1);

for iAtom = 1:N
    cIdx = cellLinIdx(iAtom);
    cellList{cIdx} = [cellList{cIdx}, iAtom];
end

%% 3. Identify neighbor cells
% For each cell, we want to check itself + neighboring cells.
% Because we might do periodic boundaries, we need to wrap around.
%
% One approach:
%   neighborShifts = [-1, 0, 1] for each dimension -> 27 neighbors
% In triclinic PBC, it’s typical to stay consistent with the fractional
% domain. So we do a 3D shift in the cell grid with wrap-around (mod).

neighborOffsets = -1:1;
neighborIndices = zeros(27,3);
count = 0;
for ix = neighborOffsets
    for iy = neighborOffsets
        for iz = neighborOffsets
            count = count + 1;
            neighborIndices(count,:) = [ix, iy, iz];
        end
    end
end

%% 4. Loop over all cells, then over neighbor cells
pairList = [];
distList = [];

for cID = 1:numCells
    atomListC = cellList{cID};
    if isempty(atomListC), continue; end

    % Convert linear index cID back to 3D subscript:
    [cx, cy, cz] = ind2sub(nCells, cID);

    for nIdx = 1:size(neighborIndices,1)
        dx = neighborIndices(nIdx,1);
        dy = neighborIndices(nIdx,2);
        dz = neighborIndices(nIdx,3);

        nx = mod(cx-1+dx, nCells(1)) + 1;
        ny = mod(cy-1+dy, nCells(2)) + 1;
        nz = mod(cz-1+dz, nCells(3)) + 1;

        neighborCellLin = sub2ind(nCells, nx, ny, nz);
        atomListN = cellList{neighborCellLin};
        if isempty(atomListN), continue; end

        % To avoid duplicating pairs, we can impose
        % that we only handle neighborCell >= currentCell in some sense.
        % Alternatively, we handle all pairs but only keep i<j.
        % The simplest approach is: if neighborCellLin < cID, skip
        if (neighborCellLin < cID)
            continue;
        end

        % Now compute pairwise distances for atomListC vs atomListN
        for i = 1:length(atomListC)
            iAtom = atomListC(i);
            for j = 1:length(atomListN)
                jAtom = atomListN(j);

                % If cID == neighborCellLin, skip j<=i to avoid double counting
                if (neighborCellLin == cID && jAtom <= iAtom)
                    continue;
                end

                % 5. Minimum-image distance in triclinic:
                % coords(iAtom,:) is x_i
                % coords(jAtom,:) is x_j
                % We do:
                %   dVec = x_j - x_i in fractional space, then wrap to [-0.5,0.5)
                diffFrac = (Hinv*(coords(jAtom,:)' - coords(iAtom,:)'))';
                % wrap to the "nearest" image:
                diffFrac = diffFrac - round2dec(diffFrac);

                % convert back to Cartesian
                dVec = (H * diffFrac')';
                d = norm(dVec);

                if d <= cutoff
                    pairList = [pairList; iAtom, jAtom];
                    distList = [distList; d];
                end
            end
        end
    end
end
end
