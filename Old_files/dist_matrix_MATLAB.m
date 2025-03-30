%% dist_matrix_atom.m
% * This function calculates the distance matrix from the atom struct, or
% the distances between atom1 and atom2
% * Note that there is also and cell_list_dist_matrix_atom function that
% might be faster for large orthogonal systems...
%
%% Version
% 3.00
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% # dist_matrix = dist_matrix_atom(atom1,Box_dim) % Basic input arguments
% # dist_matrix = dist_matrix_atom(atom1,atom2,Box_dim) % Calculates the distance matrix between sites in atom1 and in atom2
%
function [distances,dx,dy,dz] = dist_matrix_MATLAB(atoms,varargin) % ,atom2,Box_dim); % or % ,Box_dim);

format compact;

nAtoms1=size(atoms,2);

XYZ1=single([[atoms.x]' [atoms.y]' [atoms.z]']);
XYZ2=single([[atoms.x]' [atoms.y]' [atoms.z]']);

lx=Box_dim(1);
ly=Box_dim(2);
lz=Box_dim(3);
if length(Box_dim) == 6
    lx=Box_dim(4)-Box_dim(1);
    ly=Box_dim(5)-Box_dim(2);
    lz=Box_dim(6)-Box_dim(3);
    xy=0;
    xz=0;
    yz=0;
elseif length(Box_dim) == 9
    xy=Box_dim(6);
    xz=Box_dim(8);
    yz=Box_dim(9);
else
    xy=0;
    xz=0;
    yz=0;
end

distances = single(zeros(nAtoms1,nAtoms1)); % use of single instead of double
dx = distances;
dy = distances;
dz = distances;
i=1;
if size(Box_dim,2)==3
    while i<size(XYZ1,1)+1
        %Calculate Distance Components
        rx = XYZ1(i,1) - XYZ2(:,1);
        x_gt_ind=find(rx > lx/2); x_lt_ind=find(rx < - lx/2);
        rx(x_gt_ind) = rx(x_gt_ind) - lx;
        rx(x_lt_ind) = rx(x_lt_ind) + lx;

        ry = XYZ1(i,2) - XYZ2(:,2);
        y_gt_ind=find(ry > ly/2); y_lt_ind=find(ry < - ly/2);
        ry(y_gt_ind) = ry(y_gt_ind) - ly;
        ry(y_lt_ind) = ry(y_lt_ind) + ly;

        rz = XYZ1(i,3) - XYZ2(:,3);
        z_gt_ind=find(rz > lz/2); z_lt_ind=find(rz < - lz/2);
        rz(z_gt_ind) = rz(z_gt_ind) - lz;
        rz(z_lt_ind) = rz(z_lt_ind) + lz;

        r = sqrt( rx(:,1).^2 + ry(:,1).^2 + rz(:,1).^2 ); % distance calc.
        distances(:,i)=r;
        dx(:,i)=-rx;
        dy(:,i)=-ry;
        dz(:,i)=-rz;

        i=i+1;
    end
else % if cell is triclinic, and this part is not actually tested yet...
    while i<size(XYZ1,1)+1
        % Calculate Distance Components %%%%%%%%%%%%%%%%%%%%
        rx = XYZ1(i,1) - XYZ2(:,1);
        ry = XYZ1(i,2) - XYZ2(:,2);
        rz = XYZ1(i,3) - XYZ2(:,3);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        z_gt_ind=find(rz > lz/2); z_lt_ind=find(rz < - lz/2);
        rz(z_gt_ind) = rz(z_gt_ind) - lz;
        rz(z_lt_ind) = rz(z_lt_ind) + lz;
        rx(z_gt_ind) = rx(z_gt_ind) - xz;
        rx(z_lt_ind) = rx(z_lt_ind) + xz;
        ry(z_gt_ind) = ry(z_gt_ind) - yz;
        ry(z_lt_ind) = ry(z_lt_ind) + yz;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        y_gt_ind=find(ry > ly/2); y_lt_ind=find(ry < - ly/2);
        ry(y_gt_ind) = ry(y_gt_ind) - ly;
        ry(y_lt_ind) = ry(y_lt_ind) + ly;
        rx(y_gt_ind) = rx(y_gt_ind) - xy;
        rx(y_lt_ind) = rx(y_lt_ind) + xy;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        x_gt_ind=find(rx > lx/2); x_lt_ind=find(rx < - lx/2);
        rx(x_gt_ind) = rx(x_gt_ind) - lx;
        rx(x_lt_ind) = rx(x_lt_ind) + lx;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        r = sqrt( rx(:,1).^2 + ry(:,1).^2 + rz(:,1).^2 ); % distance calc.
        distances(:,i)=r;
        dx(:,i)=-rx;
        dy(:,i)=-ry;
        dz(:,i)=-rz;
        i=i+1;
    end
end

% New transposed output
distances=distances';
assignin('caller','analyzed_Box_dim',Box_dim);
assignin('caller','X_dist',(dx)');
assignin('caller','Y_dist',(dy)');
assignin('caller','Z_dist',(dz)');

end

