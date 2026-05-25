%% neigh_atom.m
% * This function checks which neighbors each atom has and ouputs their info
%
%% Version
% 3.00
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% * atom = neigh_atom(atom,Box_dim,rmax)
% * atom = neigh_atom(atom,Box_dim,rmax,101)

function atom = new_neigh_atom(atom,Box_dim,rmax,varargin)

if length(rmax) == 1
    % Simple way to set the rmax for each atom
    rmax = rmax(1)*ones(size(atom,2),1);
    % Special thing for H's
    rmax(strncmpi([atom.type],'H',1))=1.25;
end

% if we want to exclude some particles from the distance matrix calculation
% This is currenlty not implemented here in this new function
% if nargin==4
%     skip_ind=varargin{1};
% else
%     skip_ind=[];
% end

if numel(Box_dim)==1
    Box_dim(1)=Box_dim(1);
    Box_dim(2)=Box_dim(1);
    Box_dim(3)=Box_dim(1);
    xy=0;xz=0;yz=0;
elseif numel(Box_dim)==3
    xy=0;xz=0;yz=0;
else
    xy=Box_dim(6);xz=Box_dim(8);yz=Box_dim(9);
end

XYZ_data=[[atom.x]' [atom.y]' [atom.z]'];

%% calculate the distances of all pairs
dX = bsxfun(@minus,XYZ_data(:,1),XYZ_data(:,1)');
dY = bsxfun(@minus,XYZ_data(:,2),XYZ_data(:,2)');
dZ = bsxfun(@minus,XYZ_data(:,3),XYZ_data(:,3)');

if numel(Box_dim)==3
    dX = dX - Box_dim(1)*round2dec(dX./Box_dim(1));
    dY = dY - Box_dim(2)*round2dec(dY./Box_dim(2));
    dZ = dZ - Box_dim(3)*round2dec(dZ./Box_dim(3));
elseif numel(Box_dim)==9
    dX = dX - Box_dim(1)*round2dec(dX./Box_dim(1));
    
    dY = dY - Box_dim(2)*round2dec(dY./Box_dim(2));
    dX = dX - xy*round2dec(dY./Box_dim(2));
    
    dZ = dZ - Box_dim(3)*round2dec(dZ./Box_dim(3));
    dX = dX - xz*round2dec(dZ./Box_dim(3));
    dY = dY - yz*round2dec(dZ./Box_dim(3));
end

% index = find(triu(dist < rcut, 1));
% [pair1, pair2] = ind2sub(size(dist), index);
% pair = [pair1 pair2];
% dist = dist(index);
dist = sqrt(dX.^2 + dY.^2 + dZ.^2);
indmatrix = bsxfun(@lt,dist,rmax');
dist=bsxfun(@times,dist,indmatrix);
index = find(triu(indmatrix, 1));
[type1, type2] = ind2sub(size(dist), index);
bond_index = [type1 type2 dist(index)];
bond_index = sortrows(bond_index);

%% non-pbc way of doing it
% XYZ_data=[[atom.x]' [atom.y]' [atom.z]'];
% XYZ_inner_product = XYZ_data * XYZ_data';
% dist_matrix = sqrt(bsxfun(@plus,diag(XYZ_inner_product),diag(XYZ_inner_product)')-2*XYZ_inner_product);

%% How to speed this up? Using some cool struct fieldname assigment?
for i=1:size(XYZ_data,1)

    in=find(dist(i,:)>0);
    atom(i).neigh.index = in;
    atom(i).neigh.type = deal([atom([atom(i).neigh.index]).type]);
    atom(i).neigh.dist = [dist(in,1)];
    atom(i).neigh.coords = [XYZ_data(in,1) XYZ_data(in,2) XYZ_data(in,3)];
    %     atom(i).neigh.r_vec = [rx(in) ry(in) rz(in)]; % This does not work in this new function
    if mod(i,100)== 0
        i
    end
end

end