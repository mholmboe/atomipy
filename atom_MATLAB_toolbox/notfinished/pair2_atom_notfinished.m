%% pair2_atom.m
% * This function is adapted from the Matlab MDtoolbox by Yasuhiro 
% * Matsunaga, link to github
% * http://github.com/ymatsunaga/mdtoolbox/
% * Not finished yet
% * Please report problems/bugs to michael.holmboe@umu.se

function [pair,dist,num,rx,ry,rz] = pair2_atom(atom1,atom2,Box_dim,rdist,varargin)

if nargin>4
    rshort=rdist;
    rdist=varargin{1};
end

XYZ_data1=[[atom1.x]' [atom1.y]' [atom1.z]'];
XYZ_data2=[[atom2.x]' [atom2.y]' [atom2.z]'];

if numel(Box_dim)==3
    rx = bsxfun(@minus,XYZ_data1(:,1),XYZ_data2(:,1)');
    rx = rx - Box_dim(1)*round2dec(rx./Box_dim(1));
    ry = bsxfun(@minus,XYZ_data1(:,2),XYZ_data2(:,2)');
    ry = ry - Box_dim(2)*round2dec(ry./Box_dim(2));
    rz = bsxfun(@minus,XYZ_data1(:,3),XYZ_data2(:,3)');
    rz = rz - Box_dim(3)*round2dec(rz./Box_dim(3));
else %if numel(Box_dim)==9 && numel(nonzeros(Box_dim(1,4:end))) > 0 % Untested
    rx = bsxfun(@minus,XYZ_data1(:,1),XYZ_data2(:,1)');
    ry = bsxfun(@minus,XYZ_data1(:,2),XYZ_data2(:,2)');
    rz = bsxfun(@minus,XYZ_data1(:,3),XYZ_data2(:,3)');
    
    rx = rx - xz*round2dec(rz./Box_dim(3));
    ry = ry - yz*round2dec(rz./Box_dim(3));
    rz = rz - Box_dim(3)*round2dec(rz./Box_dim(3));
    
    rx = rx - xy*round2dec(ry./Box_dim(2));
    ry = ry - Box_dim(2)*round2dec(ry./Box_dim(2));
    
    rx = rx - Box_dim(1)*round2dec(rx./Box_dim(1));
end
dist = sqrt(rx.^2 + ry.^2 + rz.^2);
index = find(dist < rdist);
if isempty(index)
    pair = [];
    dist = [];
    num = 0;
    rx = [];
    ry = [];
    rz = [];
else
    [pair2, pair1] = ind2sub(size(dist), index);
    pair = zeros(numel(pair1), 2);
    pair(:,1) = pair1;
    pair(:,2) = pair2;
    dist = dist(index);
    num = numel(dist);
end

if nargin>4
    ind_H=find(strncmpi([atom.type],'H',1));
    pair_ind_H=find(ismember(pair(:,1),ind_H));
    pair_ind_H=sort([pair_ind_H; find(ismember(pair(:,2),ind_H))]);
    H_dist=find(dist<rshort);
    rm_ind=setdiff(pair_ind_H,H_dist);
    pair(rm_ind,:)=[];
    dist(rm_ind)=[];
    num = numel(dist);
end

end
