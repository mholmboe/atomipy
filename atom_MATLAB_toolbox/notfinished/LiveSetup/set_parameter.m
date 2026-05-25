function set_parameter(oldparam,newparam,infilename,varargin)


if iscell(oldparam)
    oldparam=char(oldparam);
end
   
if iscell(newparam)
    newparam=char(newparam);
end

if nargin>3
   outfilename=varargin{1};
else
   outfilename=infilename; % 'tempfile.txt';
end
    
fid  = fopen(infilename,'r');
f=fread(fid,'*char')';
fclose(fid);

f = strrep(f,oldparam,newparam);

fid  = fopen(outfilename,'w');
fprintf(fid,'%s',f);
fclose(fid);

% if nargin==3
%    copyfile('tempfile.txt',infilename);
%    delete('tempfile.txt');
% end