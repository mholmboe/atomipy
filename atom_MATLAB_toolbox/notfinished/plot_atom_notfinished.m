function plot_atom(plotstr,XYZ_labels,XYZ_data,xylimit,zlimit,Resolution,alfa)


%%
% Sets plot limits for the data
x_lo = 0; x_hi = xylimit;
y_lo = x_lo; y_hi = x_hi;
z_lo = 0; z_hi = zlimit;
%limits = [x_lo x_hi y_lo y_hi z_lo z_hi];

% Set limits of plot
xypadding = 5;
zpadding = 5;
xlo = x_lo-xypadding; xhi = x_hi+xypadding;
ylo = xlo; yhi = xhi;
zlo = z_lo-xypadding; zhi = z_hi + zpadding;
limits = [xlo xhi ylo yhi zlo zhi];
% Set plot size and colors

% Set plot size and colors
FigColor = [1 1 1];
Part_size = [27, 27, 27, 27, 24, 24, 24, 24, 24, 24, 15, 24, 15, 21, 15, 24, 30, 18, 30]; %(0.5-gray/5)*
% RGB color codes % Si     % Al    % Alt     % Mgo      % O       % Oh      % Omg      % Ohmg
Part_color = {[1,1,0];[1,0.5,0.5];[1,0.5,0.5];[0,.9,.9];[1,0,0];[1,.4,.4];[1,.2,.2];[1,0,.2];
    [1,0,.1];[1,0,0];[.9,.9,.9];[1,0,0];[.9,.9,.9];[0,1,0];[.2,.2,.2];[0,0,1];[0,.5,.5];[.2,1,0];[0,1,1]};
     %Oalt     %Odsub      % H     % Ow     % Hw      % K       %Li      % Na      %Cs       %Ca       %Cl


if size(XYZ_data,2)==4;
    XYZ_data(:,1) = [];
end
%1  2 (3) : Al (Mg)
%2  1     : Si (Al)
%3  1     : Si (Al)
%4  4 (6) : O  (Omg)
%5  4 (6) : O  (Omg)
%6  5 (7) : Oh (Ohmg)
%7  4 (6) : O  (Omg)
%8  4 (6) : O  (Omg)
%9  4 (6) : O  (Omg)
%10 10    : H
%countstr = {'Al','Alt','Si','Mg', 'O', 'Oh', 'Omg', 'Ohmg','Odsub','Oalt', 'H', 'Ow', 'Hw', 'Na', 'Ca','Cl'};
countstr = {'Si','Al','Alt','Mgo','O','Oh','Omg','Ohmg','Oalt','Odsub','H','Ow','Hw','K','Li','Na','Cs','Ca','Cl'}; % Octahedral and tetrahedral substitusionsplotstr = {'Alt', 'Mgo','Omg','Ohmg','Odsub','Oalt'};
% plotstr = {'Al','Alt','Si','Mg', 'O', 'Oh', 'Omg', 'Ohmg','Odsub','Oalt', 'H', 'Ow', 'Hw', 'Na'};

count = zeros(size(countstr,2));
XYZ_labels=strtrim(XYZ_labels);

Init_coord = cell(size(XYZ_data,1),size(countstr,2));
for i=1:size(countstr,2);
    if find(strcmp(countstr(i),XYZ_labels)) > 0; %strmatch(countstr(i), textdata,'exact') > 0;
        XYZ = Atom_molID_func(XYZ_labels,XYZ_data,countstr(i));
        Init_coord(:,i) = {XYZ(:,2:4)};
        count(i) = size(XYZ,1);
    end
end
assignin('base','Element_count',count);
assignin('base','Element_countstr',countstr);

hold on;
rotate3d on;
camlight(220,210,'infinite');
set(gca,'PlotBoxAspectRatio',[(xhi-xlo)/(zhi-zlo) (yhi-ylo)/(zhi-zlo) (zhi-zlo)/(zhi-zlo)],'FontSize',21);
set(gcf,'Color',[1,1,1]);

xlabel('X'); ylabel('Y'); zlabel('Z')
axis([xlo xhi ylo yhi zlo zhi]);
view([0,90]);
%Label_coord_func(xhi-5,0,0,15);

for i=1:size(plotstr,2);
    for j=1:size(Init_coord,2);
        if find(strcmp(plotstr(i),XYZ_labels)) > 0; %strmatch(plotstr(i), textdata,'exact') > 0;
%             if strcmp(countstr(j),plotstr(i))==0;
%                 Draw_atoms_func(Init_coord{1,j},Resolution/2,Part_size(j)/30, Part_color{j},alfa); %(Coords,Resolution,Atom_size,color,alfa)
%             end
            if strcmp(plotstr(i),countstr(j))==1;
               Draw_atoms_func(Init_coord{1,j},Resolution,Part_size(j)/30, Part_color{j},1); %(Coords,Resolution,Atom_size,color,alfa)
            end
        end
    end
end
disp('Finished!!!');
hold off;

