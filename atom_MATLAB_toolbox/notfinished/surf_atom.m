function surf_atom(atom,Box_dim,dimension,binsize,flip_n_mirror,limits,min_z,norm_conc)
%% This function draws prefered site plots


Data_OneColumn=Data_OneColumn+rand(size(Data_OneColumn))/100;
[X,Y] = meshgrid(0:binsize:(limits(2)), 0:binsize:(limits(4)));
edges{1}=[1:1:size(Y,1)];
edges{2}=[1:1:size(X,2)];
Coords=hist3([Data_OneColumn(:,2)/binsize Data_OneColumn(:,1)/binsize],'Edges',edges);
%Coords=hist3([Data_OneColumn(:,2)/binsize Data_OneColumn(:,1)/binsize],[size(Y,1),size(X,2)]);

if flip_n_mirror == 1 && dimension == 2
    disp('Flipping and mirroring the data')
    Coords=(Coords+flip(fliplr(Coords)))/2;
end
% 
h = fspecial('gaussian',[3 3], .5); %( 3 3 .5)
Coords=filter2(h, Coords);
sum(not(isnan(Data_OneColumn(:,1))))
size(Data_OneColumn,1)

% Commented 2015
%nAtoms=nAtoms*sum(not(isnan(Data_OneColumn(:,1))))/size(Data_OneColumn,1) % Because the data may contain nan's
records=sum(not(isnan(Data_OneColumn(:,1))))/nAtoms %size(Data_OneColumn,1)/nAtoms
dV=binsize^2*(limits(6)-limits(5))*1E-27 % dm^3
Coords = Coords/(records*dV*6.022e23);
Coords(Coords == 0) = NaN;

% if 7*median(Coords(isfinite(Coords))) < 10;
%     norm_conc=ceil(6*median(Coords(isfinite(Coords))))
% elseif 7*median(Coords(isfinite(Coords))) > 500;
%     norm_conc=500;
% else
%     norm_conc=100*ceil(6/100*median(Coords(isfinite(Coords))))
% end

Coords(Coords > norm_conc) = norm_conc-.01;
% Coords = Coords/norm_conc;

[x_hi y_hi] = size(Coords);
x_hi=ceil(x_hi*binsize/10)*10;
y_hi=ceil(y_hi*binsize/10)*10;
%axis_limits=[-2 y_hi/(1+2*flip_n_mirror) -2 x_hi 0 norm_conc];
axis_limits=[-2 y_hi -2 x_hi 0 norm_conc];

assignin('caller','Coords',Coords);
assignin('caller','max_Conc',norm_conc);
assignin('caller','X',X);
assignin('caller','Y',Y);

% Plots the prefered site surface
hold on;
rotate3d on; FigColor = [1 1 1];
set(gca,'PlotBoxAspectRatio',[(axis_limits(2)-axis_limits(1))/(axis_limits(6)-axis_limits(5)) (axis_limits(4)-axis_limits(3))/(axis_limits(6)-axis_limits(5)) 1],'Color',FigColor,'FontSize',28, 'XTick',0:10:axis_limits(2),'YTick',0:5:axis_limits(4), 'ZTick',0:ceil(norm_conc/5):ceil(norm_conc/5)*5); % ,   ,'fontweight','b');
pbaspect([(axis_limits(2)-axis_limits(1))/(axis_limits(6)-axis_limits(5)) (axis_limits(4)-axis_limits(3))/(axis_limits(6)-axis_limits(5)) 1])
set(gcf,'Position',[0 0 1280 800],'Color',FigColor);
xlabel('X (Å)'); ylabel('Y (Å)'); zlabel('Z (Å)');
% Commented 2015
%surf(X,Y,Coords*norm_conc,'EdgeColor','none');
surf(X,Y,Coords,'EdgeColor','none');
cbr=colormap(Cmap_func);
brighten(0.3);
caxis([0, ceil(norm_conc)])
axis(axis_limits);
colorbar('YTick',0:ceil(norm_conc/5):ceil(norm_conc/5)*5,'YTicklabel',0:ceil(norm_conc/5):ceil(norm_conc/5)*5,'FontSize',28);
view([0,90]);
hold off;

