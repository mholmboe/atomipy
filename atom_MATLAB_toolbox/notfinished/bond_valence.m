%% bond_valence.m
% * This function tries to calculate the bond valence values according to
% * http://www.iucr.org/resources/data/datasets/bond-valence-parameters
% * compiled by I. David Brown, McMaster University, Ontario, Canada
% * idbrown@mcmaster.ca
% * Data set bvparm2016.cif: 2016 version, (posted 2016-11-03)
% 
%
%% Version
% 3.00
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% # atom = bond_valence()
% # atom = bond_valence()
%
function atom = bond_valence(varargin)

valence_ion1=-100; % Dummy value
if nargin>4
    valence_ion1=varargin{3};
end

valence_ion2=-100; % Dummy value
if nargin>5
    valence_ion2=varargin{4};
end

if nargin==3
    load('bond_valence_values.mat');
    valence_ion1=-100; % Dummy value
    valence_ion2=-100; % Dummy value
else
    Ion_1=varargin{1};
    Ion_2=varargin{2};
    R0=varargin{3};
    b=varargin{4};
    Valence_1=varargin{5};
    Valence_2=varargin{6};
end

if nargin>9
    valence_ion1=varargin{7};
else
    valence_ion1=-100; % Dummy value
end
if nargin>10
    valence_ion2=varargin{8};
else
    valence_ion2=-100; % Dummy value
end

if strncmpi(ion1,'Hw',2)
    ion1='H';
end

if strncmpi(ion1,'Ow',2)
    ion1='O';
end

if strncmpi(ion2,'Hw',2)
    ion2='H';
end

if strncmpi(ion2,'Ow',2)
    ion2='O';
end

ind1=find(strcmpi(Ion_1,ion1));
ind2=find(strcmpi(Ion_2,ion2));

if valence_ion1>-50
    valence1_ind=find(Valence_1==valence_ion1);
    ind1=intersect(ind1,valence1_ind);
end

if valence_ion2>-50
    valence2_ind=find(Valence_2==valence_ion2);
    ind2=intersect(ind2,valence2_ind);
end

ind=find(ismember(ind1,ind2));
ind=ind1(ind);

if numel(ind)==0
    ind1=find(ismember(Ion_2,ion1));
    ind2=find(ismember(Ion_1,ion2));
    ind=find(ismember(ind1,ind2));
    ind=ind1(ind);
end

ion_1='Fe';
ion_2='O';
valence_ion1=3;
valence_ion2=-2;
load('bond_valence_values.mat');
bondvalence_data

[mean_bv,std_bv,bv,bvalue]=bond_valence_data(ion_1,ion_2,R,Ion_1,Ion_2,R0,b,Valence_1,Valence_2,valence_ion1,valence_ion2);
valence=sum(bv);
Rdiff=bvalue*log(valence/round2dec([valence])); % Rdiff calc the average valence and from R0 - R i.e. the ideal bond minus the actual bond distance

load('Revised_Shannon_radii.mat');
RevShannonRadii
% disp('    Mean   |  Median  |  std ')
% [mean([atom.valence]-round2dec([atom.valence])) median([atom.valence]-round2dec([atom.valence])) std([atom.valence]-round2dec([atom.valence]))]

