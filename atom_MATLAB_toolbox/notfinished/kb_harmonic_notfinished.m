%% kb_harmonic.m
% * This function 
% * Tested 15/04/2018
% * Please report problems/bugs to michael.holmboe@umu.se

%% Not tested!!!

%% Examples
% * kb = kb_harmonic(kb,Atom_label1,varargin)

if nargin == 2
    Atom_label2=Atom_label1;
else
    Atom_label2=varargin{1};
end

mass1=mass_atom(Atom_label1);
mass2=mass_atom(Atom_label2);

reduced_mass=(mass1*mass2)/(mass1+mass2);

