%
%	File 'mdout.mdp' was generated
%	By user: miho0052 (501)
%	On host: miho0052s_MacBook_Pro.local
%	At date: Mon Jan 18 12:18:33 2021
%
%	Created by:
%	                     :_) GROMACS _ gmx grompp, 2020.3 (_:
%	
%	Executable:   /usr/local/gromacs_2020.3/bin/gmx
%	Data prefix:  /usr/local/gromacs_2020.3
%	Working dir:  /Users/miho0052/Dropbox/MATLAB/Matlab_scripts/functions/atom/notfinished/LiveSetup/NVT
%	Command line:
%	  gmx grompp _f ./nvt.mdp _c prenvt.gro _r prenvt.gro _p /Users/miho0052/Dropbox/MATLAB/Matlab_scripts/functions/atom/notfinished/LiveSetup/MMT_clayff/topol.top _n /Users/miho0052/Dropbox/MATLAB/Matlab_scripts/functions/atom/notfinished/LiveSetup/MMT_clayff/index.ndx _o nvt.tpr _pp

% VARIOUS PREPROCESSING OPTIONS
% Preprocessor information: use cpp syntax.
% e.g.: _I/home/joe/doe _I/home/mary/roe

% e.g.: _DPOSRES _DFLEXIBLE (note these variable names are case sensitive)
mdp.define                   = '-DFLEXIBLE'

% RUN CONTROL PARAMETERS
mdp.integrator               = 'md'
% Start time and timestep in ps
mdp.tinit                    = 0
mdp.dt                       = 0.0001
mdp.nsteps                   = 50000
% For exact run continuation or redoing part of a run
mdp.init_step                = 0
% Part index is updated automatically on checkpointing (keeps files separate)
mdp.simulation_part          = 1
% mode for center of mass motion removal
mdp.comm_mode                = 'Linear'
% number of steps for center of mass motion removal
mdp.nstcomm                  = 10
% group(s) for center of mass motion removal
mdp.comm_grps                = simCOMgrps

% LANGEVIN DYNAMICS OPTIONS
% Friction coefficient (amu/ps) and random seed
mdp.bd_fric                  = 0
mdp.ld_seed                  = -1

% ENERGY MINIMIZATION OPTIONS
% Force tolerance and initial step_size
mdp.emtol                    = 10
mdp.emstep                   = 0.01
% Max number of iterations in relax_shells
mdp.niter                    = 20
% Step size (ps^2) for minimization of flexible constraints
mdp.fcstep                   = 0
% Frequency of steepest descents steps when doing CG
mdp.nstcgsteep               = 1000
mdp.nbfgscorr                = 10

% TEST PARTICLE INSERTION OPTIONS
mdp.rtpi                     = 0.05

% OUTPUT CONTROL OPTIONS
% Output frequency for coords (x), velocities (v) and forces (f)
mdp.nstxout                  = 0
mdp.nstvout                  = 0
mdp.nstfout                  = 0
% Output frequency for energies to log file and energy file
mdp.nstlog                   = 100
mdp.nstcalcenergy            = 100
mdp.nstenergy                = 100
% Output frequency and precision for .xtc file
mdp.nstxout_compressed       = 1000
mdp.compressed_x_precision   = 1000
% This selects the subset of atoms for the compressed
% trajectory file. You can select multiple groups. By
% default, all atoms will be written.
mdp.compressed_x_grps        = 'System'
% Selection of energy groups
mdp.energygrps               = 'System'

% NEIGHBORSEARCHING PARAMETERS
% cut_off scheme (Verlet: particle based cut_offs)
mdp.cutoff_scheme            = 'Verlet'
% nblist update frequency
mdp.nstlist                  = 20
% Periodic boundary conditions: xyz, no, xy
mdp.pbc                      = 'xyz'
mdp.periodic_molecules       = 'yes'
% Allowed energy error due to the Verlet buffer in kJ/mol/ps per atom,
% a value of _1 means: use rlist
mdp.verlet_buffer_tolerance  = 0.005
% nblist cut_off        
mdp.rlist                    = 1.0
% long_range cut_off for switched potentials

% OPTIONS FOR ELECTROSTATICS AND VDW
% Method for doing electrostatics
mdp.coulombtype              = 'PME'
mdp.coulomb_modifier         = 'Potential_shift_Verlet'
mdp.rcoulomb_switch          = 0
mdp.rcoulomb                 = 1.0
% Relative dielectric constant for the medium and the reaction field
mdp.epsilon_r                = 1
mdp.epsilon_rf               = 0
% Method for doing Van der Waals
mdp.vdw_type                 = 'Cut_off'
mdp.vdw_modifier             = 'Potential_shift_Verlet'
% cut_off lengths       
mdp.rvdw_switch              = 0
mdp.rvdw                     = 1.0
% Apply long range dispersion corrections for Energy and Pressure
mdp.DispCorr                 = 'EnerPres'

% Spacing for the PME/PPPM FFT grid
mdp.fourierspacing           = 0.16
% FFT grid size, when a value is 0 fourierspacing will be used
mdp.fourier_nx               = 0
mdp.fourier_ny               = 0
mdp.fourier_nz               = 0
% EWALD/PME/PPPM parameters
mdp.pme_order                = 4
mdp.ewald_rtol               = 1e-05
mdp.ewald_rtol_lj            = 0.001
mdp.lj_pme_comb_rule         = 'Geometric'
mdp.ewald_geometry           = '3d'
mdp.epsilon_surface          = 0
mdp.implicit_solvent         = 'no'

% OPTIONS FOR WEAK COUPLING ALGORITHMS
% Temperature coupling  
mdp.tcoupl                   = 'V-rescale'
mdp.nsttcouple               = -1
mdp.nh_chain_length          = 10
% mdp.print_nose_hoover_chain_variables = no
% Groups to couple separately
mdp.tc_grps                  = 'System'
% Time constant (ps) and reference temperature (K)
mdp.tau_t                    = 2
mdp.ref_t                    = 298
% pressure coupling     
mdp.pcoupl                   = 'no'
mdp.pcoupltype               = 'Isotropic'
mdp.nstpcouple               = -1
% Time constant (ps), compressibility (1/bar) and reference P (bar)
mdp.tau_p                    = 1
% mdp.compressibility          = 
% mdp.ref_p                    = 
% Scaling of reference coordinates, No, All or COM
% mdp.refcoord_scaling         = com

% % OPTIONS FOR QMMM calculations
% mdp.QMMM                     = no
% % Groups treated Quantum Mechanically
% mdp.QMMM_grps                = 
% % QM method             
% mdp.QMmethod                 = 
% % QMMM scheme           
% mdp.QMMMscheme               = normal
% % QM basisset           
% mdp.QMbasis                  = 
% % QM charge             
% mdp.QMcharge                 = 
% % QM multiplicity       
% mdp.QMmult                   = 
% % Surface Hopping       
% mdp.SH                       = 
% % CAS space options     
% mdp.CASorbitals              = 
% mdp.CASelectrons             = 
% mdp.SAon                     = 
% mdp.SAoff                    = 
% mdp.SAsteps                  = 
% % Scale factor for MM charges
% mdp.MMChargeScaleFactor      = 1
% 
% % SIMULATED ANNEALING  
% % Type of annealing for each temperature group (no/single/periodic)
% mdp.annealing                = 
% % Number of time points to use for specifying annealing in each group
% mdp.annealing_npoints        = 
% % List of times at the annealing points for each group
% mdp.annealing_time           = 
% % Temp. at each annealing point, for each group.
% mdp.annealing_temp           = 

% GENERATE VELOCITIES FOR STARTUP RUN
mdp.gen_vel                  = 'no'
mdp.gen_temp                 = 300
mdp.gen_seed                 = -1

% OPTIONS FOR BONDS    
mdp.constraints              = 'none'
% Type of constraint algorithm
mdp.constraint_algorithm     = 'Lincs'
% Do not constrain the start configuration
mdp.continuation             = 'no'
% Use successive overrelaxation to reduce the number of shake iterations
mdp.Shake_SOR                = 'no'
% Relative tolerance of shake
mdp.shake_tol                = '0.0001'
% Highest order in the expansion of the constraint coupling matrix
mdp.lincs_order              = 4
% Number of iterations in the final step of LINCS. 1 is fine for
% normal simulations, but use 2 to conserve energy in NVE runs.
% For energy minimization with constraints it should be 4 to 8.
mdp.lincs_iter               = 1
% Lincs will write a warning to the stderr if in one step a bond
% rotates over more degrees than
mdp.lincs_warnangle          = 30
% Convert harmonic bonds to morse potentials
mdp.morse                    = 'no'


