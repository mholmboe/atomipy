"""
"Frozen dummy mineral" parameters for inorganics NOT supported by MINFF.

For a material whose framework atom types MINFF can't assign (e.g. MnO, NiO,
Cr2O3, …), this builds a crude-but-usable model so it can still interact with
water and solutes in a qualitative simulation:

  * Partial charges  = ``charge_scale`` (default 0.5) × the guessed oxidation
    state (see :func:`atomipy.guess_oxidation_states`). Half the formal charge
    mirrors the reduced charges of CLAYFF / MINFF and keeps a neutral lattice
    neutral.
  * Lennard-Jones    = borrowed from MINFF: every oxygen gets the MINFF/OPC3
    oxygen (σ=0.31743 nm, ε=0.68369 kJ/mol — MINFF already uses OPC3-O for all
    mineral oxygens), and every metal/cation gets a small buried metal site
    (default ``Alo``, σ=0.14410 nm). Hydrogens get zero LJ (MINFF convention).
  * The framework is **frozen** (atoms flagged ``frozen=True``; the OpenMM layer
    sets their mass to 0). Freezing removes the need for ANY bonded parameters —
    exactly what MINFF lacks for the unsupported mineral — so only nonbonded
    (electrostatics + LJ) terms are needed. Run **EM and NVT only** (a frozen
    rigid body is incompatible with an NPT barostat).

The resulting per-atom params are written into a self-contained GROMACS ``.itp``
(its own ``[ atomtypes ]`` + a bond-free ``[ moleculetype ]``) by
:func:`write_dummy_mineral_itp`, which rides the same ``#include`` path as an
organic (GAFF) itp.

This is intentionally approximate — a *dummy* for qualitative questions
(wetting, ion adsorption, interfacial structuring), not quantitative energetics.
"""
import re

from .oxidation import guess_oxidation_states, _norm_element

# Borrowable MINFF Lennard-Jones sites: (sigma_nm, epsilon_kJ_per_mol).
MINFF_LJ_SITES = {
    'O_opc3': (0.31743, 0.68369),  # MINFF/OPC3 oxygen — used for every dummy O
    'Alo':    (0.14410, 0.47389),  # Al octahedral — smallest *safe* metal site
    'Sit':    (0.08223, 0.42242),  # Si tetrahedral — smallest cation site overall
    'Mgo':    (0.19555, 0.73412),  # a roomier divalent-ish option
    'H':      (0.00000, 0.00000),  # MINFF hydrogens carry no LJ
}

# Elements MINFF can type as a framework site (everything else → "dummy needed").
MINFF_FRAMEWORK_ELEMENTS = {'Si', 'Al', 'Mg', 'Fe', 'Ca', 'Ti', 'Li', 'O', 'H'}

# Atomic masses (g/mol) for the topology. Frozen atoms get mass 0 in OpenMM, but
# the .itp keeps real masses so the same topology is valid for a GROMACS export.
_ATOMIC_MASS = {
    'H': 1.008, 'Li': 6.941, 'B': 10.811, 'C': 12.011, 'N': 14.007, 'O': 15.999,
    'F': 18.998, 'Na': 22.990, 'Mg': 24.305, 'Al': 26.982, 'Si': 28.086,
    'P': 30.974, 'S': 32.06, 'Cl': 35.45, 'K': 39.098, 'Ca': 40.078,
    'Sc': 44.956, 'Ti': 47.867, 'V': 50.942, 'Cr': 51.996, 'Mn': 54.938,
    'Fe': 55.845, 'Co': 58.933, 'Ni': 58.693, 'Cu': 63.546, 'Zn': 65.38,
    'Ga': 69.723, 'Ge': 72.63, 'As': 74.922, 'Se': 78.971, 'Br': 79.904,
    'Rb': 85.468, 'Sr': 87.62, 'Y': 88.906, 'Zr': 91.224, 'Nb': 92.906,
    'Mo': 95.95, 'Ag': 107.868, 'Cd': 112.414, 'In': 114.818, 'Sn': 118.71,
    'Sb': 121.76, 'Te': 127.6, 'I': 126.904, 'Cs': 132.905, 'Ba': 137.327,
    'La': 138.905, 'W': 183.84, 'Pt': 195.084, 'Au': 196.967, 'Hg': 200.592,
    'Pb': 207.2, 'Bi': 208.98,
}

# GROMACS atomic numbers for the [ atomtypes ] at.num column.
_ATOMIC_NUMBER = {
    'H': 1, 'Li': 3, 'B': 5, 'C': 6, 'N': 7, 'O': 8, 'F': 9, 'Na': 11, 'Mg': 12,
    'Al': 13, 'Si': 14, 'P': 15, 'S': 16, 'Cl': 17, 'K': 19, 'Ca': 20, 'Sc': 21,
    'Ti': 22, 'V': 23, 'Cr': 24, 'Mn': 25, 'Fe': 26, 'Co': 27, 'Ni': 28,
    'Cu': 29, 'Zn': 30, 'Ga': 31, 'Ge': 32, 'As': 33, 'Se': 34, 'Br': 35,
    'Rb': 37, 'Sr': 38, 'Y': 39, 'Zr': 40, 'Nb': 41, 'Mo': 42, 'Ag': 47,
    'Cd': 48, 'In': 49, 'Sn': 50, 'Sb': 51, 'Te': 52, 'I': 53, 'Cs': 55,
    'Ba': 56, 'La': 57, 'W': 74, 'Pt': 78, 'Au': 79, 'Hg': 80, 'Pb': 82, 'Bi': 83,
}


def assign_dummy_mineral_params(atoms, charge_scale=0.5, metal_site='Alo',
                                resname='DUM', freeze=True, verbose=True):
    """Assign dummy (non-MINFF) parameters to a mineral framework in place.

    Sets on every atom: ``charge`` (= charge_scale × oxidation state),
    ``sigma``/``epsilon`` (nm, kJ/mol), ``_dummy_type`` (a synthetic GROMACS
    atomtype name), ``mass``, and ``frozen`` (if ``freeze``). Also sets
    ``resname`` so the whole framework is one moleculetype.

    Parameters
    ----------
    atoms : list of dict
        The mineral atoms (need 'element' or a typeable 'type').
    charge_scale : float
        Multiplier on the oxidation state to get the partial charge (default 0.5).
    metal_site : str
        Which MINFF LJ site to borrow for metal/cation atoms
        ('Alo' default, 'Sit', 'Mgo', …; see ``MINFF_LJ_SITES``).
    resname : str
        Residue / moleculetype name for the dummy framework (default 'DUM').
    freeze : bool
        If True (default), flag atoms ``frozen=True`` (OpenMM sets mass 0).
    verbose : bool
        Print a warning naming the non-MINFF elements and a summary.

    Returns
    -------
    (atoms, report) : tuple
        report = {'non_minff_elements', 'net_charge', 'n_atoms', 'metal_site'}.
    """
    if metal_site not in MINFF_LJ_SITES or metal_site in ('O_opc3', 'H'):
        raise ValueError(f"metal_site must be one of "
                         f"{[k for k in MINFF_LJ_SITES if k not in ('O_opc3', 'H')]}")

    metal_sigma, metal_eps = MINFF_LJ_SITES[metal_site]
    o_sigma, o_eps = MINFF_LJ_SITES['O_opc3']

    # Oxidation states → charges (ionic engine, neutral lattice by default).
    ox = guess_oxidation_states(atoms, method='ionic', write=True)

    non_minff = set()
    net_charge = 0.0
    for atom, o in zip(atoms, ox):
        el = _norm_element(atom)
        if el not in MINFF_FRAMEWORK_ELEMENTS:
            non_minff.add(el)
        q = round(charge_scale * float(o), 6)
        atom['charge'] = q
        net_charge += q
        if el == 'O':
            atom['sigma'], atom['epsilon'] = o_sigma, o_eps
        elif el == 'H':
            atom['sigma'], atom['epsilon'] = MINFF_LJ_SITES['H']
        else:
            atom['sigma'], atom['epsilon'] = metal_sigma, metal_eps
        atom['_dummy_type'] = f"{el}d"            # synthetic atomtype, e.g. 'Mnd'
        atom['mass'] = _ATOMIC_MASS.get(el, 0.0)
        atom['resname'] = resname
        atom['molid'] = 1
        if freeze:
            atom['frozen'] = True

    report = {
        'non_minff_elements': sorted(non_minff),
        'net_charge': round(net_charge, 6),
        'n_atoms': len(atoms),
        'metal_site': metal_site,
    }
    if verbose:
        if non_minff:
            print(f"WARNING: not MINFF-compatible — element(s) {sorted(non_minff)} "
                  f"have no MINFF framework type. Building a FROZEN DUMMY model "
                  f"(charges = {charge_scale}×oxidation state, O LJ = OPC3, metal "
                  f"LJ = {metal_site}). Qualitative only; run EM/NVT, not NPT.")
        print(f"[dummy mineral] {len(atoms)} atoms, net charge {report['net_charge']:+.3f} e, "
              f"metal site '{metal_site}'.")
        if abs(report['net_charge']) > 1e-3:
            print(f"  NOTE: net charge is non-zero ({report['net_charge']:+.3f} e) — "
                  f"add neutralizing ions before production.")
    return atoms, report


def _unique_dummy_types(atoms):
    """Ordered map {dummy_type_name: (element, sigma, epsilon, mass, at_num)}."""
    types = {}
    for a in atoms:
        name = a.get('_dummy_type')
        if not name or name in types:
            continue
        el = _norm_element(a)
        types[name] = (el, float(a.get('sigma', 0.0)), float(a.get('epsilon', 0.0)),
                       _ATOMIC_MASS.get(el, float(a.get('mass', 0.0))),
                       _ATOMIC_NUMBER.get(el, 0))
    return types


def write_dummy_mineral_itp(atoms, file_path, mol_name='DUM'):
    """Write a self-contained GROMACS .itp for a dummy mineral framework.

    Emits its own ``[ atomtypes ]`` (the borrowed LJ sites), a bond-free
    ``[ moleculetype ]`` and ``[ atoms ]`` (per-atom dummy type + scaled charge).
    No bonds/angles/dihedrals — the framework is held rigid by freezing, so the
    topology needs no MINFF bonded parameters. #include this like a GAFF itp.

    Call :func:`assign_dummy_mineral_params` first.
    """
    types = _unique_dummy_types(atoms)
    with open(file_path, 'w', encoding='utf-8') as f:
        f.write(f"; Dummy (non-MINFF) mineral topology written by atomipy\n")
        f.write(f"; FROZEN framework — nonbonded only (EM/NVT). Qualitative model.\n\n")

        f.write("[ atomtypes ]\n")
        f.write(";name  at.num   mass      charge  ptype     sigma       epsilon\n")
        for name, (el, sig, eps, mass, atn) in types.items():
            f.write(f"{name:<6} {atn:>4}  {mass:>9.5f}  {0.0:>8.5f}  A  "
                    f"{sig:>11.6f} {eps:>11.6f}\n")
        f.write("\n")

        f.write("[ moleculetype ]\n")
        f.write(";name            nrexcl\n")
        f.write(f" {mol_name:<15} 1\n\n")

        f.write("[ atoms ]\n")
        f.write(";   nr  type  resnr resnm  atom  cgnr     charge       mass\n")
        for i, a in enumerate(atoms, start=1):
            el = _norm_element(a)
            typ = a.get('_dummy_type', f"{el}d")
            q = float(a.get('charge', 0.0))
            mass = _ATOMIC_MASS.get(el, float(a.get('mass', 0.0)))
            aname = (a.get('type') or el)[:5]
            f.write(f"{i:>6} {typ:<5} {1:>5} {mol_name:<5} {aname:<5} {i:>5} "
                    f"{q:>11.6f} {mass:>11.5f}\n")
        f.write("\n")
    return file_path
