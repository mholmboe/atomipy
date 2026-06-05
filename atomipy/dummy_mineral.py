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
import math
import re

from .oxidation import guess_oxidation_states, _norm_element, PAULING_EN

# Borrowable MINFF Lennard-Jones sites: (sigma_nm, epsilon_kJ_per_mol).
MINFF_LJ_SITES = {
    'O_opc3': (0.31743, 0.68369),  # MINFF/OPC3 oxygen — used for every dummy O
    'Alo':    (0.14410, 0.47389),  # Al octahedral — smallest *safe* metal site
    'Sit':    (0.08223, 0.42242),  # Si tetrahedral — smallest cation site overall
    'Mgo':    (0.19555, 0.73412),  # a roomier divalent-ish option
    'H':      (0.00000, 0.00000),  # MINFF hydrogens carry no LJ
    'F_ion':  (0.30985, 0.69345),  # MINFF F- — fluoride anion (CaF2, F-hectorite)
}

# Elements MINFF can type as a framework site (everything else → "dummy needed").
# F is supported (fluorite CaF2, F-hectorite); O/H are the framework non-metals.
MINFF_FRAMEWORK_ELEMENTS = {'Si', 'Al', 'Mg', 'Fe', 'Ca', 'Ti', 'Li', 'O', 'H', 'F'}

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


def pauling_effective_charge(oxidation_state, element, ref_en=None):
    """Pauling effective charge: q_formal × [1 − exp(−¼ (χ_ref − χ_M)²)].

    The bracket is Pauling's fractional ionic character of the bond between the
    cation and the reference anion (oxygen by default). Returns ``oxidation_state``
    scaled by that ionicity. Only meaningful for cations (oxidation > 0); anions
    don't fit this form and are set by charge balance instead.
    """
    chi_ref = PAULING_EN['O'] if ref_en is None else ref_en
    chi_m = PAULING_EN.get(element)
    if chi_m is None:
        return 0.5 * oxidation_state           # fallback: half-formal
    f = 1.0 - math.exp(-0.25 * (chi_ref - chi_m) ** 2)
    return oxidation_state * f


def _anion_charges_minff(atoms, ox, anion_idx, Box, verbose, rmaxlong=2.45, rmaxH=1.2):
    """Coordination-resolved anion charges (the MINFF formula):

        q_anion = oxidation_anion + Σ_j (oxidation_j − partial_j) / CN_j

    where j runs over the charge donors (cations and H) coordinating the anion,
    and CN_j is donor j's number of coordinating anions. Each donor's charge
    deficit is shared equally over the anions it coordinates, so a hydroxyl O
    picks up its H's +0.6 deficit and a bridging O picks up shares from every
    metal it bridges. Exact neutrality follows for a neutral lattice.

    Mutates atoms[i]['charge'] for i in anion_idx. Coordination is computed via
    bond_angle with the metal (``rmaxlong``) and hydrogen (``rmaxH``) cutoffs —
    the same global cutoffs MINFF typing uses.
    """
    from .bond_angle import bond_angle
    # same_molecule_only=False: framework coordination is purely geometric and
    # must NOT depend on molid. An imported/uploaded structure may carry per-atom
    # or per-residue molids (and replication multiplies them), which would
    # otherwise make the default same-molecule filter find ZERO bonds and orphan
    # every cation. (MINFF sidesteps this by resetting molids first; we don't.)
    bond_angle(atoms, Box, rmaxM=rmaxlong, rmaxH=rmaxH,
               same_molecule_only=False, verbose=False)  # populate 'neigh' (0-based)
    anion_set = set(anion_idx)

    # CN of each donor = how many anions it coordinates.
    cn = {}
    for j, o in enumerate(ox):
        if float(o) > 0:                          # donor (cation or H)
            cn[j] = sum(1 for nb in (atoms[j].get('neigh') or []) if nb in anion_set)

    orphaned = 0
    for i in anion_idx:
        q = float(ox[i])
        for j in (atoms[i].get('neigh') or []):
            if float(ox[j]) > 0 and cn.get(j, 0) > 0:
                q += (float(ox[j]) - atoms[j]['charge']) / cn[j]
        atoms[i]['charge'] = round(q, 6)
    for j, c in cn.items():
        if c == 0:
            orphaned += 1
    if orphaned and verbose:
        print(f"[dummy mineral] {orphaned} cation(s) have no coordinating anion "
              f"within 2.45 Å — their charge deficit could not be redistributed "
              f"(net charge will be non-zero).")


def assign_dummy_mineral_params(atoms, Box=None, charge_mode='pauling', charge_scale=0.5,
                                h_charge=0.4, metal_site='Alo', resname='DUM',
                                rmaxlong=2.45, rmaxH=1.2, freeze=True, verbose=True):
    """Assign dummy (non-MINFF) parameters to a mineral framework in place.

    Sets on every atom: ``charge``, ``sigma``/``epsilon`` (nm, kJ/mol),
    ``_dummy_type`` (a synthetic GROMACS atomtype name), ``mass``, ``resname``,
    and ``frozen`` (if ``freeze``).

    Charge modes
    ------------
    'pauling' (default)
        Cations get a Pauling effective charge q_eff = oxidation ×
        [1 − exp(−¼(χ_O − χ_M)²)] (e.g. Si +1.79, Al +1.70, Mg +1.36,
        Ti +2.38, Fe²⁺ +0.95 / Fe³⁺ +1.43). H is fixed at ``h_charge`` (+0.4,
        the MINFF value). Anion (O, F) charges follow the MINFF coordination
        formula when ``Box`` is given — q_anion = oxidation + Σ_j (oxidation_j −
        partial_j)/CN_j over the coordinating donors j (cations and H), with
        CN_j the donor's anion coordination number. Without ``Box`` it falls
        back to an |oxidation|-weighted charge balance. Either way the framework
        is neutral.
    'half'
        Legacy: charge = ``charge_scale`` × oxidation state for every atom.

    Parameters
    ----------
    atoms : list of dict
    Box : list or None
        Simulation box (Cell or Box_dim). Enables the MINFF coordination-resolved
        anion charges in 'pauling' mode. Strongly recommended.
    charge_mode : {'pauling', 'half'}
    charge_scale : float
        Multiplier for the 'half' mode (default 0.5).
    h_charge : float
        Fixed hydrogen charge in 'pauling' mode (default 0.4).
    metal_site : str
        MINFF LJ site borrowed for metal/cation atoms ('Alo' default, 'Sit', 'Mgo').
    rmaxlong : float
        Metal–anion coordination cutoff in Å for the MINFF charge formula
        (default 2.45 — the MINFF global cutoff).
    rmaxH : float
        Hydrogen bond cutoff in Å (default 1.2 — the MINFF global cutoff).
    resname, freeze, verbose : see above.

    Returns
    -------
    (atoms, report) : tuple
        report = {'non_minff_elements', 'net_charge', 'n_atoms', 'metal_site',
                  'charge_mode'}.
    """
    if metal_site not in MINFF_LJ_SITES or metal_site in ('O_opc3', 'H', 'F_ion'):
        raise ValueError(f"metal_site must be one of "
                         f"{[k for k in MINFF_LJ_SITES if k not in ('O_opc3', 'H', 'F_ion')]}")

    metal_sigma, metal_eps = MINFF_LJ_SITES[metal_site]

    # Oxidation states (ionic engine, neutral lattice by default).
    ox = guess_oxidation_states(atoms, method='ionic', write=True)
    elements = [_norm_element(a) for a in atoms]

    # --- Charges ---
    anion_idx = []
    for i, (atom, o) in enumerate(zip(atoms, ox)):
        el = elements[i]
        if charge_mode == 'half':
            atom['charge'] = round(charge_scale * float(o), 6)
        elif el == 'H':
            atom['charge'] = round(float(h_charge), 6)
        elif o > 0:                                   # cation → Pauling effective charge
            atom['charge'] = round(pauling_effective_charge(float(o), el), 6)
        else:                                         # anion (O, F, …) → balance later
            atom['charge'] = 0.0
            anion_idx.append(i)
    if charge_mode != 'half' and anion_idx:
        used_minff = False
        if Box is not None:
            try:
                # MINFF coordination-resolved anion charges (q_O = -2 + Σ deficits).
                _anion_charges_minff(atoms, ox, anion_idx, Box, verbose,
                                     rmaxlong=rmaxlong, rmaxH=rmaxH)
                used_minff = True
            except Exception as exc:
                if verbose:
                    print(f"[dummy mineral] could not compute coordination ({exc}); "
                          f"falling back to |oxidation|-weighted charge balance.")
        if not used_minff:
            # Fallback (no Box / no coordination): anions absorb the remainder,
            # split by |oxidation| so the framework is neutral.
            assigned = sum(a['charge'] for a in atoms)
            weights = [abs(float(ox[i])) or 1.0 for i in anion_idx]
            wsum = sum(weights) or 1.0
            for i, w in zip(anion_idx, weights):
                atoms[i]['charge'] = round(-assigned * w / wsum, 6)

    # --- LJ, type, mass, freeze ---
    non_minff = set()
    for i, atom in enumerate(atoms):
        el = elements[i]
        if el not in MINFF_FRAMEWORK_ELEMENTS:
            non_minff.add(el)
        if el == 'O':
            atom['sigma'], atom['epsilon'] = MINFF_LJ_SITES['O_opc3']
        elif el == 'F':
            atom['sigma'], atom['epsilon'] = MINFF_LJ_SITES['F_ion']
        elif el == 'H':
            atom['sigma'], atom['epsilon'] = MINFF_LJ_SITES['H']
        else:
            atom['sigma'], atom['epsilon'] = metal_sigma, metal_eps
        atom['_dummy_type'] = f"{el}d"
        atom['mass'] = _ATOMIC_MASS.get(el, 0.0)
        atom['resname'] = resname
        atom['molid'] = 1
        if freeze:
            atom['frozen'] = True

    net_charge = round(sum(a['charge'] for a in atoms), 6)
    report = {
        'non_minff_elements': sorted(non_minff),
        'net_charge': net_charge,
        'n_atoms': len(atoms),
        'metal_site': metal_site,
        'charge_mode': charge_mode,
    }
    if verbose:
        if non_minff:
            _q = (f"Pauling effective charges (H={h_charge:+g})" if charge_mode != 'half'
                  else f"{charge_scale}×oxidation state")
            print(f"WARNING: not MINFF-compatible — element(s) {sorted(non_minff)} "
                  f"have no MINFF framework type. Building a FROZEN DUMMY model "
                  f"(charges = {_q}, O LJ = OPC3, metal LJ = {metal_site}). "
                  f"Qualitative only; run EM/NVT, not NPT.")
        print(f"[dummy mineral] {len(atoms)} atoms, net charge {net_charge:+.3f} e, "
              f"mode '{charge_mode}', metal site '{metal_site}'.")
        if abs(net_charge) > 1e-3:
            print(f"  NOTE: net charge is non-zero ({net_charge:+.3f} e) — "
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


_SOLVENT_RES = {'SOL', 'WAT', 'HOH', 'TIP3', 'OPC', 'OPC3', 'SPC', 'SPCE',
                'TIP4', 'TIP5'}
_WATER_FILE = {
    'spce': 'spce', 'spc': 'spc', 'tip3p': 'tip3p', 'opc3': 'opc3',
    'tip3pfb': 'tip3p-fb', 'opc': 'opc', 'tip4p': 'tip4p',
    'tip4pew': 'tip4pew', 'tip4pfb': 'tip4p-fb',
}


def _parse_itp_moltype(itp_path):
    """Return (moleculetype_name, n_atoms_per_molecule) from a GROMACS .itp."""
    name, natoms, section = None, 0, None
    with open(itp_path, 'r', encoding='utf-8', errors='ignore') as f:
        for line in f:
            s = line.split(';')[0].strip()
            if not s:
                continue
            if s.startswith('['):
                section = s.strip('[] ').lower()
                continue
            if section == 'moleculetype' and name is None:
                name = s.split()[0]
            elif section == 'atoms':
                natoms += 1
    return name, natoms


def write_dummy_system_top(atoms, box, out_top, out_gro, water_model='spce',
                           mol_name='DUM', dummy_itp='dummy.itp', organic_itps=None,
                           ion_model='SPCE_HFE_LM'):
    """Write a self-contained GROMACS .top (+ .gro) for a frozen dummy mineral
    plus optional organic molecules, water and monatomic ions — bypassing the
    MINFF mineral machinery.

    The dummy framework (atoms carrying ``_dummy_type``) is described by its own
    ``[ atomtypes ]`` + bond-free ``[ moleculetype ]``; organic molecules reuse
    their GAFF/OpenFF ``.itp`` (``organic_itps``); water and ions reuse the MINFF
    ``min.ff`` includes. Atoms are reordered framework → organic → water → ions
    so each ``[ molecules ]`` block is contiguous.

    Parameters
    ----------
    atoms : list of dict
        Framework (dummy) + organic + solvent/ion atoms.
    box : list
        Box (Cell or Box_dim).
    out_top, out_gro : str
        Output paths.
    water_model : str
        Water model for the include/define (default 'spce').
    mol_name : str
        Dummy moleculetype name (default 'DUM').
    dummy_itp : str
        Filename for the generated dummy ``.itp``.
    organic_itps : list of str, optional
        Organic ``.itp`` filenames to ``#include`` (each defines a moleculetype).

    Returns
    -------
    (ordered_atoms, n_frozen) : tuple
        The atoms in the order written to the ``.gro`` (framework first), and the
        number of leading frozen framework atoms — pass ``range(n_frozen)`` to
        ``system.setParticleMass(i, 0)`` after loading.
    """
    import os as _os
    from collections import OrderedDict
    from . import write_conf as _wc

    organic_itps = list(organic_itps or [])
    out_dir = _os.path.dirname(_os.path.abspath(out_top))

    frame = [a for a in atoms if a.get('_dummy_type')]
    rest = [a for a in atoms if not a.get('_dummy_type')]
    water = [a for a in rest if str(a.get('resname', '')).upper() in _SOLVENT_RES]
    nonwater = [a for a in rest if str(a.get('resname', '')).upper() not in _SOLVENT_RES]

    # Map each organic .itp to (moltype name, atoms/molecule), and write the dummy itp.
    org_info = []
    for oi in organic_itps:
        p = oi if _os.path.isabs(oi) else _os.path.join(out_dir, _os.path.basename(oi))
        try:
            mt, nat = _parse_itp_moltype(p)
            if mt and nat:
                org_info.append((_os.path.basename(oi), mt, nat))
        except OSError:
            pass
    natoms_to_moltype = {nat: mt for (_f, mt, nat) in org_info}

    itp_path = _os.path.join(out_dir, _os.path.basename(dummy_itp))
    write_dummy_mineral_itp(frame, itp_path, mol_name=mol_name)

    # Bucket the non-water solute by molid group: multi-atom groups are organic
    # molecules (keyed by their moleculetype), single-atom groups are monatomic
    # ions (keyed by resname). Atoms of each kind are kept CONTIGUOUS so the .gro
    # order matches [ molecules ] even with several different organics / ion types.
    groups = OrderedDict()
    for a in nonwater:
        groups.setdefault(a.get('molid'), []).append(a)
    org_by_type = OrderedDict()      # moltype -> atoms
    org_mol_counts = OrderedDict()   # moltype -> n molecules
    ion_by_res = OrderedDict()       # resname -> atoms
    for _mid, g in groups.items():
        if len(g) > 1:
            mt = natoms_to_moltype.get(len(g)) or (org_info[0][1] if org_info else None)
            if mt:
                org_by_type.setdefault(mt, []).extend(g)
                org_mol_counts[mt] = org_mol_counts.get(mt, 0) + 1
        else:
            rn = str(g[0].get('resname', '')).strip() or 'ION'
            ion_by_res.setdefault(rn, []).extend(g)
    organic = [a for v in org_by_type.values() for a in v]
    ions = [a for v in ion_by_res.values() for a in v]

    n_water = len(water) // 3
    water_define = water_model.lower().replace('/', '').replace('-', '').upper()
    wm_file = _WATER_FILE.get(water_model.lower().replace('/', '').replace('-', ''), 'spce')

    # Resolve the MINFF ion #define and remap names to the ion block's moleculetype
    # spelling (e.g. 'Na+' -> 'Na'), reusing write_merged_top's logic.
    ion_seq = [(rn, len(v)) for rn, v in ion_by_res.items()]
    ion_define = None
    ion_seq_out = ion_seq
    if ion_seq:
        from .merge_top import _resolve_ion_define, _remap_ion_molnames
        ion_define = _resolve_ion_define(ion_model, water_model,
                                         ion_molnames=[n for n, _ in ion_seq])
        ion_seq_out = _remap_ion_molnames(ion_seq, ion_define)

    with open(out_top, 'w', encoding='utf-8') as f:
        f.write('; Frozen dummy-mineral system written by atomipy (qualitative; EM/NVT only)\n\n')
        f.write('[ defaults ]\n')
        f.write('; nbfunc   comb-rule   gen-pairs   fudgeLJ   fudgeQQ\n')
        f.write('  1        2           yes         0.5       0.8333333333\n\n')
        f.write(f'#include "{_os.path.basename(dummy_itp)}"\n')   # dummy [atomtypes] + DUM moltype
        for oi in organic_itps:
            f.write(f'#include "{_os.path.basename(oi)}"\n')       # organic [atomtypes] + moltype
        f.write('\n')
        if n_water or ion_seq:
            f.write(f'#define {water_define}\n')
            if ion_define:
                f.write(f'#define {ion_define}\n')   # activates ion atomtypes + moleculetypes
            f.write('#include "min.ff/ffnonbonded.itp"\n')
            if n_water:
                f.write(f'#include "min.ff/{wm_file}.itp"\n')
            if ion_seq:
                f.write('#include "min.ff/ions.itp"\n')
            f.write('\n')
        f.write('[ system ]\n')
        f.write('Frozen dummy mineral + solvent\n\n')
        f.write('[ molecules ]\n')
        f.write(';Compound          nmols\n')
        # NOTE: molecule names MUST start at column 0 (no leading space). The
        # OpenMM MINFF loader text-translates ' Na '->' Na+ ' (to fix MINFF ion
        # atomtypes); a leading space here would mangle ' Na' in [ molecules ]
        # into 'Na+' while the ions.itp moleculetype stays 'Na' -> "Unknown
        # molecule type: Na+". Writing at column 0 (like write_merged_top) avoids it.
        # Order: framework, organics, ions, then water (SOL) last — the common
        # convention. Each block is contiguous and matches the .gro atom order.
        f.write(f'{mol_name:<18} 1\n')
        for mt, cnt in org_mol_counts.items():
            f.write(f'{mt:<18} {cnt}\n')
        for name, cnt in ion_seq_out:
            f.write(f'{name:<18} {cnt}\n')
        if n_water:
            f.write(f'{"SOL":<18} {n_water}\n')

    ordered = frame + organic + ions + water
    _wc.gro(ordered, box, out_gro)
    return ordered, len(frame)
