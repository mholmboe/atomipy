"""
"Frozen dummy mineral" parameters (the Dummy FF) for inorganics not covered by
the built-in force fields (CLAYFF / MINFF).

For a material whose framework atom types the built-in force fields can't assign
(e.g. MnO, NiO, Cr2O3, …), this builds a crude-but-usable model so it can still
interact with water and solutes in a qualitative simulation:

  * Partial charges  = ``charge_scale`` (default 0.5) × the guessed oxidation
    state (see :func:`atomipy.guess_oxidation_states`). Half the formal charge
    mirrors the reduced charges of CLAYFF / MINFF and keeps a neutral lattice
    neutral.
  * Lennard-Jones    = borrowed: every oxygen gets the OPC3 water oxygen
    (σ=0.31743 nm, ε=0.68369 kJ/mol), and every metal/cation gets a small buried
    metal site (default ``Alo``, σ=0.14410 nm). Hydrogens get zero LJ.
  * The framework is **frozen** (atoms flagged ``frozen=True``; the OpenMM layer
    sets their mass to 0). Freezing removes the need for ANY bonded parameters —
    exactly what's missing for an unsupported framework — so only nonbonded
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
from ._provenance import provenance_string

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

# Universal Force Field (UFF) nonbonded parameters (Rappé et al., JACS 1992):
# element -> (x_i, D_i) with x_i the vdW distance (Å, = the LJ r_min) and D_i the
# well depth (kcal/mol). These are derived from atomic data and cover the whole
# periodic table, so the dummy FF can compute its OWN element-specific LJ instead
# of borrowing a single MINFF site for every metal. Conversion to GROMACS units:
#   sigma_nm = x_i / 2^(1/6) / 10 ;  epsilon_kJ = D_i * 4.184
UFF_VDW = {
    'H': (2.886, 0.044), 'He': (2.362, 0.056), 'Li': (2.451, 0.025), 'Be': (2.745, 0.085),
    'B': (4.083, 0.180), 'C': (3.851, 0.105), 'N': (3.660, 0.069), 'O': (3.500, 0.060),
    'F': (3.364, 0.050), 'Ne': (3.243, 0.042), 'Na': (2.983, 0.030), 'Mg': (3.021, 0.111),
    'Al': (4.499, 0.505), 'Si': (4.295, 0.402), 'P': (4.147, 0.305), 'S': (4.035, 0.274),
    'Cl': (3.947, 0.227), 'Ar': (3.868, 0.185), 'K': (3.812, 0.035), 'Ca': (3.399, 0.238),
    'Sc': (3.295, 0.019), 'Ti': (3.175, 0.017), 'V': (3.144, 0.016), 'Cr': (3.023, 0.015),
    'Mn': (2.961, 0.013), 'Fe': (2.912, 0.013), 'Co': (2.872, 0.014), 'Ni': (2.834, 0.015),
    'Cu': (3.495, 0.005), 'Zn': (2.763, 0.124), 'Ga': (4.383, 0.415), 'Ge': (4.280, 0.379),
    'As': (4.230, 0.309), 'Se': (4.205, 0.291), 'Br': (4.189, 0.251), 'Kr': (4.141, 0.220),
    'Rb': (4.114, 0.040), 'Sr': (3.641, 0.235), 'Y': (3.345, 0.072), 'Zr': (3.124, 0.069),
    'Nb': (3.165, 0.059), 'Mo': (3.052, 0.056), 'Tc': (2.998, 0.048), 'Ru': (2.963, 0.056),
    'Rh': (2.929, 0.053), 'Pd': (2.899, 0.048), 'Ag': (3.148, 0.036), 'Cd': (2.848, 0.228),
    'In': (4.463, 0.599), 'Sn': (4.392, 0.567), 'Sb': (4.420, 0.449), 'Te': (4.470, 0.398),
    'I': (4.500, 0.339), 'Xe': (4.404, 0.332), 'Cs': (4.517, 0.045), 'Ba': (3.703, 0.364),
    'La': (3.522, 0.017), 'Ce': (3.556, 0.013), 'Pr': (3.606, 0.010), 'Nd': (3.575, 0.010),
    'Sm': (3.520, 0.008), 'Eu': (3.496, 0.008), 'Gd': (3.368, 0.009), 'Tb': (3.451, 0.007),
    'Dy': (3.428, 0.007), 'Ho': (3.409, 0.007), 'Er': (3.391, 0.007), 'Tm': (3.374, 0.006),
    'Yb': (3.355, 0.228), 'Lu': (3.640, 0.041), 'Hf': (3.141, 0.072), 'Ta': (3.170, 0.081),
    'W': (3.069, 0.067), 'Re': (2.954, 0.066), 'Os': (3.120, 0.037), 'Ir': (2.840, 0.073),
    'Pt': (2.754, 0.080), 'Au': (3.293, 0.039), 'Hg': (2.705, 0.385), 'Tl': (4.347, 0.680),
    'Pb': (4.297, 0.663), 'Bi': (4.370, 0.518), 'Th': (3.396, 0.026), 'U': (3.395, 0.022),
}
_TWO_POW_NEG_SIXTH = 2.0 ** (-1.0 / 6.0)   # x_i (r_min) -> sigma
_KCAL_TO_KJ = 4.184


def uff_lj(element):
    """UFF Lennard-Jones (sigma_nm, epsilon_kJ/mol) for an element, computed from
    the UFF vdW distance/well-depth (sigma = x_i / 2^(1/6), epsilon = D_i·4.184).
    Returns None if the element is not in the UFF table."""
    entry = UFF_VDW.get(element)
    if entry is None:
        return None
    x_i, d_i = entry
    return (round(x_i * _TWO_POW_NEG_SIXTH / 10.0, 6), round(d_i * _KCAL_TO_KJ, 6))

# Element-appropriate metal Lennard-Jones (sigma_nm, epsilon_kJ/mol), for *pure
# metals / alloys* where the borrowed buried-cation site is inappropriate. The
# face-centred-cubic metals use the well-validated 12-6 parameters of Heinz et
# al., J. Phys. Chem. C 2008, 112, 17281 (good for metal–water/biomolecule
# interfaces with Lorentz-Berthelot mixing); the remaining metals use UFF
# (Rappe et al., JACS 1992) converted as sigma = x1 / 2^(1/6), eps = D1·4.184.
# Generic/approximate — intended for a qualitative frozen wall, not energetics.
ELEMENT_LJ = {
    # fcc metals — Heinz 2008 12-6
    'Al': (0.25527, 16.82), 'Ni': (0.22200, 23.64), 'Cu': (0.22770, 19.75),
    'Pd': (0.24510, 25.73), 'Ag': (0.25740, 19.08), 'Pt': (0.24720, 32.64),
    'Au': (0.25690, 22.13), 'Pb': (0.31190, 12.26),
    # other common metals — UFF
    'Li': (0.19457, 0.10460), 'Na': (0.23681, 0.12552), 'K': (0.30255, 0.14644),
    'Mg': (0.23975, 0.46442), 'Ca': (0.26977, 0.99579), 'Sc': (0.26157, 0.07950),
    'Ti': (0.25204, 0.07113), 'V': (0.24954, 0.06694), 'Cr': (0.23992, 0.06276),
    'Mn': (0.23502, 0.05439), 'Fe': (0.23110, 0.05439), 'Co': (0.22798, 0.05858),
    'Zn': (0.21934, 0.51882), 'Zr': (0.24793, 0.28886), 'Mo': (0.24224, 0.23430),
    'Cd': (0.25373, 0.95395), 'W': (0.24358, 0.28033), 'Hg': (0.24099, 1.61084),
}

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


def _element_lj(element):
    """Best self-calculated LJ (sigma_nm, epsilon_kJ) for an element: the curated
    ELEMENT_LJ value (Heinz metals + selected UFF) if present, else UFF computed
    from vdW data, else None."""
    if element in ELEMENT_LJ:
        return ELEMENT_LJ[element]
    return uff_lj(element)


def assign_dummy_mineral_params(atoms, Box=None, charge_mode='pauling', charge_scale=0.5,
                                h_charge=0.4, lj_mode='element', metal_site='Alo',
                                resname='DUM', rmaxlong=2.45, rmaxH=1.2,
                                freeze=True, verbose=True):
    """Assign Dummy FF parameters to a mineral framework in place.

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
    lj_mode : {'element', 'minff'}
        Lennard-Jones source. 'element' (default): the Dummy FF computes its OWN
        per-element LJ from vdW data — ELEMENT_LJ (Heinz metals + selected UFF)
        where available, else UFF (σ = x_i/2^(1/6), ε = D_i) for every element
        including O/F/H. 'minff': borrow from MINFF (O→OPC3, F→F⁻, H→none,
        metals→`metal_site`), which gives stronger O–water attraction.
    metal_site : str
        MINFF LJ site used in 'minff' mode (and as the fallback for elements in
        no LJ table): 'Alo' (default), 'Sit', or 'Mgo'.
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

    # Only the (non-MINFF) mineral FRAMEWORK is dummy-parameterised. Water, ions
    # and organic molecules are left EXACTLY untouched — they keep the standard
    # MINFF/min.ff water & ion parameters and their own GAFF/OpenFF itps, i.e. the
    # identical force field they'd have in a mixed MINFF simulation.
    from .composition import classify_atom
    _skip = {'water', 'ion', 'organic'}
    framework = [a for a in atoms if classify_atom(a) not in _skip]
    n_skipped = len(atoms) - len(framework)
    if not framework:
        if verbose:
            print("[dummy mineral] no framework atoms to parameterise "
                  "(all atoms are water/ions/organics — left untouched).")
        return atoms, {'non_minff_elements': [], 'net_charge': 0.0, 'n_atoms': 0,
                       'n_skipped': n_skipped, 'lj_mode': lj_mode,
                       'lj_source': 'n/a', 'metal_site': metal_site,
                       'charge_mode': charge_mode}

    # Oxidation states (ionic engine, neutral lattice by default) — framework only.
    ox = guess_oxidation_states(framework, method='ionic', write=True)
    elements = [_norm_element(a) for a in framework]

    # --- Charges ---
    anion_idx = []
    for i, (atom, o) in enumerate(zip(framework, ox)):
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
                _anion_charges_minff(framework, ox, anion_idx, Box, verbose,
                                     rmaxlong=rmaxlong, rmaxH=rmaxH)
                used_minff = True
            except Exception as exc:
                if verbose:
                    print(f"[dummy mineral] could not compute coordination ({exc}); "
                          f"falling back to |oxidation|-weighted charge balance.")
        if not used_minff:
            # Fallback (no Box / no coordination): anions absorb the remainder,
            # split by |oxidation| so the framework is neutral.
            assigned = sum(a['charge'] for a in framework)
            weights = [abs(float(ox[i])) or 1.0 for i in anion_idx]
            wsum = sum(weights) or 1.0
            for i, w in zip(anion_idx, weights):
                framework[i]['charge'] = round(-assigned * w / wsum, 6)

    # --- LJ, type, mass, freeze ---
    # lj_mode controls where Lennard-Jones parameters come from:
    #   'element' (default) — the Dummy FF's OWN per-element LJ, computed from vdW
    #       data: ELEMENT_LJ (Heinz metals + selected UFF) where available, else
    #       UFF (sigma = x_i/2^(1/6), epsilon = D_i) for every element incl. O/F/H.
    #       No MINFF borrowing; each element gets its own size.
    #   'minff' — borrow from MINFF: O→OPC3 oxygen, F→F⁻, H→none, metals→the
    #       small buried-cation site (`metal_site`). Stronger O–water attraction.
    # The curated metallic ELEMENT_LJ (Heinz 12-6) only makes sense for a PURE
    # metal/alloy. In an ionic framework (anions present) those deep metallic
    # wells are wrong for a cation, so use the UFF vdW LJ for every element.
    _pure_metal = not any(float(o) < 0 for o in ox)
    non_minff = set()
    fallback_elems = set()
    for i, atom in enumerate(framework):
        el = elements[i]
        if el not in MINFF_FRAMEWORK_ELEMENTS:
            non_minff.add(el)
        if lj_mode == 'minff':
            if el == 'O':
                atom['sigma'], atom['epsilon'] = MINFF_LJ_SITES['O_opc3']
            elif el == 'F':
                atom['sigma'], atom['epsilon'] = MINFF_LJ_SITES['F_ion']
            elif el == 'H':
                atom['sigma'], atom['epsilon'] = MINFF_LJ_SITES['H']
            else:
                atom['sigma'], atom['epsilon'] = metal_sigma, metal_eps
        else:  # 'element' — self-calculated per-element LJ from vdW data
            lj = _element_lj(el) if _pure_metal else uff_lj(el)
            if lj is None:                      # exotic element not in any table
                lj = (metal_sigma, metal_eps)
                fallback_elems.add(el)
            atom['sigma'], atom['epsilon'] = lj
        atom['_dummy_type'] = f"{el}d"
        atom['mass'] = _ATOMIC_MASS.get(el, 0.0)
        atom['resname'] = resname
        atom['molid'] = 1
        if freeze:
            atom['frozen'] = True

    net_charge = round(sum(a['charge'] for a in framework), 6)
    if lj_mode == 'minff':
        _lj_desc = f"MINFF-borrowed (O=OPC3, F=F⁻, metals=site '{metal_site}')"
    else:
        _lj_desc = "self-calculated per-element (ELEMENT_LJ/UFF)"
        if fallback_elems:
            _lj_desc += f"; borrowed '{metal_site}' for {sorted(fallback_elems)}"
    report = {
        'non_minff_elements': sorted(non_minff),
        'net_charge': net_charge,
        'n_atoms': len(framework),
        'n_skipped': n_skipped,
        'lj_mode': lj_mode,
        'lj_source': _lj_desc,
        'metal_site': metal_site,
        'charge_mode': charge_mode,
    }
    if verbose:
        if non_minff:
            _q = (f"Pauling effective charges (H={h_charge:+g})" if charge_mode != 'half'
                  else f"{charge_scale}×oxidation state")
            print(f"WARNING: Dummy FF — element(s) {sorted(non_minff)} are not "
                  f"covered by the built-in force fields. Building a FROZEN DUMMY "
                  f"model (charges = {_q}, LJ = {_lj_desc}). "
                  f"Qualitative only; run EM/NVT, not NPT.")
        _skip_note = (f"; {n_skipped} water/ion/organic atom(s) left untouched (MINFF)"
                      if n_skipped else "")
        print(f"[dummy mineral] {len(framework)} framework atoms, net charge "
              f"{net_charge:+.3f} e, charge '{charge_mode}', LJ = {_lj_desc}{_skip_note}.")
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
        f.write(f"; {provenance_string()} — Dummy FF mineral topology\n")
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
        f.write(f'; {provenance_string()} — frozen dummy-mineral system (qualitative; EM/NVT only)\n\n')
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
