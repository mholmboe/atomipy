"""GROMACS .mdp generation following MINFF conventions.

Templates mirror the MINFF lecture examples (mholmboe/atom, ATOM_scripts_lecture
/MINFF): Verlet + PME, rcoulomb = rvdw = 1.2 nm, ``periodic-molecules = yes``
(infinite clay sheets), and a ``define`` line selecting the angle-k variant,
ion model, flexible water, and position restraints.

Prototype simplification: thermostat / COM groups default to a single ``System``
group (no custom index file). Splitting into mineral vs water/ion groups is a
later (Phase-2) addition that needs an auto-generated ``.ndx``.
"""

# Defaults common to all dynamics stages (nm, ps, K, bar).
_COMMON = {
    "cutoff-scheme": "Verlet",
    "nstlist": "20",
    "rlist": "1.2",
    "coulombtype": "PME",
    "vdw-type": "Cut-off",
    "rcoulomb": "1.2",
    "rvdw": "1.2",
    "fourierspacing": "0.12",
    "pme-order": "4",
    "ewald-rtol": "1e-05",
    "pbc": "xyz",
    "periodic-molecules": "yes",
    "constraints": "none",
    "constraint_algorithm": "lincs",
    "lincs_order": "4",
    "lincs_iter": "1",
}


def build_defines(minff_variant="GMINFF_k500", ion_model=None, flexible=True, posres=None):
    """Assemble the ``define`` flag list, e.g. ['-DGMINFF_k500','-DOPC3_IOD_LM']."""
    flags = []
    if minff_variant:
        flags.append(f"-D{minff_variant}")
    if ion_model:
        flags.append(f"-D{ion_model}")
    if flexible:
        flags.append("-DFLEXIBLE")
    if posres:
        flags.append(f"-D{posres}")
    return flags


def _render(pairs):
    return "\n".join(f"{k:<22}= {v}" for k, v in pairs) + "\n"


def mdp(stage, *, defines=None, nsteps=None, dt=0.001, temperature=298.0,
        pressure=1.0, nstxtc=1000, nstenergy=100, nstlog=100,
        emtol=1000.0, emstep=0.01, gen_vel=None, continuation=None):
    """Return .mdp text for a stage: 'em', 'nvt', 'npt', or 'md'.

    Parameters
    ----------
    stage : str
        'em' (steepest-descent minimization), 'nvt', 'npt', or 'md'.
    defines : list of str, optional
        ``define`` flags (see :func:`build_defines`). Default GMINFF_k500+flexible.
    nsteps : int, optional
        Step count (stage defaults: em 5000, nvt/npt 20000, md 50000).
    dt : float
        Timestep in ps (default 0.001 = 1 fs).
    temperature, pressure : float
        Reference T (K) and P (bar) for coupling.
    """
    if defines is None:
        defines = build_defines()
    define_line = " ".join(defines)
    pairs = [("define", define_line)]

    if stage == "em":
        pairs += [
            ("integrator", "steep"),
            ("nsteps", str(nsteps if nsteps is not None else 5000)),
            ("emtol", str(emtol)),
            ("emstep", str(emstep)),
            # steepest-descent writes the minimization trajectory to .trr (it ignores
            # nstxout-compressed/.xtc), so set nstxout to capture the EM progression.
            ("nstxout", str(nstxtc)),
            ("nstxout-compressed", str(nstxtc)),
        ]
        pairs += list(_COMMON.items())
        pairs += [("DispCorr", "No")]
        return _render(pairs)

    # dynamics stages
    pairs += [
        ("integrator", "md"),
        ("nsteps", str(nsteps if nsteps is not None else (50000 if stage == "md" else 20000))),
        ("dt", str(dt)),
        ("nstcomm", "100"),
        ("comm-mode", "Linear"),
        ("nstxout-compressed", str(nstxtc)),
        ("nstenergy", str(nstenergy)),
        ("nstlog", str(nstlog)),
        ("continuation", continuation if continuation is not None else ("yes" if stage != "nvt" else "no")),
    ]
    pairs += list(_COMMON.items())
    # temperature coupling (single group for the prototype)
    pairs += [
        ("tcoupl", "V-rescale"),
        ("tc-grps", "System"),
        ("tau_t", "1.0"),
        ("ref_t", str(temperature)),
    ]
    # pressure coupling
    if stage == "npt":
        pairs += [
            ("pcoupl", "C-rescale"),
            ("pcoupltype", "semiisotropic"),
            ("tau_p", "2.0"),
            ("ref_p", f"{pressure} {pressure}"),
            ("compressibility", "4.5e-5 4.5e-5"),
            ("refcoord-scaling", "all"),
        ]
    else:
        pairs += [("pcoupl", "no")]
    # velocity generation: on for a fresh NVT start
    gv = gen_vel if gen_vel is not None else (stage == "nvt")
    if gv:
        pairs += [("gen_vel", "yes"), ("gen_temp", str(temperature)), ("gen_seed", "-1")]
    pairs += [("DispCorr", "EnerPres")]
    return _render(pairs)
