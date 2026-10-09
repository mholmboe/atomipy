"""
Preprocessor defines that select a MINFF parameter set in GROMACS.

Since MINFF v1.0 the angle force constant is a define of its own, and the mineral
and the force constant are separate concerns:

* **General (GMINFF)** parameters are selected by the angle force constant alone::

      define = -DMINFF_k500 -DOPC3 -DOPC3_IOD_LM

* **Tailored (TMINFF)** parameters need two defines, the mineral and the force
  constant, both required::

      define = -DMontmorillonite -DMINFF_k500 -DOPC3 -DOPC3_IOD_LM

The old ``-DGMINFF_k500`` and ``-DMontmorillonite_k500`` no longer select anything.
A define that selects nothing is a bad failure: ``grompp`` then stops on undeclared
atomtypes and not on the define, so the error does not point at the cause. This
module builds the right defines and turns the old spellings into the new ones
(with a ``FutureWarning``), so that they cannot silently select nothing.

Do not mix this up with the *JSON block keys* of the bundled parameter files (see
:mod:`atomipy.ffparams`): the general keys were renamed in the same way
(``GMINFF_k500`` -> ``MINFF_k500``), but the tailored keys are unchanged
(``Montmorillonite_k500`` stays as it is, and selecting that block is equivalent
to passing both defines).
"""
from __future__ import annotations

import re
import warnings
from typing import Iterable, List, Optional, Sequence, Union

#: Angle force constants (kJ/mol/rad2) that the parameter sets are made for.
FORCE_CONSTANTS = ("k0", "k250", "k500", "k1500")

_K = r"k(?:0|250|500|1500)"
# <Mineral>_k500, as in the tailored JSON block keys and the old combined defines.
# The mineral may itself contain underscores and hyphens (cis_Oct_Fe2_cis, Hectorite-F).
_COMBINED = re.compile(rf"^(?P<mineral>[A-Za-z][A-Za-z0-9_-]*)_(?P<k>{_K})$")
_MINERAL = re.compile(r"^[A-Za-z][A-Za-z0-9_-]*$")
_GENERAL = re.compile(rf"^MINFF_{_K}$")

#: The default parameter set: the general MINFF at the 500 kJ/mol/rad2 angle force constant.
DEFAULT_VARIANT = "MINFF_k500"


def _warn_old(old: str, new: Sequence[str]) -> None:
    flags = " ".join(f"-D{n}" for n in new)
    warnings.warn(
        f"'{old}' is the pre-v1.0 MINFF spelling and no longer selects any parameters in "
        f"GROMACS; using {flags} instead. Use the new spelling to silence this warning.",
        FutureWarning,
        stacklevel=4,
    )


def _normalize(name: str) -> List[str]:
    """Return the defines (without ``-D``) that one variant string stands for."""
    name = name.strip()
    if name.startswith("-D"):
        name = name[2:]
    if _GENERAL.match(name):
        return [name]
    m = _COMBINED.match(name)
    if m:
        mineral, k = m.group("mineral"), m.group("k")
        if mineral == "GMINFF":  # old general spelling: GMINFF_k500
            new = [f"MINFF_{k}"]
        elif mineral == "MINFF":
            return [name]
        else:  # old tailored spelling: Montmorillonite_k500 -> two defines
            new = [mineral, f"MINFF_{k}"]
        _warn_old(name, new)
        return new
    return [name]  # CLAYFF_EXT, OPC3, ... are defines of their own


def minff_defines(variant: Union[str, Iterable[str], None] = DEFAULT_VARIANT,
                  mineral: Optional[str] = None) -> List[str]:
    """Return the preprocessor defines (without ``-D``) that select a MINFF parameter set.

    Parameters
    ----------
    variant : str or sequence of str, optional
        The angle force constant define of the set, e.g. ``'MINFF_k500'`` (default), or
        ``None`` for none. A sequence is taken as several defines. The pre-v1.0
        spellings ``'GMINFF_k500'`` and ``'<Mineral>_k500'`` are accepted and turned
        into the new ones, with a ``FutureWarning``.
    mineral : str, optional
        The mineral of a tailored (TMINFF) set, e.g. ``'Montmorillonite'``. Needs a
        ``variant`` that is an angle force constant define.

    Returns
    -------
    list of str
        ``['MINFF_k500']`` for a general set and ``['Montmorillonite', 'MINFF_k500']``
        for a tailored one, in the order that they are best written.

    Raises
    ------
    ValueError
        For a tailored set that lacks the force constant, for a mineral that is not a
        valid define name, or if the mineral is given twice.
    """
    if variant is None or variant == "":
        names: List[str] = []
    elif isinstance(variant, str):
        names = _normalize(variant)
    else:
        names = [d for v in variant for d in _normalize(v)]

    if mineral:
        mineral = mineral.strip()
        if not _MINERAL.match(mineral) or _COMBINED.match(mineral):
            raise ValueError(
                f"mineral={mineral!r} is not a valid MINFF mineral name (give the mineral "
                "only, for example 'Montmorillonite', and the force constant as the variant)")
        if mineral in names:
            raise ValueError(f"mineral {mineral!r} was given twice")
        if len(names) != 1 or not _GENERAL.match(names[0]):
            raise ValueError(
                "a tailored (TMINFF) parameter set needs both a mineral and an angle force "
                f"constant define such as 'MINFF_k500'; got variant={variant!r}, "
                f"mineral={mineral!r}")
        names = [mineral] + names

    out: List[str] = []
    for n in names:  # keep the order, drop repeats
        if n not in out:
            out.append(n)
    return out
