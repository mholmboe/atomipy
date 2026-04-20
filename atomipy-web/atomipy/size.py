"""
Atomic size utilities backed by the Revised Shannon radii table.


This module provides helper functions to retrieve crystal/ionic radii for single
atoms and to estimate bond distances for pairs of atoms using oxidation state
and coordination number.
"""

from __future__ import annotations

from typing import Dict, Optional, Tuple

from .bond_valence import load_shannon_radii

ShannonEntry = Dict[str, Optional[float]]


def get_radius(
    element: str,
    ox_state: int,
    coordination: int,
    radii: Optional[Dict[Tuple[str, int, int], ShannonEntry]] = None,
    prefer: str = "ionic",
) -> Optional[float]:
    """
    Return the ionic or crystal radius for an element/oxidation/coordination.

    Parameters
    ----------
    element : str
        Element symbol, e.g., "Al".
    ox_state : int
        Oxidation state, e.g., 3.
    coordination : int
        Coordination number, e.g., 6.
    radii : dict, optional
        Parsed Shannon radii lookup (defaults to loading from data file).
    prefer : {'ionic', 'crystal'}
        Which radius to return when both are available.

    Returns
    -------
    float or None
        Radius in Angstroms, or None if not found.

    Examples
    --------
    # Ionic radius for Al3+ in CN=6
    r_al = get_radius("Al", 3, 6)
    # Crystal radius for Si4+ in CN=4
    r_si = get_radius("Si", 4, 4, prefer="crystal")
    """
    if radii is None:
        radii = load_shannon_radii()
    key = (element, ox_state, coordination)
    entry = radii.get(key)
    if not entry:
        # Fallback: any coordination match with same element/ox
        for (el, ox, _cn), candidate in radii.items():
            if el == element and ox == ox_state:
                entry = candidate
                break
    if not entry:
        return None
    if prefer == "crystal":
        return entry.get("crystal_radius") or entry.get("ionic_radius")
    return entry.get("ionic_radius") or entry.get("crystal_radius")


def bond_distance(
    element1: str,
    ox1: int,
    coord1: int,
    element2: str,
    ox2: int,
    coord2: int,
    radii: Optional[Dict[Tuple[str, int, int], ShannonEntry]] = None,
    use_crystal: bool = False,
) -> Optional[float]:
    """
    Estimate a bond distance by summing radii for a pair of atoms.

    Parameters
    ----------
    element1, element2 : str
        Element symbols.
    ox1, ox2 : int
        Oxidation states.
    coord1, coord2 : int
        Coordination numbers.
    radii : dict, optional
        Parsed Shannon radii lookup (defaults to loading from data file).
    use_crystal : bool, optional
        If True, prefer crystal radii; otherwise prefer ionic radii.

    Returns
    -------
    float or None
        Estimated bond distance in Angstroms, or None if either radius is missing.

    Examples
    --------
    # Al-O bond distance estimate (ionic radii)
    d_al_o = bond_distance("Al", 3, 6, "O", -2, 4)
    # Si-O bond using crystal radii
    d_si_o = bond_distance("Si", 4, 4, "O", -2, 4, use_crystal=True)
    """
    prefer = "crystal" if use_crystal else "ionic"
    r1 = get_radius(element1, ox1, coord1, radii=radii, prefer=prefer)
    r2 = get_radius(element2, ox2, coord2, radii=radii, prefer=prefer)
    if r1 is None or r2 is None:
        return None
    return r1 + r2
