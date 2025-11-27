"""
Import GROMACS .itp topology files.

This module mirrors the MATLAB `import_itp.m` utility and parses sections into
plain Python dictionaries (column name -> list of values) to make subsequent
editing straightforward.

Examples
--------
from atomipy import import_itp_topology
data = import_itp_topology("topology.itp")
data["atomtypes"]["name"][:3]  # ['opls_135', 'opls_136', 'opls_137']
"""

import os
from typing import Dict, List, Optional, Sequence, Tuple, Union

Number = Union[int, float]


def _strip_comment(line: str) -> str:
    """Remove inline comments starting with ';' or '#' and trim whitespace."""
    for marker in (";", "#"):
        if marker in line:
            line = line.split(marker, 1)[0]
    return line.strip()


def _to_number(value: str) -> Union[str, Number]:
    """Convert a string to int/float when possible, otherwise return the string."""
    try:
        number = float(value)
    except (TypeError, ValueError):
        return value

    # Preserve floating representation (avoid casting to int so we keep decimal precision,
    # especially for charges read from .itp files).
    return number


def _parse_table(
    lines: Sequence[str],
    columns: Sequence[str],
    numeric_columns: Optional[Sequence[str]] = None,
) -> Dict[str, List[Union[str, Number]]]:
    """Parse whitespace-separated table lines into a column dictionary."""
    numeric_columns = set(numeric_columns or [])
    rows: List[List[Union[str, Number]]] = []

    for raw in lines:
        tokens = _strip_comment(raw).split()
        if not tokens:
            continue

        padded = tokens + [""] * (len(columns) - len(tokens))
        row: List[Union[str, Number]] = []
        for name, token in zip(columns, padded):
            if name in numeric_columns and token != "":
                row.append(_to_number(token))
            else:
                row.append(token)
        rows.append(row)

    # Transpose into dict-of-lists
    column_data: Dict[str, List[Union[str, Number]]] = {name: [] for name in columns}
    for row in rows:
        for name, value in zip(columns, row):
            column_data[name].append(value)
    return column_data


def _choose_atomtype_columns(first_line: str) -> Tuple[List[str], List[str]]:
    """Select atomtype column labels based on the first data row."""
    tokens = _strip_comment(first_line).split()
    # Two common layouts:
    # 1) type name atnum charge ptype v w
    # 2) name atnum mass charge ptype sigma epsilon
    if len(tokens) >= 7 and tokens[2].replace(".", "", 1).isdigit():
        return (
            ["type", "name", "atnum", "charge", "ptype", "v", "w"],
            ["atnum", "charge", "v", "w"],
        )
    return (
        ["name", "atnum", "mass", "charge", "ptype", "sigma", "epsilon"],
        ["atnum", "mass", "charge", "sigma", "epsilon"],
    )


def import_itp(filename: str, itp_prev: Optional[Dict[str, Dict]] = None) -> Dict[str, Dict]:
    """
    Import a GROMACS .itp file into a dictionary of section dictionaries.

    Parameters
    ----------
    filename : str
        Path to the .itp file.
    itp_prev : dict, optional
        Previous itp dictionary to merge into the result (fields in `itp_prev`
        take precedence).

    Returns
    -------
    dict
        Mapping of section name to column dictionary (for tabular sections) or
        raw string (for ``enddata``).

    Examples
    --------
    data = import_itp("molecule.itp")
    list(data.keys())[:3]  # ['defaults', 'atomtypes', 'atoms']
    merged = import_itp("child.itp", itp_prev=data)  # merge with previous
    """
    if not os.path.isfile(filename):
        raise FileNotFoundError(f"Could not find itp file: {filename}")

    with open(filename, "r", encoding="utf-8") as fh:
        lines = [line.rstrip("\n") for line in fh]

    section_lines: Dict[str, List[str]] = {}
    current_section: Optional[str] = None
    section_counts: Dict[str, int] = {}

    for line in lines:
        stripped = line.strip()
        if not stripped or stripped.startswith(("#", ";")):
            continue
        if stripped.startswith("[") and stripped.endswith("]"):
            name = stripped.strip("[]").strip()
            count = section_counts.get(name, 0) + 1
            section_counts[name] = count
            if count > 1:
                name = f"{name}{count}"
            current_section = name
            section_lines[current_section] = []
            continue
        if current_section:
            section_lines[current_section].append(line)

    itp: Dict[str, Dict] = {}

    # moleculetype
    if "moleculetype" in section_lines:
        cols = ["moleculetype", "nrexcl"]
        itp["moleculetype"] = _parse_table(section_lines["moleculetype"], cols, ["nrexcl"])

    # atoms
    if "atoms" in section_lines:
        atom_cols = [
            "nr",
            "type",
            "resnr",
            "residue",
            "atom",
            "cgnr",
            "charge",
            "mass",
            "typeB",
            "chargeB",
            "massB",
            "comment",
        ]
        numeric = {"nr", "resnr", "cgnr", "charge", "mass", "chargeB", "massB"}
        itp["atoms"] = _parse_table(section_lines["atoms"], atom_cols, numeric)

    # atomtypes
    if "atomtypes" in section_lines:
        first_data = next((ln for ln in section_lines["atomtypes"] if _strip_comment(ln)), "")
        atomtype_cols, numeric_cols = _choose_atomtype_columns(first_data)
        itp["atomtypes"] = _parse_table(section_lines["atomtypes"], atomtype_cols, numeric_cols)
        # Normalize: expose 'type' consistently when only 'name' exists
        if "type" not in itp["atomtypes"] and "name" in itp["atomtypes"]:
            itp["atomtypes"]["type"] = list(itp["atomtypes"]["name"])

    # bondtypes / bonds
    bond_cols = ["ai", "aj", "funct", "c0", "c1", "c2", "c3"]
    bond_numeric = set(bond_cols)
    if "bondtypes" in section_lines:
        itp["bondtypes"] = _parse_table(section_lines["bondtypes"], bond_cols, bond_numeric)
    if "bonds" in section_lines:
        itp["bonds"] = _parse_table(section_lines["bonds"], bond_cols, bond_numeric)

    # angletypes / angles
    angle_cols = ["ai", "aj", "ak", "funct", "c0", "c1", "c2", "c3"]
    angle_numeric = set(angle_cols)
    if "angletypes" in section_lines:
        itp["angletypes"] = _parse_table(section_lines["angletypes"], angle_cols, angle_numeric)
    if "angles" in section_lines:
        itp["angles"] = _parse_table(section_lines["angles"], angle_cols, angle_numeric)

    # pairtypes / pairs
    pair_cols = ["ai", "aj", "funct", "c0", "c1", "c2", "c3"]
    pair_numeric = set(pair_cols)
    if "pairtypes" in section_lines:
        itp["pairtypes"] = _parse_table(section_lines["pairtypes"], pair_cols, pair_numeric)
    if "pairs" in section_lines:
        itp["pairs"] = _parse_table(section_lines["pairs"], pair_cols, pair_numeric)

    # exclusions
    if "exclusions" in section_lines:
        exc_cols = ["ai", "aj", "ak", "funct"]
        itp["exclusions"] = _parse_table(section_lines["exclusions"], exc_cols, exc_cols)

    # dihedral blocks
    dihedral_cols = ["ai", "aj", "ak", "al", "funct", "c0", "c1", "c2", "c3", "c4", "c5"]
    if "dihedrals" in section_lines:
        itp["dihedrals"] = _parse_table(section_lines["dihedrals"], dihedral_cols, dihedral_cols)
    if "dihedrals2" in section_lines:
        itp["dihedrals2"] = _parse_table(section_lines["dihedrals2"], dihedral_cols, dihedral_cols)
    if "dihedraltypes" in section_lines:
        itp["dihedraltypes"] = _parse_table(section_lines["dihedraltypes"], dihedral_cols, dihedral_cols)
    if "dihedraltypes2" in section_lines:
        itp["dihedraltypes2"] = _parse_table(section_lines["dihedraltypes2"], dihedral_cols, dihedral_cols)

    # impropers
    if "impropers" in section_lines:
        itp["impropers"] = _parse_table(section_lines["impropers"], dihedral_cols, dihedral_cols)

    # position restraints and trailing text
    if "position_restraints" in section_lines:
        itp["position_restraints"] = _parse_table(
            section_lines["position_restraints"],
            ["ai", "aj", "ak", "al", "funct", "c0", "c1", "c2", "c3", "c4", "c5"],
            dihedral_cols,
        )

    # Preserve raw tail if present
    if "enddata" in section_lines:
        itp["enddata"] = "\n".join(section_lines["enddata"])

    if itp_prev:
        itp.update(itp_prev)

    return itp
