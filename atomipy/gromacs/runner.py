"""Local GROMACS runner: detect gmx, stage a run dir, run grompp/mdrun, stream log.

This is the engine layer for the local GROMACS path. The structure + topology
are produced by the usual atomipy writers (write_gro / write_gmx_top); this
module copies the ``min.ff`` force field into the run dir, writes the .mdp, and
drives ``gmx grompp``/``gmx mdrun`` as subprocesses whose merged stdout/stderr
is yielded line-by-line (for SSE streaming in the web app).
"""

import os
import re
import shutil
import subprocess
from pathlib import Path

from .mdp import mdp as _mdp_text, build_defines

# Canonical force-field dir shipped with atomipy: atomipy/ffparams/min.ff
_MINFF_SRC = Path(__file__).resolve().parent.parent / "ffparams" / "min.ff"


def detect_gmx(gmx="gmx"):
    """Return {'path','version','command'} for a usable gmx, or None if absent."""
    path = shutil.which(gmx)
    if not path:
        return None
    try:
        out = subprocess.run([gmx, "--version"], capture_output=True, text=True,
                             encoding="utf-8", errors="replace", timeout=30)
        blob = (out.stdout or "") + (out.stderr or "")
        m = re.search(r"GROMACS version:?\s*(\S+)", blob)
        version = m.group(1) if m else "unknown"
    except Exception:
        version = "unknown"
    return {"path": path, "version": version, "command": gmx}


def _active_atomtypes(ffdir, defines):
    """Set of [atomtypes] names active in ffnonbonded.itp under the given defines.

    Honors #ifdef/#ifndef/#else/#endif so the result reflects which mineral
    parameter block (GMINFF_k* vs CLAYFF_EXT, ion sets, …) is selected. Returns
    None if ffnonbonded.itp is absent.
    """
    fn = Path(ffdir) / "ffnonbonded.itp"
    if not fn.exists():
        return None
    defset = {(d[2:] if d.startswith("-D") else d) for d in (defines or [])}
    active, in_atomtypes, stack = set(), False, []  # stack: (kind, name)

    def ok():
        for kind, name in stack:
            if kind == "ifdef" and name not in defset:
                return False
            if kind == "ifndef" and name in defset:
                return False
        return True

    for raw in fn.read_text(encoding="utf-8", errors="ignore").splitlines():
        s = raw.strip()
        if s.startswith("#ifdef"):
            stack.append(("ifdef", s.split()[1])); continue
        if s.startswith("#ifndef"):
            stack.append(("ifndef", s.split()[1])); continue
        if s.startswith("#else") and stack:
            k, n = stack[-1]; stack[-1] = ("ifndef" if k == "ifdef" else "ifdef", n); continue
        if s.startswith("#endif"):
            if stack: stack.pop()
            continue
        if s.startswith("["):
            in_atomtypes = s.lower().startswith("[ atomtypes")
            continue
        if in_atomtypes and s and not s.startswith(";") and ok():
            active.add(s.split()[0])
    return active


def _sanitize_minff(ffdir, defines=None):
    """Comment out ffbonded.itp lines that reference an undeclared atom type.

    The shipped min.ff has a few bonded/angle entries whose atom types aren't
    declared under every parameter set — e.g. generic water Ow/Hw (no model
    declares them) or Feo2/Feo3 (GMINFF-only; CLAYFF_EXT declares generic Feo).
    grompp fatally rejects these even when the types are unused. When ``defines``
    is given we strip lines referencing any type not active under those defines
    (define-aware); otherwise we fall back to the legacy Ow/Hw check. We patch
    only the *run-dir copy* — the canonical FF is untouched. Returns #lines cut.
    """
    fb = Path(ffdir) / "ffbonded.itp"
    if not fb.exists():
        return 0
    active = _active_atomtypes(ffdir, defines) if defines else None
    defset = {(d[2:] if d.startswith("-D") else d) for d in (defines or [])}
    out, n, stack = [], 0, []  # stack: (kind, name) for #ifdef/#ifndef

    def guards_ok():
        # Mirror grompp's preprocessor: only consider a line if its enclosing
        # #ifdef/#ifndef guards are satisfied by the active defines. So a properly
        # #ifdef-guarded orphan is left alone (grompp already skips it).
        for kind, name in stack:
            if kind == "ifdef" and name not in defset:
                return False
            if kind == "ifndef" and name in defset:
                return False
        return True

    for ln in fb.read_text(encoding="utf-8", errors="ignore").splitlines():
        s = ln.strip()
        if s.startswith("#ifdef"):
            stack.append(("ifdef", s.split()[1])); out.append(ln); continue
        if s.startswith("#ifndef"):
            stack.append(("ifndef", s.split()[1])); out.append(ln); continue
        if s.startswith("#else") and stack:
            k, nm = stack[-1]; stack[-1] = ("ifndef" if k == "ifdef" else "ifdef", nm); out.append(ln); continue
        if s.startswith("#endif"):
            if stack: stack.pop()
            out.append(ln); continue
        if s and not s.startswith(";") and not s.startswith("[") and not s.startswith("#") and guards_ok():
            toks = re.split(r"\s+", s)
            type_toks = []
            for t in toks:                       # leading non-integer tokens = atom types
                if re.fullmatch(r"-?\d+", t):
                    break
                type_toks.append(t)
            if type_toks:
                bad = (any(t not in active for t in type_toks) if active is not None
                       else ("Ow" in type_toks or "Hw" in type_toks))
                if bad:
                    out.append("; [atomipy-sanitized orphan type] " + ln)
                    n += 1
                    continue
        out.append(ln)
    if n:
        fb.write_text("\n".join(out) + "\n", encoding="utf-8")
    return n


def stage_run_dir(atoms, Box, workdir, *, write_top, write_gro,
                  minff_src=None, top_name="topol.top", gro_name="conf.gro",
                  sanitize=True):
    """Create the run dir: copy min.ff, write topology + coordinates.

    ``write_top`` / ``write_gro`` are passed in (the atomipy writers) to keep
    this module import-light and testable. ``write_top(atoms, Box, path)`` and
    ``write_gro(atoms, Box, path)`` are called.
    Returns {'workdir','top','gro','sanitized'}.
    """
    wd = Path(workdir)
    wd.mkdir(parents=True, exist_ok=True)
    src = Path(minff_src) if minff_src else _MINFF_SRC
    dst = wd / "min.ff"
    if dst.exists():
        shutil.rmtree(dst)
    shutil.copytree(src, dst)
    sanitized = _sanitize_minff(dst) if sanitize else 0

    top_path = wd / top_name
    gro_path = wd / gro_name
    write_top(atoms, Box, str(top_path))
    write_gro(atoms, Box, str(gro_path))
    return {"workdir": str(wd), "top": str(top_path), "gro": str(gro_path), "sanitized": sanitized}


def stage_minff(workdir, *, minff_src=None, sanitize=True, defines=None):
    """Copy min.ff into an existing run dir (which already has top/gro/organic itps).

    Lighter than stage_run_dir: only the force field. Use when the topology and
    coordinates were already written into ``workdir`` (e.g. by the web-module
    codegen). ``defines`` (the active FF defines, e.g. ['GMINFF_k500','OPC3'] or
    ['CLAYFF_EXT', …]) makes sanitization define-aware. Returns #lines sanitized.
    """
    wd = Path(workdir)
    src = Path(minff_src) if minff_src else _MINFF_SRC
    dst = wd / "min.ff"
    if dst.exists():
        shutil.rmtree(dst)
    shutil.copytree(src, dst)
    return _sanitize_minff(dst, defines=defines) if sanitize else 0


def run_local_gmx(workdir, top, gro, stages, *, defines=None, gmx="gmx",
                  minff_src=None, do_stage_minff=True, san_defines=None,
                  on_line=print, **mdp_kwargs):
    """High-level local-GROMACS run for codegen: stage min.ff, run stages, log.

    The .top/.gro (and any organic .itp) must already exist in ``workdir``.
    Copies min.ff in, runs the stage pipeline chaining each .gro forward, and
    calls ``on_line(str)`` for every log line (so a plain print() streams to SSE).
    Returns the list of per-stage status dicts.
    """
    if do_stage_minff:
        n = stage_minff(workdir, minff_src=minff_src, defines=(san_defines if san_defines is not None else defines))
        if n:
            on_line(f"[atomipy] sanitized {n} orphan FF type line(s) in run-dir min.ff")
    statuses = []
    for item in run_pipeline(workdir, stages, gro, defines=defines, gmx=gmx, top=top, **mdp_kwargs):
        if isinstance(item, dict):
            statuses.append(item)
        else:
            on_line(item)
    return statuses


def trjconv_to_pdb(workdir, *, tpr, xtc, out, gmx="gmx", pbc="atom", skip=1,
                   group="System", on_line=None):
    """Convert a GROMACS trajectory (.xtc/.trr) to a multi-frame PDB via trjconv.

    Browser viewers (3Dmol/JSmol) can't read binary GROMACS trajectories, so we
    convert to a multi-MODEL PDB — the same format the OpenMM path produces.
    ``group`` is the trjconv output group (fed on stdin), e.g. 'System' or
    'non-Water'. Returns the output path, or None if the .xtc is missing/empty or
    trjconv fails (caller can fall back to a single final frame).
    """
    wd = Path(workdir)
    if not (wd / xtc).exists() or (wd / xtc).stat().st_size == 0:
        return None
    env = dict(os.environ)
    env.setdefault("GMX_MAXBACKUP", "-1")
    cmd = [gmx, "trjconv", "-s", tpr, "-f", xtc, "-o", out, "-pbc", pbc]
    if skip and int(skip) > 1:
        cmd += ["-skip", str(int(skip))]
    proc = subprocess.run(
        cmd, cwd=str(wd), env=env, input=f"{group}\n",
        capture_output=True, text=True, encoding="utf-8", errors="replace",
    )
    if on_line:
        for ln in (proc.stdout + proc.stderr).splitlines():
            on_line(ln)
    return str(wd / out) if proc.returncode == 0 and (wd / out).exists() else None


def _stream(cmd, cwd, env):
    """Run cmd, yielding merged stdout/stderr lines; finally yield a status dict."""
    proc = subprocess.Popen(
        cmd, cwd=cwd, env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        text=True, encoding="utf-8", errors="replace", bufsize=1,
    )
    for line in proc.stdout:
        yield line.rstrip("\n")
    proc.wait()
    yield {"returncode": proc.returncode}


def run_stage(workdir, stage, struct_in, *, defines=None, gmx="gmx", maxwarn=2,
              restraint=None, top="topol.top", ntmpi=1, ntomp=None, **mdp_kwargs):
    """Generate the .mdp, run grompp then mdrun for one stage; yield log lines.

    Yields str log lines, then a final dict {'stage','returncode','gro'} where
    'gro' is the output coordinate file (``{stage}.gro``) if mdrun succeeded.
    """
    wd = Path(workdir)
    env = dict(os.environ)
    env.setdefault("GMX_MAXBACKUP", "-1")

    mdp_path = wd / f"{stage}.mdp"
    mdp_path.write_text(_mdp_text(stage, defines=defines, **mdp_kwargs), encoding="utf-8")

    grompp = [gmx, "grompp", "-f", f"{stage}.mdp", "-c", struct_in, "-p", top,
              "-o", f"{stage}.tpr", "-maxwarn", str(maxwarn)]
    if restraint:
        grompp += ["-r", restraint]
    rc = 0
    yield f"$ {' '.join(grompp)}"
    for item in _stream(grompp, str(wd), env):
        if isinstance(item, dict):
            rc = item["returncode"]
        else:
            yield item
    if rc != 0:
        yield {"stage": stage, "returncode": rc, "gro": None, "step": "grompp"}
        return

    # -ntmpi 1: single thread-MPI rank (no domain decomposition). Robust for the
    # small / periodic-molecule systems typical here — auto DD often fails to split
    # a clay sheet that spans the cell across many ranks. OpenMP still parallelizes.
    # -v: verbose step progress to stderr (merged into the stream) so the run panel
    # shows live progress instead of going silent during the mdrun compute.
    mdrun = [gmx, "mdrun", "-deffnm", stage, "-v"]
    if ntmpi:
        mdrun += ["-ntmpi", str(int(ntmpi))]
    if ntomp:
        mdrun += ["-ntomp", str(int(ntomp))]
    yield f"$ {' '.join(mdrun)}"
    for item in _stream(mdrun, str(wd), env):
        if isinstance(item, dict):
            rc = item["returncode"]
        else:
            yield item
    gro = str(wd / f"{stage}.gro") if rc == 0 and (wd / f"{stage}.gro").exists() else None
    yield {"stage": stage, "returncode": rc, "gro": gro, "step": "mdrun"}


def run_pipeline(workdir, stages, struct_in, *, defines=None, gmx="gmx", **mdp_kwargs):
    """Run a sequence of stages, chaining each stage's output .gro to the next.

    ``stages`` is a list of stage names (e.g. ['em','nvt','npt','md']). Yields
    log lines and per-stage status dicts. Stops if any stage fails.
    """
    current = struct_in
    for stage in stages:
        last = None
        for item in run_stage(workdir, stage, current, defines=defines, gmx=gmx, **mdp_kwargs):
            if isinstance(item, dict):
                last = item
            yield item
        if not last or last.get("returncode", 1) != 0 or not last.get("gro"):
            return
        current = os.path.basename(last["gro"])
