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
        out = subprocess.run([gmx, "--version"], capture_output=True, text=True, timeout=30)
        blob = (out.stdout or "") + (out.stderr or "")
        m = re.search(r"GROMACS version:?\s*(\S+)", blob)
        version = m.group(1) if m else "unknown"
    except Exception:
        version = "unknown"
    return {"path": path, "version": version, "command": gmx}


def _sanitize_minff(ffdir):
    """Comment out orphan generic water bonded types (Ow/Hw) in ffbonded.itp.

    The shipped min.ff references generic ``Ow``/``Hw`` water types in a few
    "for completeness" bonded entries, but no water model or [atomtypes] block
    declares those names (they use suffixed names like OW_opc3), so grompp
    rejects them. We patch only the *run-dir copy*; the canonical FF is untouched.
    Returns the number of lines commented.
    """
    fb = Path(ffdir) / "ffbonded.itp"
    if not fb.exists():
        return 0
    out, n = [], 0
    for ln in fb.read_text().splitlines():
        s = ln.strip()
        if s and not s.startswith(";") and not s.startswith("["):
            toks = re.split(r"\s+", s)
            if "Ow" in toks or "Hw" in toks:
                out.append("; [atomipy-sanitized orphan water type] " + ln)
                n += 1
                continue
        out.append(ln)
    if n:
        fb.write_text("\n".join(out) + "\n")
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


def stage_minff(workdir, *, minff_src=None, sanitize=True):
    """Copy min.ff into an existing run dir (which already has top/gro/organic itps).

    Lighter than stage_run_dir: only the force field. Use when the topology and
    coordinates were already written into ``workdir`` (e.g. by the web-module
    codegen). Returns the number of orphan lines sanitized.
    """
    wd = Path(workdir)
    src = Path(minff_src) if minff_src else _MINFF_SRC
    dst = wd / "min.ff"
    if dst.exists():
        shutil.rmtree(dst)
    shutil.copytree(src, dst)
    return _sanitize_minff(dst) if sanitize else 0


def run_local_gmx(workdir, top, gro, stages, *, defines=None, gmx="gmx",
                  minff_src=None, do_stage_minff=True, on_line=print, **mdp_kwargs):
    """High-level local-GROMACS run for codegen: stage min.ff, run stages, log.

    The .top/.gro (and any organic .itp) must already exist in ``workdir``.
    Copies min.ff in, runs the stage pipeline chaining each .gro forward, and
    calls ``on_line(str)`` for every log line (so a plain print() streams to SSE).
    Returns the list of per-stage status dicts.
    """
    if do_stage_minff:
        n = stage_minff(workdir, minff_src=minff_src)
        if n:
            on_line(f"[atomipy] sanitized {n} orphan water type(s) in run-dir min.ff")
    statuses = []
    for item in run_pipeline(workdir, stages, gro, defines=defines, gmx=gmx, top=top, **mdp_kwargs):
        if isinstance(item, dict):
            statuses.append(item)
        else:
            on_line(item)
    return statuses


def _stream(cmd, cwd, env):
    """Run cmd, yielding merged stdout/stderr lines; finally yield a status dict."""
    proc = subprocess.Popen(
        cmd, cwd=cwd, env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        text=True, bufsize=1,
    )
    for line in proc.stdout:
        yield line.rstrip("\n")
    proc.wait()
    yield {"returncode": proc.returncode}


def run_stage(workdir, stage, struct_in, *, defines=None, gmx="gmx", maxwarn=2,
              restraint=None, top="topol.top", **mdp_kwargs):
    """Generate the .mdp, run grompp then mdrun for one stage; yield log lines.

    Yields str log lines, then a final dict {'stage','returncode','gro'} where
    'gro' is the output coordinate file (``{stage}.gro``) if mdrun succeeded.
    """
    wd = Path(workdir)
    env = dict(os.environ)
    env.setdefault("GMX_MAXBACKUP", "-1")

    mdp_path = wd / f"{stage}.mdp"
    mdp_path.write_text(_mdp_text(stage, defines=defines, **mdp_kwargs))

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

    mdrun = [gmx, "mdrun", "-deffnm", stage]
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
