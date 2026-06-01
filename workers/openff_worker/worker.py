from fastapi import FastAPI, HTTPException
from pydantic import BaseModel
import tempfile
import os
from pathlib import Path
import uuid
import json

from rdkit import Chem
from openff.toolkit import Molecule, ForceField, Topology
from openff.interchange import Interchange

app = FastAPI(title="OpenFF Worker Microservice")

CACHE_DIR = Path("/tmp/openff_cache")
CACHE_DIR.mkdir(parents=True, exist_ok=True)

class ParametrizeRequest(BaseModel):
    smiles: str
    forcefield: str = "openff_sage"

class MixRequest(BaseModel):
    organic_top: str
    organic_gro: str
    mineral_top: str
    mineral_gro: str
    box: list[float]

@app.post("/parametrize")
async def parametrize(req: ParametrizeRequest):
    try:
        # SMILES to OpenFF Molecule
        off_mol = Molecule.from_smiles(req.smiles, allow_undefined_stereo=True)
        off_mol.generate_conformers(n_conformers=1)
        
        # Parametrize with openff forcefield
        ff_name = req.forcefield
        if ff_name == "openff_sage":
            ff_name = "openff-2.0.0.offxml"
        elif ff_name == "openff_parsley":
            ff_name = "openff-1.0.0.offxml"
        elif ff_name == "openff_rosemary":
            ff_name = "openff-2.1.0.offxml"
            
        ff = ForceField(ff_name)
        topology = off_mol.to_topology()
        interchange = Interchange.from_smirnoff(ff, topology)
        
        job_id = str(uuid.uuid4())
        out_dir = CACHE_DIR / job_id
        out_dir.mkdir(parents=True, exist_ok=True)
        
        # We export to GROMACS format for interchange later
        interchange.to_gromacs(prefix=str(out_dir / "organic"))
        
        return {
            "job_id": job_id,
            "top": str(out_dir / "organic.top"),
            "gro": str(out_dir / "organic.gro")
        }
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

@app.post("/mix")
async def mix(req: MixRequest):
    """
    Experimental Interchange-based mixing.
    Feature-flagged by OPENFF_GROMACS_READER_STABLE.
    """
    if os.environ.get("OPENFF_GROMACS_READER_STABLE") != "1":
        raise HTTPException(status_code=501, detail="Interchange GROMACS reader not yet stable enough for mixing")
    
    # In future:
    # org_ic = Interchange.from_gromacs(req.organic_top, req.organic_gro)
    # min_ic = Interchange.from_gromacs(req.mineral_top, req.mineral_gro)
    # merged = org_ic + min_ic
    # merged.to_gromacs(...)
    
    return {"status": "deferred"}

if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8001)
