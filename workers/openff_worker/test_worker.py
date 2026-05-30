import pytest
from fastapi.testclient import TestClient
from unittest.mock import patch, MagicMock
from worker import app

client = TestClient(app)

@patch("worker.Molecule")
@patch("worker.Interchange")
@patch("worker.ForceField")
def test_parametrize_endpoint(mock_ForceField, mock_Interchange, mock_Molecule):
    # Setup mocks
    mock_mol = MagicMock()
    mock_Molecule.from_smiles.return_value = mock_mol
    mock_mol.to_topology.return_value = "mock_topology"
    
    mock_interchange_inst = MagicMock()
    mock_Interchange.from_smirnoff.return_value = mock_interchange_inst
    
    # Send request
    response = client.post("/parametrize", json={
        "smiles": "CCO",
        "forcefield": "openff_sage"
    })
    
    # Assertions
    assert response.status_code == 200
    data = response.json()
    assert "job_id" in data
    assert "top" in data
    assert "gro" in data
    assert data["top"].endswith("organic.top")
    assert data["gro"].endswith("organic.gro")
    
    # Verify openff components were called
    mock_Molecule.from_smiles.assert_called_once_with("CCO", allow_undefined_stereo=True)
    mock_mol.generate_conformers.assert_called_once_with(n_conformers=1)
    mock_ForceField.assert_called_once_with("openff-2.0.0.offxml")
    mock_Interchange.from_smirnoff.assert_called_once()
    mock_interchange_inst.to_gromacs.assert_called_once()

def test_mix_endpoint_not_implemented():
    response = client.post("/mix", json={
        "organic_top": "org.top",
        "organic_gro": "org.gro",
        "mineral_top": "min.top",
        "mineral_gro": "min.gro",
        "box": [50.0, 50.0, 50.0]
    })
    # Feature flagged as not implemented
    assert response.status_code == 501
