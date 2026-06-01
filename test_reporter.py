import openmm.app as app
import openmm as mm
from openmm import unit
import io

# Create a mock topology
top = app.Topology()
top.setPeriodicBoxVectors((mm.Vec3(1,0,0), mm.Vec3(0,2,0), mm.Vec3(0,0,3))*unit.nanometers)

class MockState:
    def __init__(self, vectors):
        self.vectors = vectors
    def getPeriodicBoxVectors(self):
        return self.vectors

state1 = MockState((mm.Vec3(1,0,0), mm.Vec3(0,2,0), mm.Vec3(0,0,3))*unit.nanometers)
state2 = MockState((mm.Vec3(2,0,0), mm.Vec3(0,3,0), mm.Vec3(0,0,4))*unit.nanometers)

def get_cryst1(state):
    top.setPeriodicBoxVectors(state.getPeriodicBoxVectors())
    _tmp = io.StringIO()
    app.PDBFile.writeHeader(top, _tmp)
    _lines = _tmp.getvalue().split('\n')
    _cryst1 = next(l for l in _lines if l.startswith('CRYST1'))
    return _cryst1

print("State 1:", get_cryst1(state1))
print("State 2:", get_cryst1(state2))
