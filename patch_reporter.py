with open('atomipy-web-module/src/components/VisualBuilder.tsx', 'r') as f:
    content = f.read()

# Replace the standard PDBReporter with the DynamicBoxPDBReporter
patch = """    # Trajectory output
          if simType === "npt" or simType === "nvt":
            pythonCode += `    class DynamicBoxPDBReporter:\\n`
            pythonCode += `        def __init__(self, file, reportInterval):\\n`
            pythonCode += `            self._out = open(file, 'w')\\n`
            pythonCode += `            self._reportInterval = reportInterval\\n`
            pythonCode += `            self._model = 1\\n`
            pythonCode += `            self._topology = None\\n`
            pythonCode += `        def describeNextReport(self, simulation):\\n`
            pythonCode += `            steps = self._reportInterval - simulation.currentStep % self._reportInterval\\n`
            pythonCode += `            return (steps, True, False, False, False, True)\\n`
            pythonCode += `        def report(self, simulation, state):\\n`
            pythonCode += `            if self._topology is None:\\n`
            pythonCode += `                self._topology = simulation.topology\\n`
            pythonCode += `                PDBFile.writeHeader(self._topology, self._out)\\n`
            pythonCode += `            self._topology.setPeriodicBoxVectors(state.getPeriodicBoxVectors())\\n`
            pythonCode += `            import io\\n`
            pythonCode += `            _tmp = io.StringIO()\\n`
            pythonCode += `            PDBFile.writeHeader(self._topology, _tmp)\\n`
            pythonCode += `            _lines = _tmp.getvalue().split('\\n')\\n`
            pythonCode += `            _cryst1 = next(l for l in _lines if l.startswith('CRYST1'))\\n`
            pythonCode += `            self._out.write(_cryst1 + '\\n')\\n`
            pythonCode += `            PDBFile.writeModel(self._topology, state.getPositions(), self._out, self._model)\\n`
            pythonCode += `            self._model += 1\\n`
            pythonCode += `        def __del__(self):\\n`
            pythonCode += `            self._out.close()\\n\\n`
            pythonCode += `    _omm_simulation.reporters.append(DynamicBoxPDBReporter('traj_${simType}.pdb', ${pdbFreq}))\\n`
          else:
            pythonCode += `    _omm_simulation.reporters.append(PDBReporter('traj_${simType}.pdb', ${pdbFreq}))\\n`
"""

# The code in VisualBuilder right now is:
#        if (writePdb) {
#          pythonCode += `    # Trajectory output\n`;
#          if (simType === "minimize") {
#            pythonCode += `    _omm_simulation.reporters.append(PDBReporter('traj_${simType}.pdb', 1))\n`;
#          } else {
#            pythonCode += `    _omm_simulation.reporters.append(PDBReporter('traj_${simType}.pdb', ${pdbFreq}))\n`;
#          }
#        }

old_code = """        if (writePdb) {
          pythonCode += `    # Trajectory output\\n`;
          if (simType === "minimize") {
            pythonCode += `    _omm_simulation.reporters.append(PDBReporter('traj_${simType}.pdb', 1))\\n`;
          } else {
            pythonCode += `    _omm_simulation.reporters.append(PDBReporter('traj_${simType}.pdb', ${pdbFreq}))\\n`;
          }
        }"""

new_code = """        if (writePdb) {
""" + patch + """        }"""

if old_code in content:
    content = content.replace(old_code, new_code)
    with open('atomipy-web-module/src/components/VisualBuilder.tsx', 'w') as f:
        f.write(content)
    print("Successfully patched VisualBuilder for DynamicBoxPDBReporter")
else:
    print("Could not find the target code to replace!")
