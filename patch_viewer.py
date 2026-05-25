import re
with open('atomipy-web-module/src/components/nodes/ViewerNode.tsx', 'r') as f:
    content = f.read()

# Add a pre-processing step for the PDB string to copy CRYST1 into every MODEL
patch = """
    if (pdb) {
      let processedPdb = pdb;
      if (isMulti && renderer === "3dmol") {
        // PDB files often only have one CRYST1 at the top. 3Dmol needs it per MODEL.
        const crystMatch = pdb.match(/^CRYST1.*$/m);
        if (crystMatch) {
          const crystStr = crystMatch[0];
          // Inject CRYST1 after each MODEL if not already there
          processedPdb = pdb.replace(/^(MODEL\\s+\\d+)\\r?\\n(?!CRYST1)/gm, `$1\\n${crystStr}\\n`);
        }
      }
      pdbLoadedRef.current = { renderer: "3dmol", pdb: processedPdb };
"""

# Replace the block
content = content.replace(
    """    if (pdb) {
      pdbLoadedRef.current = { renderer: "3dmol", pdb };""",
    patch
)

# And replace `pdb` with `processedPdb` in the model addition
content = content.replace(
    """const rawModels = viewer.addModelsAsFrames(pdb, "pdb", { keepH: true });""",
    """const rawModels = viewer.addModelsAsFrames(processedPdb, "pdb", { keepH: true });"""
)

content = content.replace(
    """model = viewer.addModel(pdb, "pdb", { keepH: true });""",
    """model = viewer.addModel(processedPdb, "pdb", { keepH: true });"""
)

with open('atomipy-web-module/src/components/nodes/ViewerNode.tsx', 'w') as f:
    f.write(content)

print("Patched ViewerNode.tsx")
