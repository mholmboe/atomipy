import re
with open('atomipy-web-module/src/components/nodes/SimulateNode.tsx', 'r') as f:
    content = f.read()

# Hide Energy Minimization section if not minimize
content = content.replace(
    """<label className="text-xs font-semibold text-muted-foreground block mb-1">Energy Minimization</label>""",
    """{simType === "minimize" && (
          <>
            <label className="text-xs font-semibold text-muted-foreground block mb-1">Energy Minimization</label>"""
)

content = content.replace(
    """onChange={(e) => updateNodeData(id, { ...data, miniSteps: Number(e.target.value) })}
              />
            </div>
          </div>
        </div>

        {/* MD Steps */}`,
    """onChange={(e) => updateNodeData(id, { ...data, miniSteps: Number(e.target.value) })}
              />
            </div>
          </>
        )}

        {/* MD Steps */}`
)
with open('atomipy-web-module/src/components/nodes/SimulateNode.tsx', 'w') as f:
    f.write(content)

print("Patched SimulateNode")
