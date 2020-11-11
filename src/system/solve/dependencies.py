import numpy as np


"""
Compute derived variables for all the models in the order that is defined by each model
in the function "derivedVariablesDependsOn".
"""
def recurseModelDependencies(t, obj, currentDepth, models, visited):
  if obj in visited:
    return visited

  visited.append(obj)
  dependsOn = obj["model"].derivedVariablesDependsOn(models)

  for dependency in dependsOn:
    recurseModelDependencies(t, dependency, currentDepth + 1, models, visited)
  
  obj["derived"] = obj["model"].computeDerivedVariables(t, obj["state"], models)

  return visited


from solve.states import applyModelStates, collectModelStates

def applyDerivedVariablesToResult(t, y, models):
  
  for key in models:
    models[key]["derivedResult"] = np.zeros((len(models[key]["derived"]), len(t)))

    # Recover derived variables
  for i in range(len(t)):
    print(t[i])
    applyModelStates(models, y[:,i])

    # Set up derived variables
    visited = []
    for key in models:
      visited = recurseModelDependencies(t[i], models[key], 0, models, visited)

    for key in models:
      models[key]["derivedResult"][:, i] = models[key]["derived"]
      
  applyModelStates(models, y)