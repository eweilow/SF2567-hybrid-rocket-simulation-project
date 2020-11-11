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