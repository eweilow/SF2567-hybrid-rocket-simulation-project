import numpy as np
from solve.matrices import makeSystemMatrices

def initializeModel(name, model):
  state = np.array(model.initializeState())
  return {
    "name": name,
    "model": model,
    "state": state,
    "derived": [],
    "derivedResult": [],
    "matrix": np.array([]),
    "invMatrix": np.array([])
  }

def computeFullSystemLength(models):
  fullSystemLength = 0
  for key in models:
    fullSystemLength = fullSystemLength + np.size(models[key]["state"])

  return fullSystemLength

def initializeModelMatrices(models):
  fullSystemLength = computeFullSystemLength(models)
  
  nextIndex = 0
  for key in models:
    matrix, invMatrix, nextIndex = makeSystemMatrices(nextIndex, np.size(models[key]["state"]), fullSystemLength)
    models[key]["matrix"] = matrix
    models[key]["invMatrix"] = invMatrix

  return fullSystemLength