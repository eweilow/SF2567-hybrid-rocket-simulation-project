import numpy as np

def applyModelStates(models, y):
  for key in models:
    models[key]["state"] = np.dot(models[key]["invMatrix"], y)

def collectModelStates(models, fullSystemLength):
  fullSystem = np.zeros(fullSystemLength)
  for key in models:
    fullSystem = fullSystem + np.dot(models[key]["matrix"], models[key]["state"])

  return fullSystem
