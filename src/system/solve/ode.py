import numpy as np

from models.tank import TankModel
from models.injector import InjectorModel
from models.combustion import CombustionModel
from models.nozzle import NozzleModel

from solve.dependencies import recurseModelDependencies
from solve.models import initializeModel, initializeModelMatrices
from solve.states import applyModelStates, collectModelStates

def makeODE():
  models = {
    "tank": initializeModel("tank", TankModel()),
    "injector": initializeModel("injector", InjectorModel()),
    "combustion": initializeModel("combustion", CombustionModel()),
    "nozzle": initializeModel("nozzle", NozzleModel()),
  }

  fullSystemLength = initializeModelMatrices(models)

  def system(t, y):
    print(t)
    applyModelStates(models, y)

    # Set up derived variables
    visited = []
    for key in models:
      visited = recurseModelDependencies(t, models[key], 0, models, visited)

    fullSystem = np.zeros((1, fullSystemLength))
    for key in models:
      derivatives = np.array(models[key]["model"].computeDerivatives(t, models[key]["state"], models[key]["derived"], models))
      fullSystem = fullSystem + np.dot(models[key]["matrix"], derivatives)

    return fullSystem


  def no_oxidizer_mass(t, y): 
    models["tank"]["state"] = np.dot(models["tank"]["invMatrix"], y)
    return models["tank"]["state"][TankModel.states_oxidizerMass] - 0.25
  no_oxidizer_mass.terminal = True

  def no_oxidizer_energy(t, y): 
    models["tank"]["state"] = np.dot(models["tank"]["invMatrix"], y)
    return models["tank"]["state"][TankModel.states_totalEnergy]
  no_oxidizer_energy.terminal = True

  def no_chamber_pressure(t, y): 
    models["combustion"]["state"] = np.dot(models["combustion"]["invMatrix"], y)
    return models["combustion"]["state"][CombustionModel.states_pressure] - 200000

  no_chamber_pressure.terminal = True
  no_chamber_pressure.direction = -1

  events = (no_oxidizer_mass, no_oxidizer_energy, no_chamber_pressure)

  initialState = collectModelStates(models, fullSystemLength)

  return models, system, events, initialState