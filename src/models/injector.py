import math
from models.base import Model
from utils import constants
from equations.falloffs import combustionEfficiencyTransient, outletPhaseFalloff, inletPhaseFalloff, injectorTransientFalloff

injectorDischargeCoefficient = 0.83
injectorHoles = 38
injectorHoleRadius = constants.Lengths.mm * 0.75

def injectorFlow(Cd, density, diameter, pressureDrop):
  area = math.pow(diameter/2.0, 2) * math.pi

  if pressureDrop < 0:
    return 0
  
  return Cd * area * math.sqrt(2 * density * pressureDrop)

class InjectorModel(Model):
  def derivativesDependsOn(self, models):
    return []

  def derivedVariablesDependsOn(self, models):
    return [models["tank"]]
  
  derived_massFlow = 0

  def initializeState(self):
    return [0]

  def computeDerivatives(self, t, state, derived, models):
    return [0]

  def computeDerivedVariables(self, t, state, models):
    from models.tank import TankModel
    from models.nozzle import NozzleModel
    from models.injector import InjectorModel
    from models.combustion import CombustionModel
    
    tankPressure = models["tank"]["derived"][TankModel.derived_pressure]
    oxidizerDensity = models["tank"]["derived"][TankModel.derived_outletDensity]
    combustionChamberPressure = models["combustion"]["state"][CombustionModel.states_pressure]

    pressureDifferential = tankPressure - combustionChamberPressure

    startupTransient = injectorTransientFalloff(t)
    massFlow = startupTransient * injectorHoles * injectorFlow(injectorDischargeCoefficient, oxidizerDensity, injectorHoleRadius*2, pressureDifferential)

    return [massFlow]