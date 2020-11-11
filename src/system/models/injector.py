from models.base import Model
from utils import constants
from equations.falloffs import combustionEfficiencyTransient, outletPhaseFalloff, inletPhaseFalloff, injectorTransientFalloff
from equations.hemInjector import computeHEMInjector

import assumptions

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
      
    tankEnthalpy =  models["tank"]["derived"][TankModel.derived_hOutlet]
    tankTemperature = models["tank"]["derived"][TankModel.derived_temperature]
    tankPhase = models["tank"]["derived"][TankModel.derived_outletPhase]
    
    startupTransient = injectorTransientFalloff(t)
    massFlow = startupTransient * computeHEMInjector(
      assumptions.injectorHoleDischargeCoefficient.get(),
      assumptions.injectorHoleCount.get(),
      assumptions.injectorHoleDiameter.get(),
      oxidizerDensity,
      tankPressure,
      combustionChamberPressure,
      tankTemperature,
      tankPhase,
      tankEnthalpy
    )

    return [massFlow]