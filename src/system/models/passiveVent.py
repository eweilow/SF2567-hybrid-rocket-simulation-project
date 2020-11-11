from models.base import Model
from utils import constants
from equations.falloffs import combustionEfficiencyTransient, outletPhaseFalloff, inletPhaseFalloff, injectorTransientFalloff
from equations.hemInjector import computeHEMInjector

import assumptions

class PassiveVentModel(Model):
  def derivativesDependsOn(self, models):
    return []

  def derivedVariablesDependsOn(self, models):
    return [models["tank"], models["environment"]]
  
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
    from models.environment import EnvironmentModel

    ambientPressure = models["environment"]["derived"][EnvironmentModel.derived_ambientPressure]
    
    tankPressure = models["tank"]["derived"][TankModel.derived_pressure]
    oxidizerDensity = models["tank"]["derived"][TankModel.derived_outletTopDensity]
      
    tankEnthalpy =  models["tank"]["derived"][TankModel.derived_hTop]
    tankTemperature = models["tank"]["derived"][TankModel.derived_temperature]
    tankPhase = models["tank"]["derived"][TankModel.derived_topPhase]
    
    startupTransient = injectorTransientFalloff(t)
    massFlow = startupTransient * computeHEMInjector(
      assumptions.tankPassiveVentDischargeCoefficient.get(), 
      1,
      assumptions.tankPassiveVentDiameter.get(),
      oxidizerDensity,
      tankPressure,
      ambientPressure,
      tankTemperature,
      tankPhase,
      tankEnthalpy
    )

    return [massFlow]