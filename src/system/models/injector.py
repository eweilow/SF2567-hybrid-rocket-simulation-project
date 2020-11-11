import math
from models.base import Model
from utils import constants
from equations.falloffs import combustionEfficiencyTransient, outletPhaseFalloff, inletPhaseFalloff, injectorTransientFalloff
import CoolProp.CoolProp as CP

import assumptions

def injectorFlow(Cd, density, diameter, pressureDrop):
  area = math.pow(diameter/2.0, 2) * math.pi

  if pressureDrop < 0:
    return 0
  
  return Cd * area * math.sqrt(2 * density * pressureDrop)

def computeHEMInjector(
  dischargeCoefficient,
  numberOfHoles,
  holeDiameter,
  densityBefore,
  pressureBefore,
  pressureAfter,
  temperatureBefore,
  phaseBefore,
  enthalpyBefore
):
  try:
    entropyBefore = CP.PropsSI('S','T',temperatureBefore,'Q',phaseBefore,'N2O')
  except Exception as exc:
    print("Failed to evaluate entropyBefore", exc)
    return 0

  try:
    enthalpyAfter = CP.PropsSI('H','P',pressureAfter,'S',entropyBefore,'N2O')
  except Exception as exc:
    print("Failed to evaluate entropyBefore", exc)
    return 0

  try:
    densityAfter = CP.PropsSI('D','P',pressureAfter,'S',entropyBefore,'N2O')
  except Exception as exc:
    print("Failed to evaluate entropyBefore", exc)
    return 0

  try:
    saturatedPressureBefore = CP.PropsSI('P','T',temperatureBefore,'Q',0,'N2O')
  except Exception as exc:
    print("Failed to evaluate entropyBefore", exc)
    return 0

  # HEM injector model https://web.stanford.edu/~cantwell/Recent_publications/Zimmerman_et_al_AIAA_2013-4045.pdf
  kappa = math.sqrt((pressureBefore - pressureAfter) / (saturatedPressureBefore - pressureAfter))
  G_spi = dischargeCoefficient * math.sqrt(2 * densityBefore * (pressureBefore - pressureAfter)) if pressureBefore > pressureAfter else 0
  G_hem = dischargeCoefficient * densityAfter * math.sqrt(2 * (enthalpyBefore - enthalpyAfter)) if enthalpyBefore > enthalpyAfter else 0

  G = (kappa * G_spi + G_hem) / (1 + kappa)
  area = numberOfHoles * math.pow(holeDiameter / 2.0, 2) * math.pi

  return G * area



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