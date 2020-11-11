import math
from models.base import Model
from utils import constants
from equations.falloffs import combustionEfficiencyTransient, outletPhaseFalloff, inletPhaseFalloff, injectorTransientFalloff
import CoolProp.CoolProp as CP

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
      
    tankEnthalpy =  models["tank"]["derived"][TankModel.derived_hOutlet]
    tankTemperature = models["tank"]["derived"][TankModel.derived_temperature]
    tankPhase = models["tank"]["derived"][TankModel.derived_outletPhase]
    
    try:
      tankEntropy = CP.PropsSI('S','P',tankPressure,'Q',tankPhase,'N2O')
      combustionChamberEnthalpy = CP.PropsSI('H','P',combustionChamberPressure,'S',tankEntropy,'N2O')
      combustionChamberDensity = CP.PropsSI('D','P',combustionChamberPressure,'S',tankEntropy,'N2O')

      saturatedTankPressure = CP.PropsSI('P','T',tankTemperature,'Q',0,'N2O')
      kappa = math.sqrt((tankPressure - combustionChamberPressure) / (saturatedTankPressure - combustionChamberPressure))

      # HEM injector model https://web.stanford.edu/~cantwell/Recent_publications/Zimmerman_et_al_AIAA_2013-4045.pdf
      G_spi = injectorDischargeCoefficient * math.sqrt(2 * oxidizerDensity * (tankPressure - combustionChamberPressure)) if tankPressure - combustionChamberPressure > 0 else 0
      G_hem = injectorDischargeCoefficient * combustionChamberDensity * math.sqrt(2 * (tankEnthalpy - combustionChamberEnthalpy)) if tankEnthalpy - combustionChamberEnthalpy > 0 else 0

      G = (kappa * G_spi + G_hem) / (1 + kappa)
      area = injectorHoles * math.pow(injectorHoleRadius, 2) * math.pi

      startupTransient = injectorTransientFalloff(t)
      massFlow = startupTransient * G * area
    

      return [massFlow]
    except:
      print("Failed to evaluate tank mass flow in time {:.4f}".format(t))
      return [0]