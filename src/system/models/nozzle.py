import math
from models.base import Model
from utils import constants
from equations.falloffs import combustionEfficiencyTransient, outletPhaseFalloff, inletPhaseFalloff, injectorTransientFalloff

ambientPressure = 101300

initialThroatRadius = constants.Lengths.mm * 19.466
exhaustRadius = constants.Lengths.mm * 44.777

nozzleEfficiency = 0.9
class NozzleModel(Model):
  def derivativesDependsOn(self, models):
    return []

  def derivedVariablesDependsOn(self, models):
    return [models["combustion"]]

  states_throatRadius = 0
  states_exhaustRadius = 1
  derived_massFlow = 0
  derived_thrust = 1
  derived_specificImpulse = 2
  
  def initializeState(self):
    return [initialThroatRadius, exhaustRadius]

  def computeDerivatives(self, t, state, derived, models):
    return [0, 0]

  def computeDerivedVariables(self, t, state, models):
    from models.tank import TankModel
    from models.nozzle import NozzleModel
    from models.injector import InjectorModel
    from models.combustion import CombustionModel

    chamberPressure = models["combustion"]["state"][CombustionModel.states_pressure]
    cStar = models["combustion"]["derived"][CombustionModel.derived_cStar]
    thrustCoefficient = models["combustion"]["derived"][CombustionModel.derived_thrustCoefficient] * nozzleEfficiency
    exhaustPressure = models["combustion"]["derived"][CombustionModel.derived_exhaustPressure]

    throatRadius = state[self.states_throatRadius]
    throatArea = math.pow(throatRadius, 2) * math.pi

    # err what
    nozzleMassFlow = (chamberPressure - ambientPressure) * throatArea / cStar
    if nozzleMassFlow < 0:
      nozzleMassFlow = 0
    # print(nozzleMassFlow, chamberTemperature, chamberPressure)
    # nozzleMassFlow = chamberPressure * throatArea / math.sqrt(chamberTemperature) * math.sqrt(chamberGamma / R * math.pow(2 / (chamberGamma + 1), (chamberGamma+1)/(chamberGamma-1)))
    
    exhaustArea = pow(models["nozzle"]["state"][NozzleModel.states_exhaustRadius], 2) * math.pi

    exhaustThrust = (exhaustPressure - ambientPressure) * exhaustArea if exhaustPressure > ambientPressure else 0

    thrust = cStar * thrustCoefficient * nozzleMassFlow + exhaustThrust

    specificImpulse = cStar * thrustCoefficient / 9.80665
    
    return [nozzleMassFlow, thrust, specificImpulse]

