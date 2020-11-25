import math
from models.base import Model
from utils import constants
from equations.falloffs import combustionEfficiencyTransient, outletPhaseFalloff, inletPhaseFalloff, injectorTransientFalloff

from equations.falloffs import sigmoid
import assumptions

from utils.extend import extend

class NozzleModel(Model):
  def derivativesDependsOn(self, models):
    return []

  def derivedVariablesDependsOn(self, models):
    return [models["combustion"], models["environment"]]

  states_throatRadius = 0
  states_exhaustRadius = 1
  derived_massFlow = 0
  derived_thrust = 1
  derived_specificImpulse = 2
  
  def initializeSimplifiedModel(self, timeHistory, stateHistory, derivedVariablesHistory):
    thrust = extend(timeHistory, derivedVariablesHistory[self.derived_thrust])
    args = (thrust, )
    mask = [None, None]
    return mask, args

  def computeSimplifiedState(self, args, time):
    return [0, 0]

  def computeSimplifiedDerivedVariables(self, args, time):
    thrust, = args
    return [None, thrust(time), None]


  def initializeState(self):
    return [assumptions.nozzleThroatRadius.get(), assumptions.nozzleExhaustRadius.get()]

  def computeDerivatives(self, t, state, derived, models):
    erosionRate = assumptions.nozzleErosionConstant.get()
    erosionStart = assumptions.nozzleErosionStart.get()
    erosionStartRadius = assumptions.nozzleErosionStartRadius.get()
    erosionFalloff = sigmoid(t, erosionStart, erosionStartRadius)

    dThroatRadius_dt = erosionFalloff * erosionRate
    return [dThroatRadius_dt, 0]

  def computeDerivedVariables(self, t, state, models):
    from models.tank import TankModel
    from models.nozzle import NozzleModel
    from models.injector import InjectorModel
    from models.combustion import CombustionModel
    from models.environment import EnvironmentModel

    ambientPressure = models["environment"]["derived"][EnvironmentModel.derived_ambientPressure]
    chamberPressure = models["combustion"]["state"][CombustionModel.states_pressure]
    cStar = models["combustion"]["derived"][CombustionModel.derived_cStar]
    thrustCoefficient = models["combustion"]["derived"][CombustionModel.derived_thrustCoefficient] * assumptions.nozzleEfficiency.get()
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

