import math
from models.base import Model
from utils import constants
from equations.falloffs import combustionEfficiencyTransient, outletPhaseFalloff, inletPhaseFalloff, injectorTransientFalloff
atmosphericPressure = 101300

fuelDensity = 900
initialPortRadius = constants.Lengths.mm * 25
aConstant = 0.155e-3 #http://www.scielo.org.za/pdf/rd/v33/05.pdf
nConstant = 0.5
portLength = 33 * constants.Lengths.cm

combustionEfficiency = 0.9

class CombustionModel(Model):
  def derivativesDependsOn(self, models):
    return [models["injector", models["nozzle"]]]

  def derivedVariablesDependsOn(self, models):
    return [models["injector"], models["tank"]]
    
  def __init__(self):
    from utils.cea import NasaCEA
    self.cea = NasaCEA()

  states_pressure = 0
  states_portRadius = 1
  states_burntFuel = 2

  derived_volume = 0
  derived_temperature = 1
  derived_gamma = 2
  derived_portArea = 3
  derived_burningArea = 4
  derived_CpT = 5
  derived_rDot = 6
  derived_fuelFlow = 7
  derived_ofRatio = 8
  derived_cStar = 9
  derived_thrustCoefficient = 10
  derived_exhaustPressure = 11
  derived_molecularMass = 12
  derived_density = 12

  def initializeState(self):
    initialPressure = atmosphericPressure
    return [initialPressure, initialPortRadius, 0]

  def computeDerivatives(self, t, state, derived, models):
    from models.tank import TankModel
    from models.nozzle import NozzleModel
    from models.injector import InjectorModel
    from models.combustion import CombustionModel

    volume = derived[self.derived_volume]
    gamma = derived[self.derived_gamma]
    CpT = derived[self.derived_CpT]
    burningArea = derived[self.derived_burningArea]
    dPortRadius_dt = derived[self.derived_rDot]
    fuelMassFlow = derived[self.derived_fuelFlow]

    injectorMassFlow = models["injector"]["derived"][InjectorModel.derived_massFlow]
    nozzleMassFlow = models["nozzle"]["derived"][NozzleModel.derived_massFlow]

    dPressure_dt = (gamma - 1) / (volume) * CpT * (injectorMassFlow + fuelMassFlow - nozzleMassFlow) - gamma * state[self.states_pressure] / volume * (burningArea * dPortRadius_dt)

    return [dPressure_dt, dPortRadius_dt, fuelMassFlow]

  def computeDerivedVariables(self, t, state, models):
    from models.tank import TankModel
    from models.nozzle import NozzleModel
    from models.injector import InjectorModel
    from models.combustion import CombustionModel

    portRadius = state[self.states_portRadius]
    portArea = math.pow(portRadius, 2) * math.pi

    burningArea = portRadius * 2 * portLength * math.pi
   
    injectorMassFlow = models["injector"]["derived"][InjectorModel.derived_massFlow]
    oxidizerFlux = injectorMassFlow / portArea

    dPortRadius_dt = aConstant * math.pow(oxidizerFlux, nConstant)
    fuelMassFlow = (math.pow(portRadius + dPortRadius_dt, 2) * math.pi - math.pow(portRadius, 2) * math.pi) * portLength * fuelDensity

    ofRatio = injectorMassFlow / fuelMassFlow if fuelMassFlow > 1e-3 else 1000

    throatArea = pow(models["nozzle"]["state"][NozzleModel.states_throatRadius], 2) * math.pi
    exhaustArea = pow(models["nozzle"]["state"][NozzleModel.states_exhaustRadius], 2) * math.pi

    areaRatio = exhaustArea / throatArea
    Isp, Cp, molecularMass, cStar, temperature, gamma, density, thrustCoefficient, exhaustPressure = self.cea.getPerformanceParameters(state[self.states_pressure], atmosphericPressure, models["tank"]["derived"][TankModel.derived_temperature], areaRatio, ofRatio)
    
    CpT = temperature * Cp
    cStar = cStar * combustionEfficiency
    startupTransient = combustionEfficiencyTransient(t)
    cStar = cStar * startupTransient

    volume = 1 * constants.Volume.liter + portArea * portLength

    return [volume, temperature, gamma, portArea, burningArea, CpT, dPortRadius_dt, fuelMassFlow, ofRatio, cStar, thrustCoefficient, exhaustPressure, molecularMass, density]
