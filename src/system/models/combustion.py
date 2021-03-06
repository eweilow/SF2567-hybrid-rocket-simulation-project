import math
from models.base import Model
from utils import constants
from equations.falloffs import combustionEfficiencyTransient, outletPhaseFalloff, inletPhaseFalloff, injectorTransientFalloff

import assumptions

import options
from utils.extend import extend, extendDerivative

class CombustionModel(Model):
  def derivativesDependsOn(self, models):
    return [models["injector", models["nozzle"]]]

  def derivedVariablesDependsOn(self, models):
    return [models["injector"], models["tank"], models["environment"]]
    
  def __init__(self):
    from utils.cea import NasaCEA
    self.cea = NasaCEA()

  states_pressure = 0
  states_portRadius = 1
  states_fuelMass = 2

  states_dae_gamma = 3

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
  derived_density = 13
  derived_oxidizerDensity = 14

  def initializeSimplifiedModel(self, timeHistory, stateHistory, derivedVariablesHistory):
    rDot = extend(timeHistory, derivedVariablesHistory[self.derived_rDot])
    pDot = extendDerivative(timeHistory, stateHistory[self.states_pressure])
    fuelMassDot = extendDerivative(timeHistory, stateHistory[self.states_fuelMass])
    args = (rDot, pDot, fuelMassDot)
    mask = [True, True, True]
    return mask, args

  def computeSimplifiedState(self, args, time):
    rDot, pDot, fuelMassDot = args
    return [pDot(time), rDot(time), fuelMassDot(time), 0, 0]

  def computeSimplifiedDerivedVariables(self, args, time):
    return [None for i in range(15)]

  def initializeState(self):
    initialPressure = assumptions.initialAtmosphericPressure.get()

    fuelVolume = (math.pow(assumptions.fuelPortMaximumRadius.get(), 2) * math.pi - math.pow(assumptions.fuelPortInitialRadius.get(), 2) * math.pi) * assumptions.fuelPortLength.get()
    fuelMass = fuelVolume * assumptions.fuelDensitySolid.get()

    return [initialPressure, assumptions.fuelPortInitialRadius.get(), fuelMass, 0, 0]

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

    massFlow = injectorMassFlow + fuelMassFlow - nozzleMassFlow

    dPressure_dt = (gamma - 1) / (volume) * CpT * (massFlow)
    if not options.currentlySolvingWithDAE:
      dPressure_dt -= gamma * state[self.states_pressure] / volume * (burningArea * dPortRadius_dt)

    return [dPressure_dt, dPortRadius_dt, -fuelMassFlow, gamma, volume]

  def computeDerivedVariables(self, t, state, models):
    from models.tank import TankModel
    from models.nozzle import NozzleModel
    from models.injector import InjectorModel
    from models.combustion import CombustionModel
    from models.environment import EnvironmentModel

    ambientPressure = models["environment"]["derived"][EnvironmentModel.derived_ambientPressure]

    portRadius = state[self.states_portRadius]
    portArea = math.pow(portRadius, 2) * math.pi

    burningArea = portRadius * 2 * assumptions.fuelPortLength.get() * math.pi
   
    injectorMassFlow = models["injector"]["derived"][InjectorModel.derived_massFlow]
    oxidizerFlux = injectorMassFlow / portArea

    pressureFactor = math.sqrt(state[self.states_pressure] / assumptions.maximumRegressionRateAt.get())

    dPortRadius_dt = pressureFactor *assumptions.fuelGrainAConstant.get() * math.pow(oxidizerFlux, assumptions.fuelGrainNConstant.get())
    fuelMassFlow = (math.pow(portRadius + dPortRadius_dt, 2) * math.pi - math.pow(portRadius, 2) * math.pi) * assumptions.fuelPortLength.get() * assumptions.fuelDensitySolid.get()

    ofRatio = injectorMassFlow / fuelMassFlow if fuelMassFlow > 1e-3 else 1000

    throatArea = pow(models["nozzle"]["state"][NozzleModel.states_throatRadius], 2) * math.pi
    exhaustArea = pow(models["nozzle"]["state"][NozzleModel.states_exhaustRadius], 2) * math.pi

    areaRatio = exhaustArea / throatArea
    _, Cp, molecularMass, cStar, temperature, gamma, density, thrustCoefficient, exhaustPressure, oxidizerDensity = self.cea.getPerformanceParameters(state[self.states_pressure], ambientPressure, models["tank"]["derived"][TankModel.derived_temperature], areaRatio, ofRatio)
    
    CpT = temperature * Cp
    cStar *= assumptions.combustionEfficiency.get()
    startupTransient = combustionEfficiencyTransient(t)
    cStar *= startupTransient
    fuelMassFlow *= startupTransient

    volume = assumptions.preCombustionChamberVolume.get() + assumptions.postCombustionChamberVolume.get() + portArea * assumptions.fuelPortLength.get()

    return [volume, temperature, gamma, portArea, burningArea, CpT, dPortRadius_dt, fuelMassFlow, ofRatio, cStar, thrustCoefficient, exhaustPressure, molecularMass, density, oxidizerDensity]
