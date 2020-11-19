from models.base import Model
import CoolProp.CoolProp as CP
import scipy.optimize
import numpy as np
import scipy.interpolate
from utils import constants
from equations.falloffs import combustionEfficiencyTransient, outletPhaseFalloff, inletPhaseFalloff, injectorTransientFalloff
from equations.hemInjector import computeHEMInjector

import options
import assumptions

from utils.extend import extendDerivative, extend

hasInited = False
def init():
  global hasInited
  global liquidDensityInterpolator, gasDensityInterpolator, liquidInternalEnergyInterpolator, gasInternalEnergyInterpolator
  if options.enableTankTemperatureInterpolation and not hasInited:
    hasInited = True
    tankInterpolationValues = np.linspace(183, 309, options.tankInterpolantPointCount)
    
    interpolant = options.tankInterpolant

    print("computing liquidDensities")
    liquidDensities = CP.PropsSI('D','T',tankInterpolationValues,'Q',0,'N2O')
    liquidDensityInterpolator = scipy.interpolate.interp1d(tankInterpolationValues, liquidDensities, kind=interpolant, bounds_error=True, copy=False)

    print("computing gasDensities")
    gasDensities = CP.PropsSI('D','T',tankInterpolationValues,'Q',1,'N2O')
    gasDensityInterpolator = scipy.interpolate.interp1d(tankInterpolationValues, gasDensities, kind=interpolant, bounds_error=True, copy=False)

    print("computing liquidSpecificInternalEnergies")
    liquidSpecificInternalEnergies = CP.PropsSI('U','T',tankInterpolationValues,'Q',0,'N2O')
    liquidInternalEnergyInterpolator = scipy.interpolate.interp1d(tankInterpolationValues, liquidSpecificInternalEnergies, kind=interpolant, bounds_error=True, copy=False)
    
    print("computing gasSpecificInternalEnergies")
    gasSpecificInternalEnergies = CP.PropsSI('U','T',tankInterpolationValues,'Q',1,'N2O')
    gasInternalEnergyInterpolator = scipy.interpolate.interp1d(tankInterpolationValues, gasSpecificInternalEnergies, kind=interpolant, bounds_error=True, copy=False)

class EquilibriumTankModel(Model):
  def derivativesDependsOn(self, models):
    return [models["injector"], models["passiveVent"]]

  def derivedVariablesDependsOn(self, models):
    return []

  def deriveVaporQuality(self, totalMass, totalEnergy, temperature):
    if options.enableTankTemperatureInterpolation:
      liquidSpecificInternalEnergy = liquidInternalEnergyInterpolator(temperature)
      gasSpecificInternalEnergy = gasInternalEnergyInterpolator(temperature)
      
    else:
      liquidSpecificInternalEnergy = CP.PropsSI('U','T',temperature,'Q',0,'N2O')
      gasSpecificInternalEnergy = CP.PropsSI('U','T',temperature,'Q',1,'N2O')

    x = (totalEnergy / totalMass - liquidSpecificInternalEnergy) / (gasSpecificInternalEnergy - liquidSpecificInternalEnergy)

    return x

  def deriveTotalEnergy(self, gasMass, liquidMass, temperature):
    liquidSpecificInternalEnergy = CP.PropsSI('U','T',temperature,'Q',0,'N2O')
    gasSpecificInternalEnergy = CP.PropsSI('U','T',temperature,'Q',1,'N2O')
    
    liquidInternalEnergy = liquidSpecificInternalEnergy * liquidMass
    gasInternalEnergy = gasSpecificInternalEnergy * gasMass

    return liquidInternalEnergy + gasInternalEnergy

  def derivePropellantVolume(self, totalMass, totalEnergy, temperature):
    if options.enableTankTemperatureInterpolation:
      liquidDensity = liquidDensityInterpolator(temperature)
      gasDensity = gasDensityInterpolator(temperature)
    else:
      liquidDensity = CP.PropsSI('D','T',temperature,'Q',0,'N2O')
      gasDensity = CP.PropsSI('D','T',temperature,'Q',1,'N2O')

    x = self.deriveVaporQuality(totalMass, totalEnergy, temperature)
    volume = totalMass * ((1-x)/liquidDensity + x/gasDensity)
    return volume

  states_oxidizerMass = 0
  states_totalEnergy = 1
  derived_pressure = 0
  derived_outletDensity = 1
  derived_outletTopDensity = 2
  derived_temperature = 3
  derived_hOutlet = 4
  derived_hTop = 5
  derived_liquidMass = 6
  derived_gasMass = 7
  derived_vaporQuality = 8
  derived_liquidLevel = 9
  derived_gasDensity = 10
  derived_liquidDensity = 11
  derived_gasVolume = 12
  derived_liquidVolume = 13
  derived_outletPhase = 14
  derived_topPhase = 15

  def initializeSimplifiedModel(self, timeHistory, stateHistory, derivedVariablesHistory):
    print("mDot...")
    mDot = extendDerivative(timeHistory, stateHistory[self.states_oxidizerMass])
    pressure = extend(timeHistory, derivedVariablesHistory[self.derived_pressure])
    args = (mDot, pressure, )
    mask = [True, None]
    return mask, args

  def computeSimplifiedState(self, args, time):
    mDot, pressure = args
    return [mDot(time), 0]

  def computeSimplifiedDerivedVariables(self, args, time):
    mdot, pressure = args
    ret = [None for i in range(16)]
    ret[self.derived_pressure] = pressure(time)

    return ret

  def initializeState(self):
    init()

    filledGasDensity = CP.PropsSI('D','T',assumptions.tankFilledTemperature.get(),'Q',1,'N2O')
    filledLiquidDensity = CP.PropsSI('D','T',assumptions.tankFilledTemperature.get(),'Q',0,'N2O')

    filledLiquidVolume = assumptions.tankVolume.get() * assumptions.tankFillingGrade.get()
    filledGasVolume = assumptions.tankVolume.get() - filledLiquidVolume

    filledGasMass = filledGasDensity * filledGasVolume
    filledLiquidMass = filledLiquidDensity * filledLiquidVolume

    totalMass = filledGasMass + filledLiquidMass
    totalEnergy = self.deriveTotalEnergy(filledGasMass, filledLiquidMass, assumptions.tankFilledTemperature.get())
    return [totalMass, totalEnergy]

  def computeDerivatives(self, t, state, derived, models):
    from models.tank import TankModel
    from models.nozzle import NozzleModel
    from models.injector import InjectorModel
    from models.combustion import CombustionModel
    from models.passiveVent import PassiveVentModel

    massFlowTop = models["passiveVent"]["derived"][PassiveVentModel.derived_massFlow]
    massFlowOutlet = models["injector"]["derived"][InjectorModel.derived_massFlow]

    hOutlet = derived[self.derived_hOutlet]
    hTop = derived[self.derived_hTop]

    massFlow = massFlowTop + massFlowOutlet
    energyFlow = massFlowOutlet * hOutlet + massFlowTop * hTop
    return [-massFlow, -energyFlow]

  tankBurnoutTime = 0

  def computeDerivedVariables(self, t, state, models):
    totalMass = state[self.states_oxidizerMass]
    totalEnergy = state[self.states_totalEnergy]

    def tempFn(temperature):
      return self.derivePropellantVolume(totalMass, totalEnergy, temperature) - assumptions.tankVolume.get()

    try:
      if options.rootFindingType == "bisect":
        temperature = scipy.optimize.bisect(tempFn, 183, 309)

      if options.rootFindingType == "brentq":
        temperature = scipy.optimize.brentq(tempFn, 183, 309)

      if options.rootFindingType == "toms748":
        temperature = scipy.optimize.toms748(tempFn, 183, 309, k=2)

      if options.rootFindingType == "newton":
        temperature = scipy.optimize.newton(tempFn, 293)

    except:
      temperature = 183

    pressure = CP.PropsSI('P','T',temperature,'Q',0,'N2O')
    vaporQuality = self.deriveVaporQuality(totalMass, totalEnergy, temperature)

    if vaporQuality < 0:
      vaporQuality = 0
    if vaporQuality > 1:
      vaporQuality = 1

    gasMass = totalMass * vaporQuality
    liquidMass = totalMass - gasMass

    gasDensity = gasDensityInterpolator(temperature) if options.enableTankTemperatureInterpolation else CP.PropsSI('D','T',temperature,'Q',1,'N2O')
    liquidDensity = liquidDensityInterpolator(temperature) if options.enableTankTemperatureInterpolation else CP.PropsSI('D','T',temperature,'Q',0,'N2O')

    gasVolume = gasMass / liquidDensity
    liquidVolume = liquidMass / liquidDensity

    liquidLevel = liquidVolume / assumptions.tankVolume.get()

    # outlet is liquid unless vapor quality is 1
    # 0 if liquidLevel > bottomFalloff else (1 - max(0, liquidLevel / bottomFalloff))
    outletPhase = outletPhaseFalloff(liquidLevel)
    #top is gas if there is some vapor in the tank
    # 1 if liquidLevel < (1 - topFalloff) else min(1, (liquidLevel - (1 - topFalloff)) / topFalloff)
    topPhase =  inletPhaseFalloff(liquidLevel)

    densityOutlet = CP.PropsSI('D','T',temperature,'Q',outletPhase,'N2O')
    densityTop = CP.PropsSI('D','T',temperature,'Q',topPhase,'N2O')

    hOutlet = CP.PropsSI('H','T',temperature,'Q',outletPhase,'N2O')
    hTop = CP.PropsSI('H','T',temperature,'Q',topPhase,'N2O')

    return [pressure, densityOutlet, densityTop, temperature, hOutlet, hTop, liquidMass, gasMass, vaporQuality, liquidLevel, gasDensity, liquidDensity, gasVolume, liquidVolume, outletPhase, topPhase]
  
class TankModel(EquilibriumTankModel):
  def __init__(self):
    super().__init__()
