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

import equations.tankHeatTransfer as tankTransfer


hasInited = False
def init():
  global hasInited
  global liquidDensityInterpolator, gasDensityInterpolator, liquidInternalEnergyInterpolator, gasInternalEnergyInterpolator
  global liquidViscosityInterpolator, gasViscosityInterpolator, liquidThermalConductivityInterpolator, gasThermalConductivityInterpolator
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
  
  # https://github.com/aesirkth/mjollnir-propulsion-simulations/blob/master/combustion/properties/oxidizerProperties_create_interpolationtables.m

  propertiesInterpolant = "cubic"

  T_mu = np.concatenate(([182.33, 184.69], np.arange(185, 305, 5)))
  mu_l = 1e-3*np.array([0.4619,0.4306,0.4267,0.3705,0.3244,0.2861,0.2542,0.2272,0.2042,0.1846,0.1676,0.1528,0.1399,0.1284,0.1183,0.1093,0.1012,0.0939,0.0872,0.0810,0.0754,0.0700,0.0650,0.0601,0.0552,0.0501])
  mu_g = 1e-6*np.array([9.4,9.6,9.6,9.8,10.1,10.3,10.6,10.9,11.1,11.4,11.7,12.0,12.3,12.6,12.9,13.3,13.7,14.0,14.4,14.9,15.4,15.9,16.5,17.2,18.1,19.2])
  liquidViscosityInterpolator = scipy.interpolate.interp1d(T_mu, mu_l, kind=propertiesInterpolant, fill_value="extrapolate")
  gasViscosityInterpolator = scipy.interpolate.interp1d(T_mu, mu_g, kind=propertiesInterpolant, fill_value="extrapolate")

  T_k = np.concatenate(([182.33, 184.69], np.arange(185, 285, 5)))
  k_l = 1e-3*np.array([146.9,145.6,145.5,142.8,140.2,137.6,135.1,132.6,130.0,127.6,125.1,122.7,120.3,117.9,115.6,113.3,111.0,108.7,106.5,104.4,102.2,100.1])
  k_g = 1e-3*np.array([8.2,8.4,8.4,8.9,9.3,9.7,10.2,10.7,11.2,11.7,12.2,12.8,13.4,14.0,14.7,15.3,16.1,16.8,17.6,18.5,19.5,20.6])
  liquidThermalConductivityInterpolator = scipy.interpolate.interp1d(T_k, k_l, kind=propertiesInterpolant, fill_value="extrapolate")
  gasThermalConductivityInterpolator = scipy.interpolate.interp1d(T_k, k_g, kind=propertiesInterpolant, fill_value="extrapolate")

class EquilibriumTankModel(Model):
  def derivativesDependsOn(self, models):
    return [models["injector"], models["passiveVent"]]

  def derivedVariablesDependsOn(self, models):
    return [models["environment"]]

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
  states_gasWallTankTemperature = 2
  states_liquidWallTankTemperature = 3

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
    mDot = extendDerivative(timeHistory, stateHistory[self.states_oxidizerMass])
    pressure = extend(timeHistory, derivedVariablesHistory[self.derived_pressure])
    args = (mDot, pressure, )
    mask = [True, None]
    return mask, args

  def computeSimplifiedState(self, args, time):
    mDot, pressure = args
    return [mDot(time), 0, 0, 0]

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
    return [totalMass, totalEnergy, assumptions.tankInitialWallTemperature.get(), assumptions.tankInitialWallTemperature.get()]

  def computeDerivatives(self, t, state, derived, models):
    from models.tank import TankModel
    from models.nozzle import NozzleModel
    from models.injector import InjectorModel
    from models.combustion import CombustionModel
    from models.passiveVent import PassiveVentModel
    from models.environment import EnvironmentModel
    from models.flight import FlightModel

    massFlowTop = models["passiveVent"]["derived"][PassiveVentModel.derived_massFlow]
    massFlowOutlet = models["injector"]["derived"][InjectorModel.derived_massFlow]

    hOutlet = derived[self.derived_hOutlet]
    hTop = derived[self.derived_hTop]

    massFlow = -massFlowTop - massFlowOutlet

    gasDensity = derived[self.derived_gasDensity]
    liquidDensity = derived[self.derived_liquidDensity]
    temperature = derived[self.derived_temperature]
    liquidLevel = derived[self.derived_liquidLevel]

    airThermalConductivity = models["environment"]["derived"][EnvironmentModel.derived_ambientThermalConductivity]
    airPressure = models["environment"]["derived"][EnvironmentModel.derived_ambientPressure]
    airTemperature = models["environment"]["derived"][EnvironmentModel.derived_ambientTemperature]
    airViscosity = models["environment"]["derived"][EnvironmentModel.derived_ambientViscosity]
    airDensity = models["environment"]["derived"][EnvironmentModel.derived_ambientDensity]

    airCp = CP.PropsSI('CPMASS','T',airTemperature,'P',airPressure,'air')

    gasThermalConductivity = gasThermalConductivityInterpolator(temperature)
    gasCp = CP.PropsSI('CPMASS','T',temperature,'Q',1,'N2O')
    gasViscosity = gasViscosityInterpolator(temperature)

    liquidThermalConductivity = liquidThermalConductivityInterpolator(temperature)
    liquidCp = CP.PropsSI('CPMASS','T',temperature,'Q',0,'N2O')
    liquidViscosity = liquidViscosityInterpolator(temperature)
    
    gasWallHeight = (1 - liquidLevel) * assumptions.tankLength.get()
    liquidWallHeight = liquidLevel * assumptions.tankLength.get()

    gasWallTemperature = state[self.states_gasWallTankTemperature]
    liquidWallTemperature = state[self.states_liquidWallTankTemperature]
    
    g = models["flight"]["derived"][FlightModel.derived_perceivedGravity]
    energyFlowIntoGasPhaseFromTank = tankTransfer.tankWallHeatTransfer(
      gasThermalConductivity,
      gasCp,
      gasViscosity,
      gasDensity,
      g,
      temperature,
      gasWallTemperature,
      gasWallHeight,
      assumptions.tankInsideRadius.get(),
      2/5,
      0.021
    )

    energyFlowIntoGasPartOfTankFromAmbient = tankTransfer.tankWallHeatTransfer(
      airThermalConductivity,
      airCp,
      airViscosity,
      airDensity,
      g,
      gasWallTemperature,
      airTemperature,
      gasWallHeight,
      assumptions.tankInsideRadius.get(),
      1/4,
      0.59
    )

    energyFlowIntoLiquidPhaseFromTank = tankTransfer.tankWallHeatTransfer(
      liquidThermalConductivity,
      liquidCp,
      liquidViscosity,
      liquidDensity,
      g,
      temperature,
      liquidWallTemperature,
      liquidWallHeight,
      assumptions.tankInsideRadius.get(),
      2/5,
      0.021
    )

    energyFlowIntoLiquidPartOfTankFromAmbient = tankTransfer.tankWallHeatTransfer(
      airThermalConductivity,
      airCp,
      airViscosity,
      airDensity,
      g,
      liquidWallTemperature,
      airTemperature,
      liquidWallHeight,
      assumptions.tankInsideRadius.get(),
      1/4,
      0.59
    )

    energyFlowIntoGasPartOfTankFromLiquidPartOfTank = tankTransfer.tankWallHeatConduction(
      assumptions.tankWallThermalConductivity.get(),
      assumptions.tankInsideRadius.get(),
      assumptions.tankInsideRadius.get() + assumptions.tankThickness.get(),
      gasWallTemperature,
      liquidWallTemperature,
      assumptions.tankLength.get(),
      liquidWallHeight
    )

    massPerLength = tankTransfer.massPerLength(
      assumptions.tankInsideRadius.get(),
      assumptions.tankInsideRadius.get() + assumptions.tankThickness.get(),
      assumptions.tankWallDensity.get()
    )

    gasPartOfTankMass = gasWallHeight * massPerLength
    liquidPartOfTankMass = liquidWallHeight * massPerLength

    # Replace with DAE
    liquidLevelChangeDerivative = tankTransfer.estimate_liquidLevelChangeDerivative(
      -massFlow,
      liquidDensity,
      assumptions.tankInsideRadius.get()
    )

    boundaryTransferFromLiquidToGasPartOfTank = tankTransfer.tankWallControlVolumeBoundaryTransfer(
      assumptions.tankInsideRadius.get(),
      assumptions.tankThickness.get() + assumptions.tankInsideRadius.get(),
      assumptions.tankWallDensity.get(),
      liquidLevelChangeDerivative,
      assumptions.tankWallSpecificHeatCapacity.get(),
      gasWallTemperature,
      liquidWallTemperature
    )

    tankHeatCapacity = assumptions.tankWallSpecificHeatCapacity.get()

    dGasWallTemperature_dt = 1 / (gasPartOfTankMass * tankHeatCapacity) * (-energyFlowIntoGasPhaseFromTank + energyFlowIntoGasPartOfTankFromAmbient + energyFlowIntoGasPartOfTankFromLiquidPartOfTank - boundaryTransferFromLiquidToGasPartOfTank)
    dLiquidWallTemperature_dt = 1 / (liquidPartOfTankMass * tankHeatCapacity) * (-energyFlowIntoLiquidPhaseFromTank + energyFlowIntoLiquidPartOfTankFromAmbient - energyFlowIntoGasPartOfTankFromLiquidPartOfTank + boundaryTransferFromLiquidToGasPartOfTank)

    energyFlow = -massFlowOutlet * hOutlet - massFlowTop * hTop + energyFlowIntoGasPhaseFromTank + energyFlowIntoLiquidPhaseFromTank
    
    printMessages = False
    if printMessages:
      print("")
      print("t = {:.2f} s".format(t))
      print("masses:")
      print(" {:45s} = {:8.2f} kg".format("gas", derived[self.derived_gasMass]))
      print(" {:45s} = {:8.2f} kg".format("liquid", derived[self.derived_liquidMass]))
      print(" {:45s} = {:8.2f} kg".format("tank (gas part)", gasPartOfTankMass))
      print(" {:45s} = {:8.2f} kg".format("tank (liquid part)", liquidPartOfTankMass))
      print("energy:")
      print(" {:45s} = {:8.2f} kJ".format("oxidizer (total)", state[self.states_totalEnergy] / 1e3))
      print(" {:45s} = {:8.2f} kJ".format("tank (gas part)", state[self.states_gasWallTankTemperature] * gasPartOfTankMass * tankHeatCapacity / 1e3))
      print(" {:45s} = {:8.2f} kJ".format("tank (liquid part)", state[self.states_liquidWallTankTemperature] *  liquidPartOfTankMass * tankHeatCapacity / 1e3))
      print("heat capacity:")
      print(" {:45s} = {:8.2f} kJ/K".format("gas", derived[self.derived_gasMass] * gasCp / 1e3))
      print(" {:45s} = {:8.2f} kJ/K".format("liquid", derived[self.derived_liquidMass] * liquidCp / 1e3))
      print(" {:45s} = {:8.2f} kJ/K".format("tank (gas part)", gasPartOfTankMass * tankHeatCapacity / 1e3))
      print(" {:45s} = {:8.2f} kJ/K".format("tank (liquid part)", liquidPartOfTankMass * tankHeatCapacity / 1e3))
      print("heat flux:")
      print(" {:45s} = {:8.2f} kW".format("tank (gas part) -> gas", energyFlowIntoGasPhaseFromTank / 1e3))
      print(" {:45s} = {:8.2f} kW".format("ambient -> tank (gas part)", energyFlowIntoGasPartOfTankFromAmbient / 1e3))
      print(" {:45s} = {:8.2f} kW".format("tank (liquid part) -> liquid", energyFlowIntoLiquidPhaseFromTank / 1e3))
      print(" {:45s} = {:8.2f} kW".format("ambient -> tank (liquid part)", energyFlowIntoLiquidPartOfTankFromAmbient / 1e3))
      print(" {:45s} = {:8.2f} kW".format("tank (liquid part) -> tank (gas part)", energyFlowIntoGasPartOfTankFromLiquidPartOfTank / 1e3))
      print(" {:45s} = {:8.2f} kW".format("energy leaving oxidizer through top",  -massFlowTop * hTop / 1e3))
      print(" {:45s} = {:8.2f} kW".format("energy leaving oxidizer through bottom", -massFlowOutlet * hOutlet / 1e3))


    return [massFlow, energyFlow, dGasWallTemperature_dt, dLiquidWallTemperature_dt]

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

    # dLiquidWallTemperature_dt = tankWallHeatTransfer()
    return [
      pressure, densityOutlet, densityTop, temperature, hOutlet, hTop, liquidMass, gasMass, vaporQuality, liquidLevel, gasDensity, liquidDensity, gasVolume, liquidVolume, outletPhase, topPhase
    ]
  
class TankModel(EquilibriumTankModel):
  def __init__(self):
    super().__init__()
