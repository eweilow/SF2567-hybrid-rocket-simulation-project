from models.base import Model
import CoolProp.CoolProp as CP
import scipy.optimize
from utils import constants
from equations.falloffs import combustionEfficiencyTransient, outletPhaseFalloff, inletPhaseFalloff, injectorTransientFalloff
from equations.hemInjector import computeHEMInjector

import assumptions

class EquilibriumTankModel(Model):
  def derivativesDependsOn(self, models):
    return [models["injector"], models["passiveVent"]]

  def derivedVariablesDependsOn(self, models):
    return []

  def deriveVaporQuality(self, totalMass, totalEnergy, temperature):
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

  def initializeState(self):
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
      temperature = scipy.optimize.brentq(tempFn, 183, 309)
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

    gasDensity = CP.PropsSI('D','T',temperature,'Q',1,'N2O')
    liquidDensity = CP.PropsSI('D','T',temperature,'Q',0,'N2O')

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
