import CoolProp.CoolProp as CP
import numpy as np
import constants
import math
from scikits.odes import dae

def computeNormalInjector(
  dischargeCoefficient,
  numberOfHoles,
  holeDiameter,
  densityBefore,
  pressureBefore,
  pressureAfter
):
  area = numberOfHoles * math.pow(holeDiameter / 2.0, 2) * math.pi

  sign = 1 if pressureBefore > pressureAfter else -1
  G = dischargeCoefficient * math.sqrt(2 * densityBefore * abs(pressureBefore - pressureAfter))
  return sign * G * area



tankFilledTemperature = 293
tankFillingGrade = 0
tankVolume = 35 * constants.Volume.liter

filledGasDensity = CP.PropsSI('D','T',tankFilledTemperature,'Q',1,'N2O')
filledLiquidDensity = CP.PropsSI('D','T',tankFilledTemperature,'Q',0,'N2O')

filledLiquidVolume = tankVolume * tankFillingGrade
filledGasVolume = tankVolume - filledLiquidVolume

filledGasMass = filledGasDensity * filledGasVolume
filledLiquidMass = filledLiquidDensity * filledLiquidVolume

totalMass = filledGasMass + filledLiquidMass

liquidSpecificInternalEnergy = CP.PropsSI('U','T',tankFilledTemperature,'Q',0,'N2O')
gasSpecificInternalEnergy = CP.PropsSI('U','T',tankFilledTemperature,'Q',1,'N2O')

liquidInternalEnergy = liquidSpecificInternalEnergy * filledLiquidMass
gasInternalEnergy = gasSpecificInternalEnergy * filledGasMass

totalEnergy = liquidInternalEnergy + gasInternalEnergy

filledPressure = CP.PropsSI('P','T',tankFilledTemperature,'Q',0,'N2O')

initialState = [totalMass, totalEnergy, tankFilledTemperature, filledPressure, 1, 0, 0, totalMass, 0]
initialDerivatives = [0, 0, 0, 0, 0, 0, 0, 0, 0]


def computeResidual(t, y, ydot, residual):
  mass = y[0]
  energy = y[1]
  temperature = y[2]

  inlet_temperature = 293
  inlet_density = CP.PropsSI('D','T',inlet_temperature,'Q',0,'N2O')
  inlet_pressure = CP.PropsSI('P','T',inlet_temperature,'Q',0,'N2O')
  inlet_massFlow = computeNormalInjector(0.7, 1, 20 * constants.Lengths.mm, inlet_density, inlet_pressure, y[3])
  inlet_h = CP.PropsSI('H','T',temperature,'Q',0,'N2O')

  outlet_density = CP.PropsSI('D','T',temperature,'Q',1,'N2O')
  outlet_pressure = CP.PropsSI('P','T',temperature,'Q',1,'N2O')
  outlet_massFlow = computeNormalInjector(0.7, 1, 0.7 * constants.Lengths.mm, outlet_density, outlet_pressure, 101300)
  outlet_h = CP.PropsSI('H','T',temperature,'Q',1,'N2O')

  liquidSpecificInternalEnergy = CP.PropsSI('U','T',temperature,'Q',0,'N2O')
  gasSpecificInternalEnergy = CP.PropsSI('U','T',temperature,'Q',1,'N2O')

  vaporQuality = (energy / mass - liquidSpecificInternalEnergy) / (gasSpecificInternalEnergy - liquidSpecificInternalEnergy)
  liquidDensity = CP.PropsSI('D','T',temperature,'Q',0,'N2O')
  gasDensity = CP.PropsSI('D','T',temperature,'Q',1,'N2O')
  volume = mass * ((1-vaporQuality)/liquidDensity + vaporQuality/gasDensity)
  pressure = CP.PropsSI('P','T',temperature,'Q',0,'N2O')
  
  gasMass = mass * vaporQuality
  liquidMass = mass - gasMass
  
  liquidVolume = liquidMass / liquidDensity
  liquidLevel = liquidVolume / tankVolume

  massFlow = inlet_massFlow - outlet_massFlow
  energyFlow = inlet_massFlow * inlet_h - outlet_massFlow * outlet_h

  residual[0] = ydot[0] - massFlow
  residual[1] = ydot[1] - energyFlow
  residual[2] = tankVolume - volume
  residual[3] = y[3] - pressure
  residual[4] = y[4] - liquidLevel
  residual[5] = ydot[5] - inlet_massFlow
  residual[6] = ydot[6] - outlet_massFlow
  residual[7] = y[7] - gasMass
  residual[8] = y[8] - liquidMass

rtol = 1e-4
atol = 1e-4

algebraic_vars_idx = [2, 3, 4, 7, 8]

solver = dae('ida', computeResidual, 
    compute_initcond="yp0",
    exclude_algvar_from_error=False,
    algebraic_vars_idx=algebraic_vars_idx,
    atol=atol,
    rtol=rtol,
    max_steps=500000,
  )


times = np.linspace(0, 50, 100)
solution = solver.solve(times, initialState, initialDerivatives)

print(solution.message)

with open('/data/tankFilling.npy', 'wb') as f:
  t = np.array(solution.values.t)
  y = np.array(solution.values.y)
  np.save(f, t) 
  np.save(f, y) 