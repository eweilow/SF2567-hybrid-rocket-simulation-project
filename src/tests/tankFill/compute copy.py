import CoolProp.CoolProp as CP
import numpy as np
import constants
from scikits.odes import dae

def computeNormalInjector(
  dischargeCoefficient,
  numberOfHoles,
  holeDiameter,
  densityBefore,
  pressureBefore,
  pressureAfter,
  temperatureBefore,
  phaseBefore,
  enthalpyBefore,
):
  area = numberOfHoles * math.pow(holeDiameter / 2.0, 2) * math.pi
  G_spi = dischargeCoefficient * math.sqrt(2 * densityBefore * (pressureBefore - pressureAfter)) if pressureBefore > pressureAfter else 0
  return G_spi * area



ambientTemperature = 293
filledPressure = 1 * constants.Pressure.bar
tankVolume = 35 * constants.Volume.liter

filledGasDensity = CP.PropsSI('D','T',ambientTemperature,'P',filledPressure,'N2O')
filledGasVolume = tankVolume
filledGasMass = filledGasDensity * filledGasVolume
gasSpecificInternalEnergy = CP.PropsSI('U','T',ambientTemperature,'Q',1,'N2O')
gasInternalEnergy = gasSpecificInternalEnergy * filledGasMass



initialState = [filledGasMass, ambientTemperature, filledGasVolume, filledPressure]
initialDerivatives = [0, 0, 0, 0]

print("d(U)/d(D)|T = {:.2f}".format(CP.PropsSI('d(U)/d(D)|T','P',101325,'T|gas',293,'N2O')))
print("Cv = {:.2f}".format(CP.PropsSI("CVMASS",'P',101325,'T|gas',293,'N2O')))
print("U = {:.2f}".format(CP.PropsSI("U",'P',101325,'T|gas',293,'N2O')))

print(CP.PropsSI("PMIN", "T", 293, "Q", 0, 'N2O'))

def computeResidual(t, y, ydot, residual):
  gasMass = y[0]
  gasTemperature = y[1]
  gasVolume = y[2]
  pressure = y[3]

  gasDensity = CP.PropsSI('D','T|gas',gasTemperature,'P',pressure,'N2O')
  hInlet = CP.PropsSI('H','T|gas',gasTemperature,'P',25*constants.Pressure.bar,'N2O')


  massFlowIn = 1e-4 * (np.sqrt(25*constants.Pressure.bar - pressure) if 25*constants.Pressure.bar > pressure else -np.sqrt(pressure - 25*constants.Pressure.bar))
  energyFlowIn = massFlowIn * hInlet

  dm_dt = massFlowIn
  dU_dt = energyFlowIn
  dV_dt = ydot[2]

  print(t, dm_dt, dU_dt)

  du_drho_T = CP.PropsSI('d(U)/d(D)|T','T|gas',gasTemperature,'P',25*constants.Pressure.bar,'N2O')
  cv = CP.PropsSI('CVMASS','T|gas',gasTemperature,'P',25*constants.Pressure.bar,'N2O')
  u = CP.PropsSI('U','T|gas',gasTemperature,'P',25*constants.Pressure.bar,'N2O')
  
  drho_dt = 1 / gasVolume * dm_dt - gasMass / (gasVolume * gasVolume) * dV_dt
  dT_dt = 1 / cv * (1 / gasMass * (dU_dt - u * dm_dt) - du_drho_T * drho_dt)

  residual[0] = ydot[0] - dm_dt
  residual[1] = ydot[1] - dT_dt
  residual[2] = y[2] - tankVolume
  residual[3] = gasMass / gasDensity - tankVolume

  return 0
#  hGas = CP.PropsSI('H','T',gasTemperature,'Q',1,'N2O')
#  hLiquid = CP.PropsSI('H','T',liquidTemperature,'Q',0,'N2O')
#
#  gasDensity = CP.PropsSI('D','T',gasTemperature,'Q',1,'N2O')
#  liquidDensity = CP.PropsSI('D','T',liquidTemperature,'Q',0,'N2O')
#  
#  gasVolume = gasMass / gasDensity
#  liquidVolume = liquidMass / liquidDensity
#
#  evaporationMassFlow = 0
#  condensationMassFlow = 0
#
#  gasMassFlow = evaporationMassFlow - condensationMassFlow
#  liquidMassFlow = condensationMassFlow - evaporationMassFlow
#
#  gasEnergyFlow = evaporationMassFlow * hLiquid - condensationMassFlow * hGas - pressure * ydot[4]
#  liquidEnergyFlow = -evaporationMassFlow * hLiquid + condensationMassFlow * hGas - pressure * ydot[5]
#  

#
#  massFlowLiquid = 0
#  massFlowGas = 0
#
#  hLiquid = CP.PropsSI('H','T',temperature,'Q',0,'N2O')
#  hGas = CP.PropsSI('H','T',temperature,'Q',1,'N2O')
#
#  energyFlowLiquid = hLiquid * massFlowLiquid
#  energyFlowGas = hGas * massFlowLiquid
#
#  liquidDensity = CP.PropsSI('D','T',temperature,'Q',0,'N2O')
#  gasDensity = CP.PropsSI('D','T',temperature,'Q',1,'N2O')
#
#  volume = liquidMass / liquidDensity + gasMass / gasDensity
#
#  residual[4] = volume - tankVolume
#  residual[5] = pressure -  CP.PropsSI('P','T',temperature,'Q',0,'N2O')
  

rtol = 1e-4
atol = 1e-4

algebraic_vars_idx = [1,2,3]

solver = dae('ida', computeResidual, 
    compute_initcond="yp0",
    exclude_algvar_from_error=False,
    algebraic_vars_idx=algebraic_vars_idx,
    atol=atol,
    rtol=rtol,
    max_steps=500000,
  )


times = np.linspace(0, 15, 100)
solution = solver.solve(times, initialState, initialDerivatives)

print(solution.message)

with open('/data/tankFilling.npy', 'wb') as f:
  t = np.array(solution.values.t)
  y = np.array(solution.values.y)
  np.save(f, t) 
  np.save(f, y) 