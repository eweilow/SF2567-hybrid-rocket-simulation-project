import math
import CoolProp.CoolProp as CP

def computeHEMInjector(
  dischargeCoefficient,
  numberOfHoles,
  holeDiameter,
  densityBefore,
  pressureBefore,
  pressureAfter,
  temperatureBefore,
  phaseBefore,
  enthalpyBefore
):
  try:
    saturatedPressureBefore = CP.PropsSI('P','T',temperatureBefore,'Q',0,'N2O')
  except Exception as exc:
    print("Failed to evaluate saturatedPressureBefore", exc)
    return 0

  try:
    entropyBefore = CP.PropsSI('Smolar','T',temperatureBefore,'Q',phaseBefore,'N2O')
  except Exception as exc:
    print("Failed to evaluate entropyBefore", exc)
    return 0

  entropyBefore = max(entropyBefore, 100)
  try:
    enthalpyAfter = CP.PropsSI('H','P',pressureAfter,'Smolar',entropyBefore,'N2O')
  except Exception as exc:
    print("Failed to evaluate enthalpyAfter", exc)
    return 0

  try:
    densityAfter = CP.PropsSI('D','P',pressureAfter,'Smolar',entropyBefore,'N2O')
  except Exception as exc:
    print("Failed to evaluate densityAfter", exc)
    return 0

  # HEM injector model https://web.stanford.edu/~cantwell/Recent_publications/Zimmerman_et_al_AIAA_2013-4045.pdf
  kappa = math.sqrt((pressureBefore - pressureAfter) / (saturatedPressureBefore - pressureAfter))

  area = numberOfHoles * math.pow(holeDiameter / 2.0, 2) * math.pi
  G_spi = dischargeCoefficient * math.sqrt(2 * densityBefore * (pressureBefore - pressureAfter)) if pressureBefore > pressureAfter else 0
  G_hem = dischargeCoefficient * densityAfter * math.sqrt(2 * (enthalpyBefore - enthalpyAfter)) if enthalpyBefore > enthalpyAfter else 0
  G = (kappa * G_spi + G_hem) / (1 + kappa)

  return G * area
