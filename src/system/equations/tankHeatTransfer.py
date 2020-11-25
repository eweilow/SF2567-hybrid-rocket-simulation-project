import math

# https://web.stanford.edu/~cantwell/Recent_publications/Zimmerman_et_al_AIAA_2013-4045.pdf
def tankWallHeatTransfer(
  fluidThermalConductivity,
  fluidCp,
  fluidMu,
  fluidDensity,
  acceleration,
  fluidTemperature,
  wallTemperature,
  tankWallHeight,
  tankWallRadius,
  n = 2/5,
  c = 0.021
):
  beta = 1 / fluidTemperature # Assumption due to ideal gas

  temperatureDifference = wallTemperature - fluidTemperature
  tankWallArea = tankWallHeight * 2 * tankWallRadius * math.pi
  
  Ra = fluidCp * fluidDensity * fluidDensity * acceleration * beta * abs(temperatureDifference) * pow(tankWallHeight, 3) / (fluidMu * fluidThermalConductivity)
  Nu = c * pow(Ra, n)
  heatTransferCoefficient = Nu * fluidThermalConductivity / tankWallHeight

  return heatTransferCoefficient * tankWallArea * temperatureDifference

def tankWallHeatConduction(
  thermalConductivity,
  tankWallInnerRadius,
  tankWallOuterRadius,
  topTemperature,
  bottomTemperature,
  totalLength,
  liquidLevel
):
  tankWallCrossSectionalArea = math.pi * (pow(tankWallOuterRadius, 2) - pow(tankWallInnerRadius, 2))
  temperatureDifference = bottomTemperature - topTemperature

  topPhaseLength = (totalLength - liquidLevel)
  topPhaseCenter = liquidLevel + topPhaseLength / 2

  bottomPhaseLength = liquidLevel
  bottomPhaseCenter = bottomPhaseLength / 2

  phaseDistance = topPhaseCenter - bottomPhaseCenter

  return thermalConductivity * temperatureDifference * tankWallCrossSectionalArea / phaseDistance


def estimate_liquidLevelChangeDerivative(
  liquidMassFlow,
  liquidDensity,
  tankWallInnerRadius
):
  # https://github.com/aesirkth/mjollnir-propulsion-simulations/blob/master/combustion/models/tankHeatTransferModel.m#L140
  volumetricChange = liquidMassFlow / liquidDensity # m^3 / s
  tankInnerCrossSectionalArea = pow(tankWallInnerRadius, 2) * math.pi
  heightChange = volumetricChange / tankInnerCrossSectionalArea
  return heightChange
  
def massPerLength(
  tankWallInnerRadius,
  tankWallOuterRadius,
  tankWallDensity,
):
  tankWallCrossSectionalArea = math.pi * (pow(tankWallOuterRadius, 2) - pow(tankWallInnerRadius, 2))
  return tankWallCrossSectionalArea * tankWallDensity
  

def tankWallControlVolumeBoundaryTransfer(
  tankWallInnerRadius,
  tankWallOuterRadius,
  tankWallDensity,
  liquidLevelChangeDerivative,
  wallHeatCapacity,
  topTemperature,
  bottomTemperature
):
  tankWallCrossSectionalAreaChange = liquidLevelChangeDerivative * math.pi * (pow(tankWallOuterRadius, 2) - pow(tankWallInnerRadius, 2))
  mDotWall = tankWallCrossSectionalAreaChange * tankWallDensity
  return mDotWall * wallHeatCapacity * (topTemperature - bottomTemperature)