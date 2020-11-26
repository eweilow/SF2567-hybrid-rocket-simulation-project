import random
import math
from utils import constants

variables = []
class Variable:
  def __init__(self, value):
    self._baseValue = value
    self.value = value
    variables.append(self)

  def get(self):
    return self.value

  def set(self, value):
    self.value = value

  def randomize(self, radius = 0.1):
    self.randomizeInRange(self._baseValue * (1 - radius), self._baseValue * (1 + radius))

  def randomizeInRange(self, a, b):
    self.value = random.uniform(a, b)

  @staticmethod
  def reset():
    for variable in variables:
      variable.set(variable._baseValue)

class DerivedVariable:
  def __init__(self, compute):
    self.compute = compute

  def get(self):
    return self.compute()



# The "radius" in terms of liquid level [0..1] where the tank outlet thermodynamic phase "smoothly" transitions from liquid to gas.
# Unfortunately this is an empirical constant which have to be determined after a test.
tankOutletLiquidGasTransferRadius = Variable(0.01)

# The "radius" in terms of liquid level [0..1] where the tank top vent thermodynamic phase "smoothly" transitions from gas to liquid.
# Unfortunately this is an empirical constant.
tankTopVentLiquidGasTransferRadius = Variable(0.01)


initialAtmosphericPressure = Variable(101300)
initialAtmosphericTemperature = Variable(293)

tankFillingGrade = Variable(0.95)
tankFilledTemperature = Variable(293)
tankInitialWallTemperature = Variable(293)

tankLength = Variable(1.87 * constants.Lengths.m)
tankInsideRadius = Variable(75 * constants.Lengths.mm)
tankThickness = Variable(3.5 * constants.Lengths.mm)

tankVolume = DerivedVariable(lambda: tankLength.get() * math.pow(tankInsideRadius.get(), 2) * math.pi)

tankWallDensity = Variable(2700) # http://www.aalco.co.uk/datasheets/Aluminium-Alloy-6082-T6T651-Plate_148.ashx
tankWallThermalConductivity = Variable(180) # http://www.aalco.co.uk/datasheets/Aluminium-Alloy-6082-T6T651-Plate_148.ashx
tankWallSpecificHeatCapacity = Variable(890) # https://matmatch.com/materials/mitf374-bs-en-573-3-grade-6082-t6

tankPassiveVentDischargeCoefficient = Variable(0.7)
tankPassiveVentDiameter = Variable(constants.Lengths.mm * 0.4)

fuelDensitySolid = Variable(900) # page 23 - https://drive.google.com/file/d/1o-I2h1bSd9EcnG2-ctqdSLR2UjsFNFpF/view?usp=sharing
fuelDensityLiquid = Variable(720) # page 23 - https://drive.google.com/file/d/1o-I2h1bSd9EcnG2-ctqdSLR2UjsFNFpF/view?usp=sharing
fuelEnthalpyOfFormation = Variable(-1438.2) # page 23 - https://drive.google.com/file/d/1o-I2h1bSd9EcnG2-ctqdSLR2UjsFNFpF/view?usp=sharing

fuelPortLength = Variable(33 * constants.Lengths.cm)
fuelPortInitialRadius = Variable(constants.Lengths.mm * 25)
fuelPortMaximumRadius = Variable(constants.Lengths.mm * 70)
carbonBlackFraction = Variable(0.02)

preCombustionChamberVolume = Variable(1 * constants.Volume.liter)
postCombustionChamberVolume = Variable(1 * constants.Volume.liter)

combustionEfficiency = Variable(0.9)

fuelGrainAConstant = Variable(0.155e-3) # 0.132e-3 to 0.155e-3. page 37 - https://drive.google.com/drive/folders/1zWr8Qbn6sgkfRpe2k6DzyGvg2Ws2mxZb
fuelGrainNConstant = Variable(0.5) # 0.555 to 0.5 page 37 - https://drive.google.com/drive/folders/1zWr8Qbn6sgkfRpe2k6DzyGvg2Ws2mxZb

injectorHoleCount = Variable(30)
injectorHoleDischargeCoefficient = Variable(0.83)
injectorHoleDiameter = Variable(constants.Lengths.mm * 1.5)

nozzleEfficiency = Variable(0.9)
nozzleExhaustRadius = Variable(constants.Lengths.mm * 44.777)
nozzleThroatRadius = Variable(constants.Lengths.mm * 19.466)
nozzleErosionConstant = Variable(constants.Lengths.mm * 0) # Set to 0 for now
nozzleErosionStart = Variable(10)
nozzleErosionStartRadius = Variable(3)

launchLatitudeDegrees = Variable(61)
launchLongitudeDegrees = Variable(14)
launchSeaLevelAltitude = Variable(10000)


launchTowerLength = Variable(10)
launchTowerVerticalAngle = Variable(80) # 90 degrees = straight upward
launchTowerDirectionAngle = Variable(0) # 0 degrees = due north, 90 degrees = due west


rocketBodyDiameter = Variable(16 * constants.Lengths.cm)

rocketOnBoardRecoverySystemMass = Variable(6)
rocketOnBoardElectronicsMass = Variable(2.3)
rocketPayloadMass = Variable(2)
rocketBodyMass = Variable(10)

rocketOxidizerTankMass = Variable(6.7)
rocketFuelCasingMass = Variable(0.3)

kastrullMass = Variable(3.3)
fastenersMass = Variable(0.1)
injectorMass = Variable(1.1)
chamberWallsMass = Variable(3.6)
nozzleMass = Variable(3.8)
rocketEngineMass = DerivedVariable(lambda: kastrullMass.get() + fastenersMass.get() + injectorMass.get() + chamberWallsMass.get() + nozzleMass.get())


dragLevelAtZero = Variable(0.5)
dragLevelAtPeak = Variable(0.7)
dragPeakMachNumber = Variable(1.1)
dragLevelMachAsymptote = Variable(0.2)
dragDropoffConstant = Variable(0.6)
dragPeakSmoothingRadius = Variable(0.5)

combustionEfficiencyStartupTransientTime = Variable(0.25)
combustionEfficiencyStartupDelay = Variable(0)
injectorStartupTransientTime = Variable(0.1)
injectorStartupDelay = Variable(0)
# aerospace-06-00075.pdf
maximumRegressionRateAt = Variable(constants.Pressure.bar * 30)
