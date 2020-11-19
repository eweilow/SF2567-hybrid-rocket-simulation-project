import random
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


initialAtmosphericPressure = Variable(101300)

tankVolume = Variable(35 * constants.Volume.liter)
tankFillingGrade = Variable(0.9)
tankFilledTemperature = Variable(293)

tankPassiveVentDischargeCoefficient = Variable(0.7)
tankPassiveVentDiameter = Variable(constants.Lengths.mm * 0.4)

fuelDensity = Variable(900)
fuelPortLength = Variable(33 * constants.Lengths.cm)
fuelPortInitialRadius = Variable(constants.Lengths.mm * 25)
fuelPortMaximumRadius = Variable(constants.Lengths.mm * 70)

preCombustionChamberVolume = Variable(1 * constants.Volume.liter)
postCombustionChamberVolume = Variable(1 * constants.Volume.liter)

combustionEfficiency = Variable(0.9)

#http://www.scielo.org.za/pdf/rd/v33/05.pdf
fuelGrainAConstant = Variable(0.155e-3)
fuelGrainNConstant = Variable(0.5)


injectorHoleCount = Variable(38)
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
launchTowerVerticalAngle = Variable(88) # 90 degrees = straight upward
launchTowerDirectionAngle = Variable(0) # 0 degrees = due north, 90 degrees = due west


rocketBodyDiameter = Variable(16 * constants.Lengths.cm)

rocketOnBoardRecoverySystemMass = Variable(10)
rocketOnBoardElectronicsMass = Variable(2.3)
rocketPayloadMass = Variable(10)
rocketBodyMass = Variable(7)

rocketOxidizerTankMass = Variable(6.7)
rocketFuelCasingMass = Variable(0.3)
rocketEngineMass = Variable(10.9)


dragBaseLevel = Variable(0.2)
dragPeakAroundMach = Variable(1.1)
dragPeak = Variable(0.7)
dragPeakSmoothingRadius = Variable(0.5)





combustionEfficiencyStartupTransientTime = Variable(0.25)
injectorStartupTransientTime = Variable(0.1)






# aerospace-06-00075.pdf
maximumRegressionRateAt = Variable(constants.Pressure.bar * 30)
