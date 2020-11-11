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

  def randomize(self):
    self.randomizeInRange(self._baseValue * 0.9, self._baseValue * 1.1)

  def randomizeInRange(self, a, b):
    self.value = random.uniform(a, b)

  @staticmethod
  def reset():
    for variable in variables:
      variable.set(variable._baseValue)



tankVolume = Variable(35 * constants.Volume.liter)
tankFillingGrade = Variable(0.95)
tankFilledTemperature = Variable(293)

tankPassiveVentDischargeCoefficient = Variable(0.7)
tankPassiveVentDiameter = Variable(constants.Lengths.mm * 0.4)

fuelDensity = Variable(900)
fuelPortLength = Variable(33 * constants.Lengths.cm)
fuelPortInitialRadius = Variable(constants.Lengths.mm * 25)

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