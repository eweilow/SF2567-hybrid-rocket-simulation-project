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

  def _reset(self):
    self.set(self._baseValue)

  @staticmethod
  def reset():
    for variable in variables:
      variable._reset()

class DerivedVariable:
  def __init__(self, compute):
    self.compute = compute
    self._adjustment = 1.0
    variables.append(self)

  def randomize(self, radius = 0.1):
    self._adjustment = random.uniform((1 - radius), (1 + radius))

  def get(self):
    return self.compute() * self._adjustment

  def _reset(self):
    self._adjustment = 1.0


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
tankInitialWallTemperature = DerivedVariable(lambda: tankFilledTemperature.get())

tankVolume = Variable(30 * constants.Volume.liter)
tankThickness = Variable(3.5 * constants.Lengths.mm)
tankInsideRadius = DerivedVariable(lambda: rocketBodyDiameter.get()/2 - tankThickness.get())
tankLength = DerivedVariable(lambda: tankVolume.get() / (math.pow(tankInsideRadius.get(), 2) * math.pi))


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
launchSeaLevelAltitude = Variable(0)


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
injectorStartupTransientTime = Variable(0.25)
injectorStartupDelay = Variable(0)
# aerospace-06-00075.pdf
maximumRegressionRateAt = Variable(constants.Pressure.bar * 30)

rocketDryMass = DerivedVariable(
  lambda: rocketOnBoardRecoverySystemMass.get() + rocketOnBoardElectronicsMass.get() + rocketPayloadMass.get() + rocketBodyMass.get() + rocketOxidizerTankMass.get() + rocketFuelCasingMass.get() + rocketEngineMass.get()
)


availableToRandomize = {
#  "tankOutletLiquidGasTransferRadius": { 
#    "title": "tankOutletLiquidGasTransferRadius", 
#    "var": tankOutletLiquidGasTransferRadius, 
#    "randomize": lambda var: var.randomize(0.25) 
#  },
#  "tankTopVentLiquidGasTransferRadius": { 
#    "title": "tankTopVentLiquidGasTransferRadius", 
#    "var": tankTopVentLiquidGasTransferRadius, 
#    "randomize": lambda var: var.randomize(0.25) 
#  },
  "initialAtmosphericPressure": { 
    "unit": "bar",
    "unitScale": 1 / constants.Pressure.bar,
    "title": "Initial ambient press.", 
    "var": initialAtmosphericPressure, 
    "randomize": lambda var: var.randomize(0.1) 
  },
  "initialAtmosphericTemperature": { 
    "unit": "K",
    "title": "Initial ambient temp.", 
    "var": initialAtmosphericTemperature, 
    "randomize": lambda var: var.randomize(0.1) 
  },
  "tankFillingGrade": { 
    "unit": "%",
    "unitScale": 100,
    "title": "Filling grade", 
    "var": tankFillingGrade, 
    "randomize": lambda var: var.randomizeInRange(0.5, 0.99) 
  },
  "tankFilledTemperature": { 
    "unit": "K",
    "title": "Initial oxidizer temp.", 
    "var": tankFilledTemperature, 
    "randomize": lambda var: var.randomizeInRange(268, 298) 
  },
  "tankVolume": { 
    "unit": "l",
    "unitScale": 1 / constants.Volume.liter,
    "title": "Tank volume", 
    "var": tankVolume, 
    "randomize": lambda var: var.randomize(0.1) 
  },
#  "tankInsideRadius": { 
#    "title": "tankInsideRadius", 
#    "var": tankInsideRadius, 
#    "randomize": lambda var: var.randomize(0.25) 
#  },
  "tankThickness": { 
    "title": "Tank wall thickness", 
    "var": tankThickness, 
    "randomize": lambda var: var.randomize(0.1) 
  },
  "tankWallDensity": { 
    "title": "Tank wall density", 
    "var": tankWallDensity, 
    "randomize": lambda var: var.randomize(0.1) 
  },
  "tankWallThermalConductivity": { 
    "title": "Tank wall heat cond.", 
    "var": tankWallThermalConductivity, 
    "randomize": lambda var: var.randomize(0.1) 
  },
  "tankWallSpecificHeatCapacity": { 
    "title": "Tank wall heat capac.", 
    "var": tankWallSpecificHeatCapacity, 
    "randomize": lambda var: var.randomize(0.1) 
  },
  "tankPassiveVentDischargeCoefficient": { 
    "title": "Passive vent Cd", 
    "var": tankPassiveVentDischargeCoefficient, 
    "randomize": lambda var: var.randomize(0.1) 
  },
  "tankPassiveVentDiameter": { 
    "title": "Passive vent diam.", 
    "var": tankPassiveVentDiameter, 
    "randomize": lambda var: var.randomize(0.1) 
  },
  "fuelDensitySolid": { 
    "title": "Fuel density (sol.)", 
    "var": fuelDensitySolid, 
    "randomize": lambda var: var.randomize(0.1) 
  },
  "fuelDensityLiquid": { 
    "title": "Fuel density (liq.)", 
    "var": fuelDensityLiquid, 
    "randomize": lambda var: var.randomize(0.1) 
  },
  "fuelEnthalpyOfFormation": { 
    "title": "Fuel enthalpy of formation", 
    "var": fuelEnthalpyOfFormation, 
    "randomize": lambda var: var.randomize(0.1) 
  },
  "fuelPortLength": { 
    "unit": "cm",
    "unitScale": 1/constants.Lengths.cm,
    "title": "Fuel port length", 
    "var": fuelPortLength, 
    "randomize": lambda var: var.randomize(0.1) 
  },
  "fuelPortInitialRadius": { 
    "unit": "mm",
    "unitScale": 1/constants.Lengths.mm,
    "title": "Fuel port radius", 
    "var": fuelPortInitialRadius, 
    "randomize": lambda var: var.randomize(0.1) 
  },
#  "fuelPortMaximumRadius": { 
#    "title": "fuelPortMaximumRadius", 
#    "var": fuelPortMaximumRadius, 
#    "randomize": lambda var: var.randomize(0.25) 
#  },
  "carbonBlackFraction": { 
    "title": "Fraction Carbon Black", 
    "var": carbonBlackFraction, 
    "randomize": lambda var: var.randomize(0.1) 
  },
  "preCombustionChamberVolume": { 
    "title": "Pre-CC volume", 
    "var": preCombustionChamberVolume, 
    "randomize": lambda var: var.randomize(0.1) 
  },
  "postCombustionChamberVolume": { 
    "title": "Post-CC volume", 
    "var": postCombustionChamberVolume, 
    "randomize": lambda var: var.randomize(0.1) 
  },
  "combustionEfficiency": { 
    "title": "combustionEfficiency", 
    "var": combustionEfficiency, 
    "randomize": lambda var: var.randomize(0.1) 
  },
  "fuelGrainAConstant": { 
    "title": "Regression a constant", 
    "var": fuelGrainAConstant, 
    "randomize": lambda var: var.randomize(0.1) 
  },
  "fuelGrainNConstant": { 
    "title": "Regression n constant", 
    "var": fuelGrainNConstant, 
    "randomize": lambda var: var.randomize(0.1) 
  },
  "injectorHoleCount": { 
    "title": "Injector # holes", 
    "var": injectorHoleCount, 
    "randomize": lambda var: var.randomize(0.1) 
  },
  "injectorHoleDischargeCoefficient": { 
    "title": "Injector Cd", 
    "var": injectorHoleDischargeCoefficient, 
    "randomize": lambda var: var.randomize(0.1) 
  },
  "injectorHoleDiameter": { 
    "unit": "mm",
    "unitScale": 1/constants.Lengths.mm,
    "title": "Injector hole diam.", 
    "var": injectorHoleDiameter, 
    "randomize": lambda var: var.randomize(0.1) 
  },
  "nozzleEfficiency": { 
    "title": "nozzleEfficiency", 
    "var": nozzleEfficiency, 
    "randomize": lambda var: var.randomize(0.1) 
  },
  "nozzleExhaustRadius": { 
    "unit": "mm",
    "unitScale": 1/constants.Lengths.mm,
    "title": "Exhaust radius", 
    "var": nozzleExhaustRadius, 
    "randomize": lambda var: var.randomize(0.1) 
  },
  "nozzleThroatRadius": { 
    "unit": "mm",
    "unitScale": 1/constants.Lengths.mm,
    "title": "Throat radius", 
    "var": nozzleThroatRadius, 
    "randomize": lambda var: var.randomize(0.1) 
  },
#  "nozzleErosionConstant": { 
#    "title": "nozzleErosionConstant", 
#    "var": nozzleErosionConstant, 
#    "randomize": lambda var: var.randomize(0.25) 
#  },
#  "nozzleErosionStart": { 
#    "title": "nozzleErosionStart", 
#    "var": nozzleErosionStart, 
#    "randomize": lambda var: var.randomize(0.25) 
#  },
#  "nozzleErosionStartRadius": { 
#    "title": "nozzleErosionStartRadius", 
#    "var": nozzleErosionStartRadius, 
#    "randomize": lambda var: var.randomize(0.25) 
#  },
#  "launchLatitudeDegrees": { 
#    "title": "launchLatitudeDegrees", 
#    "var": launchLatitudeDegrees, 
#    "randomize": lambda var: var.randomize(0.25) 
#  },
#  "launchLongitudeDegrees": { 
#    "title": "launchLongitudeDegrees", 
#    "var": launchLongitudeDegrees, 
#    "randomize": lambda var: var.randomize(0.25) 
#  },
#  "launchSeaLevelAltitude": { 
#    "title": "launchSeaLevelAltitude", 
#    "var": launchSeaLevelAltitude, 
#    "randomize": lambda var: var.randomize(0.25) 
#  },
  "launchTowerLength": { 
    "unit": "m",
    "title": "Tower length", 
    "var": launchTowerLength, 
    "randomize": lambda var: var.randomize(0.1) 
  },
  "launchTowerVerticalAngle": { 
    "unit": "°",
    "title": "Launch angle", 
    "var": launchTowerVerticalAngle, 
    "randomize": lambda var: var.randomize(0.1) 
  },
#  "launchTowerDirectionAngle": { 
#    "title": "launchTowerDirectionAngle", 
#    "var": launchTowerDirectionAngle, 
#    "randomize": lambda var: var.randomize(0.25) 
#  },
   "rocketBodyDiameter": { 
     "unit": "mm",
     "unitScale": 1/constants.Lengths.mm,
     "title": "Body diameter", 
     "var": rocketBodyDiameter, 
     "randomize": lambda var: var.randomize(0.1) 
   },

   "dryMass": {
     "unit": "kg",
     "title": "Dry mass",
     "var": rocketDryMass,
     "randomize": lambda var: var.randomize(0.1)
   },
#  "rocketOnBoardRecoverySystemMass": { 
#    "title": "rocketOnBoardRecoverySystemMass", 
#    "var": rocketOnBoardRecoverySystemMass, 
#    "randomize": lambda var: var.randomize(0.25) 
#  },
#  "rocketOnBoardElectronicsMass": { 
#    "title": "rocketOnBoardElectronicsMass", 
#    "var": rocketOnBoardElectronicsMass, 
#    "randomize": lambda var: var.randomize(0.25) 
#  },
#  "rocketPayloadMass": { 
#    "title": "rocketPayloadMass", 
#    "var": rocketPayloadMass, 
#    "randomize": lambda var: var.randomizeInRange(0, 20) 
#  },
#  "rocketBodyMass": { 
#    "unit": "kg",
#    "title": "Body mass", 
#    "var": rocketBodyMass, 
#    "randomize": lambda var: var.randomize(0.1) 
#  },
#  "rocketOxidizerTankMass": { 
#    "unit": "kg",
#    "title": "Tank mass", 
#    "var": rocketOxidizerTankMass, 
#    "randomize": lambda var: var.randomize(0.1) 
#  },
#  "rocketFuelCasingMass": { 
#    "title": "rocketFuelCasingMass", 
#    "var": rocketFuelCasingMass, 
#    "randomize": lambda var: var.randomize(0.25) 
#  },
#  "kastrullMass": { 
#    "unit": "kg",
#    "title": "Kastrullen mass", 
#    "var": kastrullMass, 
#    "randomize": lambda var: var.randomize(0.1) 
#  },
#  "fastenersMass": { 
#    "title": "fastenersMass", 
#    "var": fastenersMass, 
#    "randomize": lambda var: var.randomize(0.25) 
#  },
#  "injectorMass": { 
#    "title": "injectorMass", 
#    "var": injectorMass, 
#    "randomize": lambda var: var.randomize(0.25) 
#  },
#  "chamberWallsMass": { 
#    "title": "chamberWallsMass", 
#    "var": chamberWallsMass, 
#    "randomize": lambda var: var.randomize(0.25) 
#  },
#  "nozzleMass": { 
#    "title": "nozzleMass", 
#    "var": nozzleMass, 
#    "randomize": lambda var: var.randomize(0.25) 
#  },
#  "rocketEngineMass": { 
#    "unit": "kg",
#    "title": "Chamber mass", 
#    "var": rocketEngineMass, 
#    "randomize": lambda var: var.randomize(0.1) 
#  },
#  "dragLevelAtZero": { 
#    "title": "dragLevelAtZero", 
#    "var": dragLevelAtZero, 
#    "randomize": lambda var: var.randomize(0.25) 
#  },
#  "dragLevelAtPeak": { 
#    "title": "dragLevelAtPeak", 
#    "var": dragLevelAtPeak, 
#    "randomize": lambda var: var.randomize(0.25) 
#  },
#  "dragPeakMachNumber": { 
#    "title": "dragPeakMachNumber", 
#    "var": dragPeakMachNumber, 
#    "randomize": lambda var: var.randomize(0.25) 
#  },
#  "dragLevelMachAsymptote": { 
#    "title": "dragLevelMachAsymptote", 
#    "var": dragLevelMachAsymptote, 
#    "randomize": lambda var: var.randomize(0.25) 
#  },
#  "dragDropoffConstant": { 
#    "title": "dragDropoffConstant", 
#    "var": dragDropoffConstant, 
#    "randomize": lambda var: var.randomize(0.25) 
#  },
#  "dragPeakSmoothingRadius": { 
#    "title": "dragPeakSmoothingRadius", 
#    "var": dragPeakSmoothingRadius, 
#    "randomize": lambda var: var.randomize(0.25) 
#  },
#  "combustionEfficiencyStartupTransientTime": { 
#    "title": "combustionEfficiencyStartupTransientTime", 
#    "var": combustionEfficiencyStartupTransientTime, 
#    "randomize": lambda var: var.randomize(0.25) 
#  },
#  "combustionEfficiencyStartupDelay": { 
#    "title": "combustionEfficiencyStartupDelay", 
#    "var": combustionEfficiencyStartupDelay, 
#    "randomize": lambda var: var.randomize(0.25) 
#  },
#  "injectorStartupTransientTime": { 
#    "title": "injectorStartupTransientTime", 
#    "var": injectorStartupTransientTime, 
#    "randomize": lambda var: var.randomize(0.25) 
#  },
#  "injectorStartupDelay": { 
#    "title": "injectorStartupDelay", 
#    "var": injectorStartupDelay, 
#    "randomize": lambda var: var.randomize(0.25) 
#  },
  "maximumRegressionRateAt": { 
    "title": "Regression rate pressure peak", 
    "var": maximumRegressionRateAt, 
    "randomize": lambda var: var.randomize(0.1)
  }
}

def randomizeOne(mode = None):
  if mode is None:
    import random
    
    keys = list(availableToRandomize.keys())
    mode = random.choice(keys)
    print(mode)

  Variable.reset()
  var = availableToRandomize[mode]["var"]
  baseline = var.get()
  availableToRandomize[mode]["randomize"](var)
  current = var.get()

  return mode, baseline, current