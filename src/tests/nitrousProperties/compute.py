import CoolProp.CoolProp as CP
import numpy as np

class Propellant:
  def __init__(self, propellantName):
    self.propellantName = propellantName

  """
  Value in g/mol - not kg/mol!
  """
  def getMolarMass(self, temperature, vapourQuality):
    return CP.PropsSI('M','T',temperature,'Q',vapourQuality,self.propellantName) * 1000

  def getSpecificHeat(self, temperature, vapourQuality):
    return CP.PropsSI('CP0MOLAR','T',temperature,'Q',vapourQuality,self.propellantName)

  def getSaturationPressure(self, temperature, vapourQuality):
    return CP.PropsSI('P','T',temperature,'Q',vapourQuality,self.propellantName)

  def getSpecificEnthalpy(self, temperature, vapourQuality):
    return CP.PropsSI('H','T',temperature,'Q',vapourQuality,self.propellantName)

  def getSpecificInternalEnergy(self, temperature, vapourQuality):
    return CP.PropsSI('U','T',temperature,'Q',vapourQuality,self.propellantName)

  def getDensity(self, temperature, vapourQuality):
    return CP.PropsSI('D','T',temperature,'Q',vapourQuality,self.propellantName)

nitrousOxide = Propellant("N2O")

temperatures = np.linspace(230, 310, 1000)
liquidMolarMass = nitrousOxide.getMolarMass(temperatures, 0)
gaseousMolarMass = nitrousOxide.getMolarMass(temperatures, 1)
liquidSpecificHeat = nitrousOxide.getSpecificHeat(temperatures, 0)
gaseousSpecificHeat = nitrousOxide.getSpecificHeat(temperatures, 1)
liquidSaturationPressure = nitrousOxide.getSaturationPressure(temperatures, 0)
gaseousSaturationPressure = nitrousOxide.getSaturationPressure(temperatures, 1)
liquidSpecificEnthalpy = nitrousOxide.getSpecificEnthalpy(temperatures, 0)
gaseousSpecificEnthalpy = nitrousOxide.getSpecificEnthalpy(temperatures, 1)
liquidSpecificInternalEnergy = nitrousOxide.getSpecificInternalEnergy(temperatures, 0)
gaseousSpecificInternalEnergy = nitrousOxide.getSpecificInternalEnergy(temperatures, 1)
liquidDensity = nitrousOxide.getDensity(temperatures, 0)
gaseousDensity = nitrousOxide.getDensity(temperatures, 1)

with open('/data/nitrousProperties.npy', 'wb') as f:
  np.save(f, temperatures)
  np.save(f, liquidMolarMass)
  np.save(f, gaseousMolarMass)
  np.save(f, liquidSpecificHeat)
  np.save(f, gaseousSpecificHeat)
  np.save(f, liquidSaturationPressure)
  np.save(f, gaseousSaturationPressure)
  np.save(f, liquidSpecificEnthalpy)
  np.save(f, gaseousSpecificEnthalpy)
  np.save(f, liquidSpecificInternalEnergy)
  np.save(f, gaseousSpecificInternalEnergy)
  np.save(f, liquidDensity)
  np.save(f, gaseousDensity)