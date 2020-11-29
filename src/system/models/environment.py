from models.base import Model
import math
import CoolProp.CoolProp as CP
import numpy as np
import scipy.interpolate

import assumptions

def init():
  global seaLevelTemperature, seaLevelPressure
  from ambiance import Atmosphere

  seaLevelAtmosphere = Atmosphere(0)
  seaLevelPressure = seaLevelAtmosphere.pressure[0]
  seaLevelTemperature = seaLevelAtmosphere.temperature[0]

class EnvironmentModel(Model):
  def derivativesDependsOn(self, models):
    return []

  def derivedVariablesDependsOn(self, models):
    return []
  
  derived_ambientPressure = 0
  derived_ambientDensity = 1
  derived_ambientSpeedOfSound = 2
  derived_gravityVerticalComponent = 3
  derived_gravityNorthwardComponent = 4
  derived_ambientViscosity = 5
  derived_ambientThermalConductivity = 6
  derived_ambientTemperature = 7

  derived_dynamicPressure = 8
  derived_stagnationPressure = 9
  derived_stagnationTemperature = 10

  def initializeState(self):
    init()

    return [0]

  def computeDerivatives(self, t, state, derived, models):
    return [0]

  def computeDerivedVariables(self, t, state, models):
    from models.flight import FlightModel

    from rellipsoid import earth
    from ambiance import Atmosphere

    surfaceAltitude = assumptions.launchSeaLevelAltitude.get() + models["flight"]["state"][FlightModel.states_z]
    if surfaceAltitude < 0:
      surfaceAltitude = 0
    if surfaceAltitude > 80000:
      surfaceAltitude = 80000

    atmosphere = Atmosphere(surfaceAltitude)
    pressure = atmosphere.pressure[0]
    density = atmosphere.density[0]
    speedOfSound = atmosphere.speed_of_sound[0]
    temperature = atmosphere.temperature[0]

    pressure = pressure * (assumptions.initialAtmosphericPressure.get() / seaLevelPressure)
    density = density * (assumptions.initialAtmosphericPressure.get() / seaLevelPressure)
    speedOfSound = speedOfSound * (assumptions.initialAtmosphericPressure.get() / seaLevelPressure)
    temperature = temperature * (assumptions.initialAtmosphericTemperature.get() / seaLevelTemperature)
    viscosity = atmosphere.dynamic_viscosity[0]
    thermalConductivity = atmosphere.thermal_conductivity[0]

    velocityThroughAir = np.linalg.norm([models["flight"]["state"][FlightModel.states_vx], models["flight"]["state"][FlightModel.states_vy], models["flight"]["state"][FlightModel.states_vz]])
    dynamicPressure = 0.5 * density * velocityThroughAir * velocityThroughAir

    stagnationPressure = pressure - dynamicPressure

    airCp = CP.PropsSI('CPMASS','T',temperature,'P',pressure,'air')
    airCv = CP.PropsSI('CVMASS','T',temperature,'P',pressure,'air')
    airGamma = airCp / airCv

    mach = velocityThroughAir / speedOfSound

    # This is the adiabatic stagnation temperature on the surface of the rocket
    stagnationTemperature = temperature * (1 + (airGamma - 1) / 2 * mach * mach)

    launchLatitudeRadians = assumptions.launchLatitudeDegrees.get() / 180 * math.pi

    # todo for future: take into account the present latitude, not just the starting one
    verticalGravity, northwardGravity = earth.get_analytic_gravity(launchLatitudeRadians, surfaceAltitude)

    return [pressure, density, speedOfSound, verticalGravity, northwardGravity, viscosity, thermalConductivity, temperature, dynamicPressure, stagnationPressure, stagnationTemperature]