from models.base import Model
import math
import numpy as np
import scipy.interpolate

import assumptions
import options


hasInited = False
def init():
  global hasInited
  global pressureInterpolator, densityInterpolator, speedOfSoundInterpolator, seaLevelPressure
  from ambiance import Atmosphere
  if options.enableAtmosphereInterpolation and not hasInited:
    altitudes = np.linspace(0, 20e3, options.atmosphereInterpolantPointCount)
    atmosphere = Atmosphere(altitudes)
    interpolant = options.atmosphereInterpolant
    pressureInterpolator = scipy.interpolate.interp1d(altitudes, atmosphere.pressure, kind=interpolant, bounds_error=True, copy=True)
    densityInterpolator = scipy.interpolate.interp1d(altitudes, atmosphere.density, kind=interpolant, bounds_error=True, copy=True)
    speedOfSoundInterpolator = scipy.interpolate.interp1d(altitudes, atmosphere.speed_of_sound, kind=interpolant, bounds_error=True, copy=True)
    hasInited = True
  
  seaLevelAtmosphere = Atmosphere(0)
  seaLevelPressure = seaLevelAtmosphere.pressure[0]

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

  def initializeState(self):
    init()

    return [0]

  def computeDerivatives(self, t, state, derived, models):
    return [0]

  def computeDerivedVariables(self, t, state, models):
    from models.flight import FlightModel

    from rellipsoid import earth
    from ambiance import Atmosphere

    surfaceAltitude = models["flight"]["state"][FlightModel.states_z]
    if surfaceAltitude < 0:
      surfaceAltitude = 0

    if options.enableAtmosphereInterpolation:
      pressure = float(pressureInterpolator(surfaceAltitude))
      density = float(densityInterpolator(surfaceAltitude))
      speedOfSound = float(speedOfSoundInterpolator(surfaceAltitude))
    else:
      atmosphere = Atmosphere(surfaceAltitude)
      pressure = atmosphere.pressure[0]
      density = atmosphere.density[0]
      speedOfSound = atmosphere.speed_of_sound[0]

    pressure = pressure * (assumptions.initialAtmosphericPressure.get() / seaLevelPressure)
    density = density * (assumptions.initialAtmosphericPressure.get() / seaLevelPressure)
    speedOfSound = speedOfSound * (assumptions.initialAtmosphericPressure.get() / seaLevelPressure)

    launchLatitudeRadians = assumptions.launchLatitudeDegrees.get() / 180 * math.pi

    # todo for future: take into account the present latitude, not just the starting one
    verticalGravity, northwardGravity = earth.get_analytic_gravity(launchLatitudeRadians, assumptions.launchSeaLevelAltitude.get())

    return [pressure, density, speedOfSound, verticalGravity, northwardGravity]