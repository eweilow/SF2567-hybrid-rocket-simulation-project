from models.base import Model
import math

import assumptions

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
    return [0]

  def computeDerivatives(self, t, state, derived, models):
    return [0]

  def computeDerivedVariables(self, t, state, models):
    from models.flight import FlightModel

    from rellipsoid import earth
    from ambiance import Atmosphere

    surfaceAltitude = models["flight"]["state"][FlightModel.states_z]

    atmosphere = Atmosphere(surfaceAltitude)

    launchLatitudeRadians = assumptions.launchLatitudeDegrees.get() / 180 * math.pi

    # todo for future: take into account the present latitude, not just the starting one
    verticalGravity, northwardGravity = earth.get_analytic_gravity(launchLatitudeRadians, assumptions.launchSeaLevelAltitude.get())

    return [atmosphere.pressure[0], atmosphere.density[0], atmosphere.speed_of_sound[0], verticalGravity, northwardGravity]