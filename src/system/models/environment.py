from models.base import Model
import math

import assumptions

class EnvironmentModel(Model):
  def derivativesDependsOn(self, models):
    return []

  def derivedVariablesDependsOn(self, models):
    return []
  
  derived_ambientPressure = 0

  def initializeState(self):
    return [0]

  def computeDerivatives(self, t, state, derived, models):
    return [0]

  def computeDerivedVariables(self, t, state, models):
    from rellipsoid import earth
    from ambiance import Atmosphere


    surfaceAltitude = 0

    atmosphere = Atmosphere(surfaceAltitude)

    launchLatitudeRadians = assumptions.launchLatitudeDegrees.get() / 180 * math.pi
    verticalGravity, northwardGravity = earth.get_analytic_gravity(launchLatitudeRadians, assumptions.launchSeaLevelAltitude.get())

    return [atmosphere.pressure[0]]