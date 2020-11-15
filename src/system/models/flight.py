from models.base import Model
import math
import numpy as np
from utils import constants
from equations.drag import getDragCoefficient

import assumptions

class FlightModel(Model):
  def derivativesDependsOn(self, models):
    return []

  def derivedVariablesDependsOn(self, models):
    return [models["nozzle"], models["environment"]]

  def initializeState(self):
    # x: north, y: west, z: vertical
    # x, y, z, vx, vy, vz
    return [0, 0, 0, 0, 0, 0]

  states_x = 0
  states_northwardPosition = states_x
  states_y = 1
  states_westwardPosition = states_y
  states_z = 2
  states_verticalPosition = states_z

  states_vx = 3
  states_vy = 4
  states_vz = 5

  derived_ax = 0
  derived_ay = 1
  derived_az = 2
  derived_onTower = 3
  derived_mach = 4
  derived_propellantMass = 5
  derived_rocketMass = 6
  derived_totalMass = 7

  def computeDerivatives(self, t, state, derived, models):
    return [state[self.states_vx], state[self.states_vy], state[self.states_vz], derived[self.derived_ax], derived[self.derived_ay], derived[self.derived_az]]

  def computeDerivedVariables(self, t, state, models):
    from models.nozzle import NozzleModel
    from models.environment import EnvironmentModel
    from models.tank import TankModel
    from models.combustion import CombustionModel

    thrust = models["nozzle"]["derived"][NozzleModel.derived_thrust]

    northwardGravity = models["environment"]["derived"][EnvironmentModel.derived_gravityNorthwardComponent]
    verticalGravity = models["environment"]["derived"][EnvironmentModel.derived_gravityVerticalComponent]
    airDensity = models["environment"]["derived"][EnvironmentModel.derived_ambientDensity]
    speedOfSound = models["environment"]["derived"][EnvironmentModel.derived_ambientSpeedOfSound]

    towerX = np.cos(assumptions.launchTowerDirectionAngle.get() / 180 * math.pi) * np.cos(assumptions.launchTowerVerticalAngle.get() / 180 * math.pi)
    towerY = np.sin(assumptions.launchTowerDirectionAngle.get() / 180 * math.pi) * np.cos(assumptions.launchTowerVerticalAngle.get() / 180 * math.pi)
    towerZ = np.sin(assumptions.launchTowerVerticalAngle.get() / 180 * math.pi)

    position = np.array([state[self.states_x], state[self.states_y], state[self.states_z]])
    velocity = np.array([state[self.states_vx], state[self.states_vy], state[self.states_vz]])
    launchTowerVector = np.array([towerX, towerY, towerZ])

    projection = np.dot(launchTowerVector, position)
    projectionVector = projection * launchTowerVector
    perpendicularVector = position - projectionVector
    perpendicular = np.linalg.norm(perpendicularVector)

    onTower = False
    # some tolerance
    if perpendicular < 5 * constants.Lengths.cm and projection < assumptions.launchTowerLength.get():
      onTower = True
    
    direction = launchTowerVector if onTower else velocity
    direction = direction / np.linalg.norm(direction)

    velocityNorm = np.linalg.norm(velocity)
    mach = velocityNorm / speedOfSound

    dragCoefficient = getDragCoefficient(mach)

    radius = assumptions.rocketBodyDiameter.get() / 2.0
    drag = dragCoefficient * 0.5 * math.pow(radius, 2)*math.pi * airDensity * velocityNorm*velocityNorm


    oxidizerMass = models["tank"]["state"][TankModel.states_oxidizerMass]
    fuelMass = models["combustion"]["state"][CombustionModel.states_fuelMass]
    chamberPropellantMass = models["combustion"]["state"][CombustionModel.states_propellantMass]
    propellantMass = oxidizerMass + fuelMass + chamberPropellantMass
    rocketMass = assumptions.rocketOnBoardRecoverySystemMass.get() + assumptions.rocketOnBoardElectronicsMass.get() + assumptions.rocketPayloadMass.get() + assumptions.rocketBodyMass.get() + assumptions.rocketOxidizerTankMass.get() + assumptions.rocketFuelCasingMass.get() + assumptions.rocketEngineMass.get()

    totalMass = propellantMass + rocketMass

    ## UGH masses
    dragAcceleration = -drag * direction / totalMass
    thrustAcceleration = thrust * direction / totalMass
    gravityAcceleration = np.array([northwardGravity, 0, verticalGravity])
    
    initialAcceleration = dragAcceleration + thrustAcceleration + gravityAcceleration

    if onTower:
      accelerationProjection = np.dot(launchTowerVector, initialAcceleration)
      if projection <= 0 and accelerationProjection < 0:
        accelerationProjection = 0
      acceleration = accelerationProjection * launchTowerVector
      frictionComponent = np.linalg.norm(acceleration - initialAcceleration)
    else:
      acceleration = initialAcceleration
    

      

    return [acceleration[0], acceleration[1], acceleration[2], 1 if onTower else 0, mach, propellantMass, rocketMass, totalMass]