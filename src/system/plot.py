import matplotlib.pyplot as plt
import numpy as np

from models.tank import TankModel
from models.injector import InjectorModel
from models.combustion import CombustionModel
from models.nozzle import NozzleModel
from models.passiveVent import PassiveVentModel
from models.environment import EnvironmentModel
from models.flight import FlightModel
import utils.constants as constants


with open('./tmp/simulation.npy', 'rb') as f:
  [t0, t1, t2] = np.load(f)

  t = np.load(f)
  
  models = {
    "tank": {
      "state": np.load(f),
      "derived": np.load(f),
    },
    "injector": {
      "state": np.load(f),
      "derived": np.load(f),
    },
    "passiveVent": {
      "state": np.load(f),
      "derived": np.load(f),
    },
    "combustion": {
      "state": np.load(f),
      "derived": np.load(f),
    },
    "nozzle": {
      "state": np.load(f),
      "derived": np.load(f),
    },
    "environment": {
      "state": np.load(f),
      "derived": np.load(f),
    },
    "flight": {
      "state": np.load(f),
      "derived": np.load(f),
    },
  }

plt.subplot(2,2,1)
plt.grid()
downrange = np.linalg.norm([models["flight"]["state"][FlightModel.states_x], models["flight"]["state"][FlightModel.states_y]], axis=0)
altitude = models["flight"]["state"][FlightModel.states_z]
plt.plot(downrange, altitude)
plt.gca().set_aspect('equal', adjustable='box')
plt.xlabel("Downrange distance (m)")
plt.ylabel("Altitude (m)")

plt.subplot(2,2,2)
plt.grid()
plt.plot(t, models["flight"]["state"][FlightModel.states_x])
plt.plot(t, models["flight"]["state"][FlightModel.states_y])
plt.plot(t, models["flight"]["state"][FlightModel.states_z])
plt.legend(("x", "y", "z"))
plt.xlabel("Time (s)")
plt.ylabel("Position (m)")

plt.subplot(2,2,3)
plt.grid()
plt.plot(t, models["flight"]["state"][FlightModel.states_vx])
plt.plot(t, models["flight"]["state"][FlightModel.states_vy])
plt.plot(t, models["flight"]["state"][FlightModel.states_vz])
plt.legend(("x", "y", "z"))
plt.xlabel("Time (s)")
plt.ylabel("Velocity (m/s)")

plt.subplot(2,2,4)
plt.grid()
plt.plot(t, models["flight"]["derived"][FlightModel.derived_ax])
plt.plot(t, models["flight"]["derived"][FlightModel.derived_ay])
plt.plot(t, models["flight"]["derived"][FlightModel.derived_az])
plt.legend(("x", "y", "z"))
plt.xlabel("Time (s)")
plt.ylabel("Acceleration (m s-2)")


plt.tight_layout()
plt.show()
index = 1
width = 5
height = 4

# plt.plot(models["combustion"]["derived"][CombustionModel.derived_ofRatio], models["nozzle"]["derived"][NozzleModel.derived_specificImpulse])
# plt.xlabel("Oxidizer-Fuel Ratio")
# plt.ylabel("Specific Impulse [s]")
# plt.grid()
# plt.show()

def nextSubplot():
  global index
  print(height, width, index)
  plt.subplot(height, width, index)
  plt.grid()
  index = (index) % (width * height) + 1

nextSubplot()
plt.plot(t, models["tank"]["state"][TankModel.states_oxidizerMass] + models["combustion"]["state"][CombustionModel.states_fuelMass], '--')
plt.plot(t, models["combustion"]["state"][CombustionModel.states_fuelMass], '-')
plt.plot(t, models["tank"]["derived"][TankModel.derived_liquidMass], '-')
plt.plot(t, models["tank"]["derived"][TankModel.derived_gasMass], '-')
plt.legend(("Total", "Fuel", "Oxidizer (Liquid)", "Oxidizer (Gas)"))
plt.xlabel("Time (s)")
plt.ylabel("Propellant mass (kg)")
plt.title("Propellant")
plt.xlim(t0, t1)

# nextSubplot()
# plt.plot(t, models["combustion"]["state"][CombustionModel.states_fuelMass], '-')
# plt.xlabel("Time (s)")
# plt.ylabel("Burnt mass (kg)")
# plt.title("Fuel")


nextSubplot()

plt.plot(t, models["injector"]["derived"][InjectorModel.derived_massFlow], '-')
plt.plot(t, models["combustion"]["derived"][CombustionModel.derived_fuelFlow], '-')
plt.plot(t, models["nozzle"]["derived"][NozzleModel.derived_massFlow], '--')
plt.legend(("Oxidizer", "Fuel", "Nozzle"))
plt.xlabel("Time (s)")
plt.ylabel("Mass flow (kg/s)")
plt.title("Propellant flow")
plt.xlim(t0, t1)

nextSubplot()
plt.plot(t, models["passiveVent"]["derived"][PassiveVentModel.derived_massFlow]*1e3, '-')
plt.xlabel("Time (s)")
plt.ylabel("Mass flow (g/s)")
plt.title("Passive vent")
plt.xlim(t0, t1)

nextSubplot()
plt.plot(t, models["tank"]["derived"][TankModel.derived_pressure] / constants.Pressure.bar, '-')
plt.plot(t, models["combustion"]["state"][CombustionModel.states_pressure] / constants.Pressure.bar, '-')
plt.plot(t, models["combustion"]["derived"][CombustionModel.derived_exhaustPressure] / constants.Pressure.bar, '-')
plt.plot(t, models["environment"]["derived"][EnvironmentModel.derived_ambientPressure] / constants.Pressure.bar, '-')

plt.legend(("Tank", "Chamber", "Exhaust", "Ambient"))
plt.xlabel("Time (s)")
plt.ylabel("Pressure (bar)")
plt.title("Pressures")
plt.xlim(t0, t1)

nextSubplot()
plt.plot(t, models["tank"]["derived"][TankModel.derived_outletDensity], '-')
plt.plot(t, models["tank"]["derived"][TankModel.derived_outletTopDensity], '-')
plt.plot(t, models["combustion"]["derived"][CombustionModel.derived_density], '--')
plt.plot(t, models["combustion"]["derived"][CombustionModel.derived_oxidizerDensity], '--')
plt.legend(("Main outlet", "Passive vent", "Combustion mixture", "Oxidizer in chamber"))
plt.xlabel("Time (s)")
plt.ylabel("Density (kg/m^3)")
plt.title("Densities")
plt.xlim(t0, t1)

nextSubplot()
plt.plot(t, models["tank"]["derived"][TankModel.derived_temperature], '-')
plt.plot(t, models["tank"]["state"][TankModel.states_gasWallTankTemperature], '-')
plt.plot(t, models["tank"]["state"][TankModel.states_liquidWallTankTemperature], '-')
plt.xlabel("Time (s)")
plt.ylabel("Temperature (K)")
plt.title("Tank temperature")
plt.legend(("Oxidizer", "Tank wall (gas part)", "Tank wall (liquid part)"))
plt.xlim(t0, t1)

nextSubplot()
plt.plot(t, models["tank"]["derived"][TankModel.derived_liquidLevel] * 100, '-')
plt.xlabel("Time (s)")
plt.ylabel("Liquid level (%)")
plt.title("Tank liquid level")
plt.xlim(t0, t1)

nextSubplot()
plt.plot(t, models["tank"]["derived"][TankModel.derived_vaporQuality] * 100, '-')
plt.xlabel("Time (s)")
plt.ylabel("Vapor quality (%)")
plt.title("Tank vapor quality")
plt.xlim(t0, t1)

nextSubplot()
plt.plot(t, models["tank"]["derived"][TankModel.derived_liquidDensity], '-')
plt.plot(t, models["tank"]["derived"][TankModel.derived_gasDensity], '-')
plt.legend(('Liquid', 'Gas'))
plt.xlabel("Time (s)")
plt.ylabel("Density (kg/m^3)")
plt.title("Oxidizer density")
plt.xlim(t0, t1)


nextSubplot()
plt.plot(t, models["combustion"]["state"][CombustionModel.states_portRadius] / constants.Lengths.mm, '-')
plt.hlines(75 - 2, np.min(t), np.max(t), 'r', linestyles="dashed", label="Casing")
plt.hlines(75 - 2 - 10, np.min(t), np.max(t), 'g', linestyles="dashdot", label="Margin")
plt.legend()
plt.xlabel("Time (s)")
plt.ylabel("Radius (mm)")
plt.title("Port radius")
plt.xlim(t0, t1)


nextSubplot()
plt.plot(t, models["combustion"]["derived"][CombustionModel.derived_volume] / constants.Volume.liter, '-')
plt.xlabel("Time (s)")
plt.ylabel("Volume (liter)")
plt.title("Chamber volume")
plt.xlim(t0, t1)


nextSubplot()
plt.plot(t, models["combustion"]["derived"][CombustionModel.derived_temperature], '-')
plt.xlabel("Time (s)")
plt.ylabel("Temperature (K)")
plt.title("Adiabatic chamber flame temperature")
plt.xlim(t0, t1)

nextSubplot()
plt.plot(t, models["combustion"]["derived"][CombustionModel.derived_gamma], '-')
plt.xlabel("Time (s)")
plt.ylabel("Gamma")
plt.title("Combustion gamma")
plt.xlim(t0, t1)

nextSubplot()
plt.plot(t, models["combustion"]["derived"][CombustionModel.derived_CpT], '-')
plt.xlabel("Time (s)")
plt.ylabel("CpT")
plt.title("Combustion CpT")
plt.xlim(t0, t1)

nextSubplot()
plt.plot(t, models["combustion"]["derived"][CombustionModel.derived_rDot] / constants.Lengths.mm, '-')
plt.xlabel("Time (s)")
plt.ylabel("Rate (mm/s)")
plt.title("Regression rate")
plt.xlim(t0, t1)

nextSubplot()
plt.plot(t, models["combustion"]["derived"][CombustionModel.derived_ofRatio], '-')
plt.xlabel("Time (s)")
plt.ylabel("OF")
plt.title("Oxidizer-fuel ratio")
plt.xlim(t0, t1)

nextSubplot()
plt.plot(t, models["combustion"]["derived"][CombustionModel.derived_cStar], '-')
plt.xlabel("Time (s)")
plt.ylabel("C*")
plt.title("Characteristic velocity")
plt.xlim(t0, t1)

nextSubplot()
plt.plot(t, models["combustion"]["derived"][CombustionModel.derived_thrustCoefficient], '-')
plt.xlabel("Time (s)")
plt.ylabel("Coefficient")
plt.title("Thrust coefficient")
plt.xlim(t0, t1)

nextSubplot()
plt.plot(t, models["nozzle"]["derived"][NozzleModel.derived_thrust], '-')
plt.xlabel("Time (s)")
plt.ylabel("Thrust (N)")
plt.title("Thrust")
plt.xlim(t0, t1)

nextSubplot()
plt.plot(t, models["nozzle"]["derived"][NozzleModel.derived_specificImpulse], '-')
plt.xlabel("Time (s)")
plt.ylabel("Specific impulse (s)")
plt.title("Specific impulse")
plt.xlim(t0, t1)

plt.tight_layout()
plt.subplots_adjust(
  left=0.05,
  bottom=0.05,
  top=1-0.05,
  right=1-0.05,
  wspace=0.25,
  hspace=0.35
)
plt.show()