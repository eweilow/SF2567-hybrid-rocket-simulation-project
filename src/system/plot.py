import matplotlib.pyplot as plt
import numpy as np

from models.tank import TankModel
from models.injector import InjectorModel
from models.combustion import CombustionModel
from models.nozzle import NozzleModel
from models.passiveVent import PassiveVentModel
from models.environment import EnvironmentModel
import utils.constants as constants


with open('./tmp/simulation.npy', 'rb') as f:
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
  }


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
plt.plot(t, models["tank"]["state"][TankModel.states_oxidizerMass], '--')
plt.plot(t, models["tank"]["derived"][TankModel.derived_liquidMass], '-')
plt.plot(t, models["tank"]["derived"][TankModel.derived_gasMass], '-')
plt.legend(("Total", "Liquid", "Gas"))
plt.xlabel("Time (s)")
plt.ylabel("Oxidizer mass (kg)")
plt.title("Oxidizer")

# nextSubplot()
# plt.plot(t, models["combustion"]["state"][CombustionModel.states_burntFuel], '-')
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

nextSubplot()
plt.plot(t, models["passiveVent"]["derived"][PassiveVentModel.derived_massFlow]*1e3, '-')
plt.xlabel("Time (s)")
plt.ylabel("Mass flow (g/s)")
plt.title("Passive vent")

nextSubplot()
plt.plot(t, models["tank"]["derived"][TankModel.derived_pressure] / constants.Pressure.bar, '-')
plt.plot(t, models["combustion"]["state"][CombustionModel.states_pressure] / constants.Pressure.bar, '-')
plt.plot(t, models["combustion"]["derived"][CombustionModel.derived_exhaustPressure] / constants.Pressure.bar, '-')
plt.plot(t, models["environment"]["derived"][EnvironmentModel.derived_ambientPressure] / constants.Pressure.bar, '-')
plt.legend(("Tank", "Chamber", "Exhaust", "Ambient"))
plt.xlabel("Time (s)")
plt.ylabel("Pressure (bar)")
plt.title("Pressures")

nextSubplot()
plt.plot(t, models["tank"]["derived"][TankModel.derived_outletDensity], '-')
plt.plot(t, models["tank"]["derived"][TankModel.derived_outletTopDensity], '-')
plt.legend(("Main outlet", "Passive vent"))
plt.xlabel("Time (s)")
plt.ylabel("Density (kg/m^3)")
plt.title("Tank outlets")

nextSubplot()
plt.plot(t, models["tank"]["derived"][TankModel.derived_temperature], '-')
plt.xlabel("Time (s)")
plt.ylabel("Temperature (K)")
plt.title("Tank temperature")

nextSubplot()
plt.plot(t, models["tank"]["derived"][TankModel.derived_liquidLevel] * 100, '-')
plt.xlabel("Time (s)")
plt.ylabel("Liquid level (%)")
plt.title("Tank liquid level")

nextSubplot()
plt.plot(t, models["tank"]["derived"][TankModel.derived_vaporQuality] * 100, '-')
plt.xlabel("Time (s)")
plt.ylabel("Vapor quality (%)")
plt.title("Tank vapor quality")

nextSubplot()
plt.plot(t, models["tank"]["derived"][TankModel.derived_liquidDensity], '-')
plt.plot(t, models["tank"]["derived"][TankModel.derived_gasDensity], '-')
plt.legend(('Liquid', 'Gas'))
plt.xlabel("Time (s)")
plt.ylabel("Density (kg/m^3)")
plt.title("Oxidizer density")


nextSubplot()
plt.plot(t, models["combustion"]["state"][CombustionModel.states_portRadius] / constants.Lengths.mm, '-')
plt.xlabel("Time (s)")
plt.ylabel("Radius (mm)")
plt.title("Port radius")


nextSubplot()
plt.plot(t, models["combustion"]["derived"][CombustionModel.derived_volume] / constants.Volume.liter, '-')
plt.xlabel("Time (s)")
plt.ylabel("Volume (liter)")
plt.title("Chamber volume")


nextSubplot()
plt.plot(t, models["combustion"]["derived"][CombustionModel.derived_temperature], '-')
plt.xlabel("Time (s)")
plt.ylabel("Temperature (K)")
plt.title("Adiabatic chamber flame temperature")

nextSubplot()
plt.plot(t, models["combustion"]["derived"][CombustionModel.derived_gamma], '-')
plt.xlabel("Time (s)")
plt.ylabel("Gamma")
plt.title("Combustion gamma")

nextSubplot()
plt.plot(t, models["combustion"]["derived"][CombustionModel.derived_CpT], '-')
plt.xlabel("Time (s)")
plt.ylabel("CpT")
plt.title("Combustion CpT")

nextSubplot()
plt.plot(t, models["combustion"]["derived"][CombustionModel.derived_rDot] / constants.Lengths.mm, '-')
plt.xlabel("Time (s)")
plt.ylabel("Rate (mm/s)")
plt.title("Regression rate")

nextSubplot()
plt.plot(t, models["combustion"]["derived"][CombustionModel.derived_ofRatio], '-')
plt.xlabel("Time (s)")
plt.ylabel("OF")
plt.title("Oxidizer-fuel ratio")

nextSubplot()
plt.plot(t, models["combustion"]["derived"][CombustionModel.derived_cStar], '-')
plt.xlabel("Time (s)")
plt.ylabel("C*")
plt.title("Characteristic velocity")

nextSubplot()
plt.plot(t, models["combustion"]["derived"][CombustionModel.derived_thrustCoefficient], '-')
plt.xlabel("Time (s)")
plt.ylabel("Coefficient")
plt.title("Thrust coefficient")

nextSubplot()
plt.plot(t, models["nozzle"]["derived"][NozzleModel.derived_thrust], '-')
plt.xlabel("Time (s)")
plt.ylabel("Thrust (N)")
plt.title("Thrust")

nextSubplot()
plt.plot(t, models["nozzle"]["derived"][NozzleModel.derived_specificImpulse], '-')
plt.xlabel("Time (s)")
plt.ylabel("Specific impulse (s)")
plt.title("Specific impulse")

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