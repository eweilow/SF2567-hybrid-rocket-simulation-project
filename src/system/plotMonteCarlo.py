import matplotlib.pyplot as plt
import numpy as np
import scipy.integrate

from models.tank import TankModel
from models.injector import InjectorModel
from models.combustion import CombustionModel
from models.nozzle import NozzleModel
from models.passiveVent import PassiveVentModel
from models.environment import EnvironmentModel
from models.flight import FlightModel
import utils.constants as constants

with open('./tmp/montecarlo.npy', 'rb') as f:
  [N] = np.load(f)

  initialTemperatures = []
  fillingGrades = []
  peakThrusts = []
  meanThrusts = []
  peakPressures = []
  meanPressures = []
  peakCCPressures = []
  meanCCPressures = []
  impulses = []
  velocities = []

  for i in range(N):
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

    initialTemperature = models["tank"]["derived"][TankModel.derived_temperature][0] - 273.15
    fillingGrade = models["tank"]["derived"][TankModel.derived_liquidLevel][0]
    peakThrust = np.max(models["nozzle"]["derived"][NozzleModel.derived_thrust]) / 1000
    meanThrust = np.mean(models["nozzle"]["derived"][NozzleModel.derived_thrust]) / 1000
    peakPressure = np.max(models["tank"]["derived"][TankModel.derived_pressure]) / constants.Pressure.bar
    meanPressure = np.mean(models["tank"]["derived"][TankModel.derived_pressure]) / constants.Pressure.bar
    peakCCPressure = np.max(models["combustion"]["state"][CombustionModel.states_pressure]) / constants.Pressure.bar
    meanCCPressure = np.mean(models["combustion"]["state"][CombustionModel.states_pressure]) / constants.Pressure.bar
    velocity = np.linalg.norm([models["flight"]["state"][FlightModel.states_vx][-1], models["flight"]["state"][FlightModel.states_vy][-1], models["flight"]["state"][FlightModel.states_vz][-1]])

    totalImpulse = scipy.integrate.simps(models["nozzle"]["derived"][NozzleModel.derived_thrust], t) / 1000

    initialTemperatures.append(initialTemperature)
    fillingGrades.append(fillingGrade * 100)
    peakThrusts.append(peakThrust)
    meanThrusts.append(meanThrust)
    peakPressures.append(peakPressure)
    meanPressures.append(meanPressure)
    peakCCPressures.append(peakCCPressure)
    meanCCPressures.append(meanCCPressure)
    impulses.append(totalImpulse)
    velocities.append(velocity)

    plt.subplot(3,4,1)
    plt.plot(t, models["nozzle"]["derived"][NozzleModel.derived_thrust] / 1000, '-', linewidth=1)

    plt.subplot(3,4,2)
    plt.plot(t, models["flight"]["state"][FlightModel.states_z], '-', linewidth=1)

    plt.subplot(3,4,3)
    plt.plot(t, models["tank"]["derived"][TankModel.derived_pressure] / constants.Pressure.bar, '-', linewidth=1)
  
  plt.subplot(3,4,1)
  plt.title("Thrust for all sampled series")
  plt.xlabel("Time (s)")
  plt.ylabel("Thrust (kN)")
  plt.margins(y=0.1)

  plt.subplot(3,4,2)
  plt.title("Altitude for all sampled series")
  plt.xlabel("Time (s)")
  plt.ylabel("Altitude (m)")
  plt.margins(y=0.1)

  plt.subplot(3,4,3)
  plt.title("Tank pressure for all sampled series")
  plt.xlabel("Time (s)")
  plt.ylabel("Tank pressure (bar)")
  plt.margins(y=0.1)

  times = np.load(f)

  initialTemperatures = np.array(initialTemperatures)
  fillingGrades = np.array(fillingGrades)
  peakThrusts = np.array(peakThrusts)
  meanThrusts = np.array(meanThrusts)
  peakPressures = np.array(peakPressures)
  meanPressures = np.array(meanPressures)
  impulses = np.array(impulses)
  velocities = np.array(velocities)
  
  plt.subplot(3,4,4)
  plt.plot(times, '.')
  plt.title("Computation time")
  plt.xlabel("Simulation #")
  plt.ylabel("Time (s)")
  plt.margins(y=0.1)
  plt.ylim(ymin=0)

  plt.subplot(3,4,5)
  plt.title("Thrust with varying tank filling grade (%)")
  plt.plot(fillingGrades, peakThrusts, 'o', markersize=4)
  plt.plot(fillingGrades, meanThrusts, 'x', markersize=4)
  plt.legend(('Peak', 'Mean'))
  plt.xlabel("Filling grade (%)")
  plt.ylabel("Thrust (kN)")
  plt.margins(y=0.1)
  plt.ylim(ymin=0)

  plt.subplot(3,4,6)
  plt.title("Thrust with varying tank filling temperature")
  plt.plot(initialTemperatures, peakThrusts, 'o', markersize=4)
  plt.plot(initialTemperatures, meanThrusts, 'x', markersize=4)
  plt.legend(('Peak', 'Mean'))
  plt.xlabel("Temperature (°C)")
  plt.ylabel("Thrust (kN)")
  plt.margins(y=0.1)
  plt.ylim(ymin=0)

  plt.subplot(3,4,7)
  plt.title("Pressure with varying tank filling grade (%)")
  plt.plot(fillingGrades, peakPressures, 'o', markersize=4)
  plt.plot(fillingGrades, meanPressures, 'x', markersize=4)
  plt.plot(fillingGrades, peakCCPressures, 'o', markersize=4)
  plt.plot(fillingGrades, meanCCPressures, 'x', markersize=4)
  plt.legend(('Peak (Tank)', 'Mean (Tank)', 'Peak (Chamber)', 'Mean (Chamber)'))
  plt.xlabel("Filling grade (%)")
  plt.ylabel("Pressure (bar)")
  plt.margins(y=0.1)
  plt.ylim(ymin=0)

  plt.subplot(3,4,8)
  plt.title("Pressure with varying tank filling temperature")
  plt.plot(initialTemperatures, peakPressures, 'o', markersize=4)
  plt.plot(initialTemperatures, meanPressures, 'x', markersize=4)
  plt.plot(initialTemperatures, peakCCPressures, 'o', markersize=4)
  plt.plot(initialTemperatures, meanCCPressures, 'x', markersize=4)
  plt.legend(('Peak (Tank)', 'Mean (Tank)', 'Peak (Chamber)', 'Mean (Chamber)'))
  plt.xlabel("Temperature (°C)")
  plt.ylabel("Pressure (bar)")
  plt.margins(y=0.1)
  plt.ylim(ymin=0)

  plt.subplot(3,4,9)
  plt.title("Total impulse with varying tank filling grade (%)")
  plt.plot(fillingGrades, impulses, 'o', markersize=4)
  plt.xlabel("Filling grade (%)")
  plt.ylabel("Impulse (kNs)")
  plt.margins(y=0.1)
  plt.ylim(ymin=0)

  plt.subplot(3,4,10)
  plt.title("Total impulse with varying tank filling temperature")
  plt.plot(initialTemperatures, impulses, 'o', markersize=4)
  plt.xlabel("Temperature (°C)")
  plt.ylabel("Impulse (kNs)")
  plt.margins(y=0.1)
  plt.ylim(ymin=0)


  plt.subplot(3,4,11)
  plt.title("Final velocity with varying tank filling grade (%)")
  plt.plot(fillingGrades, velocities, 'o', markersize=4)
  plt.xlabel("Filling grade (%)")
  plt.ylabel("Velocity (m/s)")
  plt.margins(y=0.1)
  plt.ylim(ymin=0)

  plt.subplot(3,4,12)
  plt.title("Final velocity with varying tank filling temperature")
  plt.plot(initialTemperatures, velocities, 'o', markersize=4)
  plt.xlabel("Temperature (°C)")
  plt.ylabel("Velocity (m/s)")
  plt.margins(y=0.1)
  plt.ylim(ymin=0)

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
  