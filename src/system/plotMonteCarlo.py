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

h = 3
w = 5

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
  peakAccelerations = []
  meanAccelerations = []
  finalPortRadii = []


  times = []

  try:
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

      [_time] = np.load(f)

      initialTemperature = models["tank"]["derived"][TankModel.derived_temperature][0] - 273.15
      fillingGrade = models["tank"]["derived"][TankModel.derived_liquidLevel][0]
      peakThrust = np.max(models["nozzle"]["derived"][NozzleModel.derived_thrust]) / 1000
      meanThrust = np.mean(models["nozzle"]["derived"][NozzleModel.derived_thrust]) / 1000
      peakPressure = np.max(models["tank"]["derived"][TankModel.derived_pressure]) / constants.Pressure.bar
      meanPressure = np.mean(models["tank"]["derived"][TankModel.derived_pressure]) / constants.Pressure.bar
      peakCCPressure = np.max(models["combustion"]["state"][CombustionModel.states_pressure]) / constants.Pressure.bar
      meanCCPressure = np.mean(models["combustion"]["state"][CombustionModel.states_pressure]) / constants.Pressure.bar
      velocity = np.linalg.norm([models["flight"]["state"][FlightModel.states_vx], models["flight"]["state"][FlightModel.states_vy], models["flight"]["state"][FlightModel.states_vz]], axis=0)

      acceleration = np.linalg.norm([models["flight"]["derived"][FlightModel.derived_ax], models["flight"]["derived"][FlightModel.derived_ay], models["flight"]["derived"][FlightModel.derived_az]], axis=0)

      totalImpulse = scipy.integrate.simps(models["nozzle"]["derived"][NozzleModel.derived_thrust], t) / 1000

      finalPortRadius = models["combustion"]["state"][CombustionModel.states_portRadius][-1] / constants.Lengths.mm

      initialTemperatures.append(initialTemperature)
      fillingGrades.append(fillingGrade * 100)
      peakThrusts.append(peakThrust)
      meanThrusts.append(meanThrust)
      peakPressures.append(peakPressure)
      meanPressures.append(meanPressure)
      peakCCPressures.append(peakCCPressure)
      meanCCPressures.append(meanCCPressure)
      impulses.append(totalImpulse)
      velocities.append(velocity[-1])
      peakAccelerations.append(np.max(acceleration))
      meanAccelerations.append(np.mean(acceleration))
      finalPortRadii.append(finalPortRadius)
      times.append(_time)

      plt.subplot(h,w,5)
      plt.plot(t, models["nozzle"]["derived"][NozzleModel.derived_thrust] / 1000, '-', linewidth=0.5)

      plt.subplot(h,w,10)
      plt.plot(t, models["flight"]["state"][FlightModel.states_z], '-', linewidth=0.5)

      plt.subplot(h,w,15)
      plt.plot(t, models["tank"]["derived"][TankModel.derived_pressure] / constants.Pressure.bar, '-', linewidth=0.5)

      # plt.subplot(h,w,14)
      # plt.plot(t, models["tank"]["derived"][TankModel.derived_temperature] - 273.15, '-', linewidth=0.5)
    
      # plt.subplot(h,w,14)
      # plt.plot(t, models["combustion"]["state"][CombustionModel.states_portRadius] / constants.Lengths.mm, '-', linewidth=0.5)
    
  except Exception:
    print("oh no")
    
  plt.subplot(h,w,5)
  plt.title("Thrust for all sampled series")
  plt.xlabel("Time (s)")
  plt.ylabel("Thrust (kN)")
  plt.margins(y=0.1)

  plt.subplot(h,w,10)
  plt.title("Altitude for all sampled series")
  plt.xlabel("Time (s)")
  plt.ylabel("Altitude (m)")
  plt.margins(y=0.1)

  plt.subplot(h,w,15)
  plt.title("Tank pressure for all sampled series")
  plt.xlabel("Time (s)")
  plt.ylabel("Tank pressure (bar)")
  plt.margins(y=0.1)

  # plt.subplot(h,w,14)
  # plt.title("Tank temperature for all sampled series")
  # plt.xlabel("Time (s)")
  # plt.ylabel("Tank temperature (°C)")
  # plt.margins(y=0.1)

  initialTemperatures = np.array(initialTemperatures)
  fillingGrades = np.array(fillingGrades)
  peakThrusts = np.array(peakThrusts)
  meanThrusts = np.array(meanThrusts)
  peakPressures = np.array(peakPressures)
  meanPressures = np.array(meanPressures)
  impulses = np.array(impulses)
  velocities = np.array(velocities)
  peakAccelerations = np.array(peakAccelerations)
  meanAccelerations = np.array(meanAccelerations)
  finalPortRadii = np.array(finalPortRadii)
  
  # plt.subplot(h,w,4)
  # plt.plot(times, '.')
  # plt.title("Computation time")
  # plt.xlabel("Simulation #")
  # plt.ylabel("Time (s)")
  # plt.margins(y=0.1)
  # plt.ylim(ymin=0)

  plt.subplot(h,w,1)
  plt.title("Thrust with varying tank filling grade (%)")
  plt.plot(fillingGrades, peakThrusts, 'o', markersize=4)
  plt.plot(fillingGrades, meanThrusts, 'x', markersize=4)
  plt.legend(('Peak', 'Mean'))
  plt.xlabel("Filling grade (%)")
  plt.ylabel("Thrust (kN)")
  plt.margins(y=0.1)
  plt.ylim(ymin=0)

  plt.subplot(h,w,2)
  plt.title("Thrust with varying tank filling temperature")
  plt.plot(initialTemperatures, peakThrusts, 'o', markersize=4)
  plt.plot(initialTemperatures, meanThrusts, 'x', markersize=4)
  plt.legend(('Peak', 'Mean'))
  plt.xlabel("Temperature (°C)")
  plt.ylabel("Thrust (kN)")
  plt.margins(y=0.1)
  plt.ylim(ymin=0)

  plt.subplot(h,w,3)
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

  plt.subplot(h,w,4)
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

  plt.subplot(h,w,6)
  plt.title("Total impulse with varying tank filling grade (%)")
  plt.plot(fillingGrades, impulses, 'o', markersize=4)
  plt.xlabel("Filling grade (%)")
  plt.ylabel("Impulse (kNs)")
  plt.margins(y=0.1)
  plt.ylim(ymin=0)

  plt.subplot(h,w,7)
  plt.title("Total impulse with varying tank filling temperature")
  plt.plot(initialTemperatures, impulses, 'o', markersize=4)
  plt.xlabel("Temperature (°C)")
  plt.ylabel("Impulse (kNs)")
  plt.margins(y=0.1)
  plt.ylim(ymin=0)


  plt.subplot(h,w,8)
  plt.title("Final velocity with varying tank filling grade (%)")
  plt.plot(fillingGrades, velocities, 'o', markersize=4)
  plt.xlabel("Filling grade (%)")
  plt.ylabel("Velocity (m/s)")
  plt.margins(y=0.1)
  plt.ylim(ymin=0)

  plt.subplot(h,w,9)
  plt.title("Final velocity with varying tank filling temperature")
  plt.plot(initialTemperatures, velocities, 'o', markersize=4)
  plt.xlabel("Temperature (°C)")
  plt.ylabel("Velocity (m/s)")
  plt.margins(y=0.1)
  plt.ylim(ymin=0)

  plt.subplot(h,w,11)
  plt.title("Acceleration with varying tank filling grade (%)")
  plt.plot(fillingGrades, peakAccelerations, 'o', markersize=4)
  plt.plot(fillingGrades, meanAccelerations, 'x', markersize=4)
  plt.legend(('Peak', 'Mean'))
  plt.xlabel("Filling grade (%)")
  plt.ylabel("Acceleration (m s-2)")
  plt.margins(y=0.1)
  plt.ylim(ymin=0)

  plt.subplot(h,w,12)
  plt.title("Acceleration with varying tank filling temperature")
  plt.plot(initialTemperatures, peakAccelerations, 'o', markersize=4)
  plt.plot(initialTemperatures, meanAccelerations, 'x', markersize=4)
  plt.legend(('Peak', 'Mean'))
  plt.xlabel("Temperature (°C)")
  plt.ylabel("Acceleration (m s-2)")
  plt.margins(y=0.1)
  plt.ylim(ymin=0)
  
  plt.subplot(h,w,13)
  plt.plot(fillingGrades, finalPortRadii, 'o', markersize=4)
  plt.title("Port radius with varying tank filling grade (%)")
  plt.xlabel("Filling grade (%)")
  plt.ylabel("Port radius (mm)")
  plt.hlines(75 - 2, np.min(fillingGrades), np.max(fillingGrades), 'r', linestyles="dashed", label="Casing")
  plt.hlines(75 - 2 - 10, np.min(fillingGrades), np.max(fillingGrades), 'g', linestyles="dashdot", label="Margin")
  plt.legend()
  plt.margins(y=0.1)
  plt.ylim(ymin=0)

  plt.subplot(h,w,14)
  plt.plot(initialTemperatures, finalPortRadii, 'o', markersize=4)
  plt.title("Port radius with varying tank filling temperature")
  plt.xlabel("Temperature (°C)")
  plt.ylabel("Port radius (mm)")
  plt.hlines(75 - 2, np.min(initialTemperatures), np.max(initialTemperatures), 'r', linestyles="dashed", label="Casing")
  plt.hlines(75 - 2 - 10, np.min(initialTemperatures), np.max(initialTemperatures), 'g', linestyles="dashdot", label="Margin")
  plt.legend()
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
  