import matplotlib.pyplot as plt
import numpy as np
import scipy.integrate
import matplotlib.tri as tri

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

def makehist(data, bins = 20):
  plt.hist(data, bins, weights=np.ones_like(data)/float(len(data)), cumulative=False)

with open('./tmp/montecarlo.npy', 'rb') as f:
  # [N] = np.load(f)
  N = 20000

  initialTemperatures = []
  fillingGrades = []
  peakThrusts = []
  meanThrusts = []
  peakPressures = []
  meanPressures = []
  peakCCPressures = []
  meanCCPressures = []
  impulses = []
  accelerations = []
  velocities = []
  altitudes = []
  peakAccelerations = []
  meanAccelerations = []
  finalPortRadii = []
  leavingTowerTimes = []
  burnoutTimes = []
  leavingTowerVelocities = []

  times = []

  try:
    for i in range(N):
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


      _time = np.load(f)

      initialTemperature = models["tank"]["derived"][TankModel.derived_temperature][0] - 273.15
      fillingGrade = models["tank"]["derived"][TankModel.derived_liquidLevel][0]
      peakThrust = np.max(models["nozzle"]["derived"][NozzleModel.derived_thrust]) / 1000
      meanThrust = np.mean(models["nozzle"]["derived"][NozzleModel.derived_thrust]) / 1000
      peakPressure = np.max(models["tank"]["derived"][TankModel.derived_pressure]) / constants.Pressure.bar
      meanPressure = np.mean(models["tank"]["derived"][TankModel.derived_pressure]) / constants.Pressure.bar
      peakCCPressure = np.max(models["combustion"]["state"][CombustionModel.states_pressure]) / constants.Pressure.bar
      meanCCPressure = np.mean(models["combustion"]["state"][CombustionModel.states_pressure]) / constants.Pressure.bar
      velocity = np.linalg.norm([models["flight"]["state"][FlightModel.states_vx], models["flight"]["state"][FlightModel.states_vy], models["flight"]["state"][FlightModel.states_vz]], axis=0)
      altitude = models["flight"]["state"][FlightModel.states_z]
      acceleration = np.linalg.norm([models["flight"]["derived"][FlightModel.derived_ax], models["flight"]["derived"][FlightModel.derived_ay], models["flight"]["derived"][FlightModel.derived_az]], axis=0)
      totalImpulse = scipy.integrate.trapz(models["nozzle"]["derived"][NozzleModel.derived_thrust], t) / 1000
      finalPortRadius = models["combustion"]["state"][CombustionModel.states_portRadius][-1] / constants.Lengths.mm
      leavingTowerIndex = np.where(models["flight"]["derived"][FlightModel.derived_onTower] == 0)
      leavingTowerIndex = leavingTowerIndex[0][0]
      leavingTowerTime = t[leavingTowerIndex]
      leavingTowerVelocity = velocity[leavingTowerIndex]

      outletPhase = np.nan_to_num(models["tank"]["derived"][TankModel.derived_outletPhase])
      
      burnoutIndex = np.where(outletPhase > 0.5)
      burnoutIndex = burnoutIndex[0][0]
      burnoutTime = t[burnoutIndex]

      leavingTowerTimes.append(leavingTowerTime)
      burnoutTimes.append(burnoutTime)
      leavingTowerVelocities.append(leavingTowerVelocity)
      initialTemperatures.append(initialTemperature)
      fillingGrades.append(fillingGrade * 100)
      peakThrusts.append(peakThrust)
      meanThrusts.append(meanThrust)
      peakPressures.append(peakPressure)
      meanPressures.append(meanPressure)
      peakCCPressures.append(peakCCPressure)
      meanCCPressures.append(meanCCPressure)
      impulses.append(totalImpulse)
      altitudes.append(np.max(altitude))
      peakAccelerations.append(np.max(acceleration))
      meanAccelerations.append(np.mean(acceleration))
      finalPortRadii.append(finalPortRadius)
      accelerations.append(np.max(acceleration))
      velocities.append(np.max(velocity))
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
    
  except Exception as e:
    print(e)
    
  plt.subplot(h,w,5)
  plt.title("Thrust for all sampled series")
  plt.xlabel("Time (s)")
  plt.ylabel("Thrust (kN)")
  plt.margins(y=0.1)
  plt.xlim(0, 20)

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
  plt.xlim(0, 20)

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
  altitudes = np.array(altitudes)
  accelerations = np.array(accelerations)
  velocities = np.array(velocities)
  peakAccelerations = np.array(peakAccelerations)
  meanAccelerations = np.array(meanAccelerations)
  finalPortRadii = np.array(finalPortRadii)
  times = np.array(times)

  leavingTowerTimes = np.array(leavingTowerTimes)
  burnoutTimes = np.array(burnoutTimes)
  leavingTowerVelocities = np.array(leavingTowerVelocities)
  
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

  triang = tri.Triangulation(fillingGrades, initialTemperatures)
  interpolator = tri.LinearTriInterpolator(triang, impulses)
  xi = np.linspace(85, 100, 150)
  yi = np.linspace(-5, 25, 150)
  Xi, Yi = np.meshgrid(xi, yi)
  zi = interpolator(Xi, Yi)

  plt.subplot(h,w,6)
  plt.title("Total impulse (kNs)")
  CS = plt.contour(Xi, Yi, zi, levels=10, linewidths=0.5, colors='k')
  plt.contourf(Xi, Yi, zi, levels=10, cmap="RdBu_r")
  plt.colorbar()
  plt.xlabel("Filling grade (%)")
  plt.ylabel("Temperature (°C)")

  # plt.subplot(h,w,6)
  # plt.title("Total impulse with varying tank filling grade (%)")
  # plt.plot(fillingGrades, impulses, 'o', markersize=4)
  # plt.xlabel("Filling grade (%)")
  # plt.ylabel("Impulse (kNs)")
  # plt.margins(y=0.1)
  # plt.ylim(ymin=0)

  triang = tri.Triangulation(fillingGrades, initialTemperatures)
  interpolator = tri.LinearTriInterpolator(triang, altitudes / 1000)
  xi = np.linspace(85, 100, 150)
  yi = np.linspace(-5, 25, 150)
  Xi, Yi = np.meshgrid(xi, yi)
  zi = interpolator(Xi, Yi)
  
  plt.subplot(h,w,7)
  plt.title("Altitude (km)")
  CS = plt.contour(Xi, Yi, zi, levels=10, linewidths=0.5, colors='k')
  plt.contourf(Xi, Yi, zi, levels=10, cmap="RdBu_r")
  plt.colorbar()
  plt.xlabel("Filling grade (%)")
  plt.ylabel("Temperature (°C)")


  plt.subplot(h,w,8)
  makehist(altitudes/1000, bins=20)
  plt.ylabel("Probability")
  plt.xlabel("Altitude (km)")

  plt.subplot(h,w,9)
  makehist(finalPortRadii, bins=20)
  plt.ylabel("Probability")
  plt.xlabel("Final port radius (mm)")
#  plt.subplot(h,w,8)
#  plt.title("Final velocity with varying tank filling grade (%)")
#  plt.plot(fillingGrades, altitudes, 'o', markersize=4)
#  plt.xlabel("Filling grade (%)")
#  plt.ylabel("Velocity (m/s)")
#  plt.margins(y=0.1)
#  plt.ylim(ymin=0)
#

#  plt.subplot(h,w,9)
#  plt.title("Final velocity with varying tank filling temperature")
#  plt.plot(initialTemperatures, altitudes, 'o', markersize=4)
#  plt.xlabel("Temperature (°C)")
#  plt.ylabel("Velocity (m/s)")
#  plt.margins(y=0.1)
#  plt.ylim(ymin=0)

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



index = 1
width = 4
height = 4

def nextSubplot():
  global index
  print(height, width, index)
  plt.subplot(height, width, index)
  index = (index) % (width * height) + 1

nextSubplot()
makehist(finalPortRadii)
plt.ylabel("Probability")
plt.xlabel("Final port radius (mm)")

nextSubplot()
makehist(peakThrusts)
plt.ylabel("Probability")
plt.xlabel("Peak thrusts (kN)")

nextSubplot()
makehist(meanThrusts)
plt.ylabel("Probability")
plt.xlabel("Mean thrusts (kN)")

nextSubplot()
makehist(peakPressures)
plt.ylabel("Probability")
plt.xlabel("Peak tank pressure (kN)")

nextSubplot()
makehist(meanPressures)
plt.ylabel("Probability")
plt.xlabel("Mean tank pressure (kN)")

nextSubplot()
makehist(peakCCPressures)
plt.ylabel("Probability")
plt.xlabel("Peak chamber pressure (kN)")

nextSubplot()
makehist(meanCCPressures)
plt.ylabel("Probability")
plt.xlabel("Mean chamber pressure (kN)")

nextSubplot()
makehist(impulses)
plt.ylabel("Probability")
plt.xlabel("Impulse (kNs)")

nextSubplot()
makehist(altitudes / 1e3)
plt.ylabel("Probability")
plt.xlabel("Altitude (km)")

nextSubplot()
makehist(velocities)
plt.ylabel("Probability")
plt.xlabel("Max velocity (m/s)")

nextSubplot()
makehist(peakAccelerations / 9.81066)
plt.ylabel("Probability")
plt.xlabel("Peak acceleration (G)")

nextSubplot()
makehist(meanAccelerations / 9.81066)
plt.ylabel("Probability")
plt.xlabel("Mean acceleration (G)")

nextSubplot()
makehist(times)
plt.ylabel("Probability")
plt.xlabel("Simulation time (s)")

nextSubplot()
makehist(leavingTowerTimes, bins=6)
plt.ylabel("Probability")
plt.xlabel("Tower leaving time (s)")

nextSubplot()
makehist(burnoutTimes, bins=6)
plt.ylabel("Probability")
plt.xlabel("Burnout time (s)")

nextSubplot()
makehist(leavingTowerVelocities, bins=6)
plt.ylabel("Probability")
plt.xlabel("Velocity when leaving tower (m/s)")

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
