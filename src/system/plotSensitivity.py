import matplotlib.pyplot as plt
import numpy as np
import scipy.integrate
import matplotlib.tri as tri
from scipy import interpolate, stats

from models.tank import TankModel
from models.injector import InjectorModel
from models.combustion import CombustionModel
from models.nozzle import NozzleModel
from models.passiveVent import PassiveVentModel
from models.environment import EnvironmentModel
from models.flight import FlightModel
import utils.constants as constants

def interpolateToFindPeak(x, y):
#  x, unique_indices = np.unique(x, return_index=True)
#  y = y[unique_indices]
#
#  peakIndex = np.argmax(y)
#  peakY = y[peakIndex]
#
#  interpolator = interpolate.interp1d(x, y, kind=5)
#
#  xvaluesLeft = max(x[0], x[peakIndex - 5])
#  xvaluesRight = min(x[-1], x[peakIndex + 5])
#
#  xv = np.linspace(xvaluesLeft, xvaluesRight, 1000)
#
#  newBest = np.max(interpolator(xv))
#
  return np.max(y)
  

unphysicalModesToSkip = ["nozzleEfficiency", "combustionEfficiency", "tankInsideRadius", "tankTopVentLiquidGasTransferRadius", "tankOutletLiquidGasTransferRadius"]

with open('./tmp/sensitivity.npy', 'rb') as f:
  [N] = np.load(f)

  modeDatas = {}
  legends = tuple()

  modes = []
  
  try:
    for i in range(N):
      mode = np.load(f)
      baseline = np.load(f)
      currentValue = np.load(f)
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
      time = np.load(f)

      mode = str(mode)
      if mode in unphysicalModesToSkip:
        continue
      
      import assumptions
      if not mode in assumptions.availableToRandomize:
        continue
      
      initialTemperature = models["tank"]["derived"][TankModel.derived_temperature][0] - 273.15
      fillingGrade = models["tank"]["derived"][TankModel.derived_liquidLevel][0]
      peakThrust = interpolateToFindPeak(t, models["nozzle"]["derived"][NozzleModel.derived_thrust]) / 1000
      meanThrust = np.mean(models["nozzle"]["derived"][NozzleModel.derived_thrust]) / 1000
      peakPressure = interpolateToFindPeak(t, models["tank"]["derived"][TankModel.derived_pressure]) / constants.Pressure.bar
      meanPressure = np.mean(models["tank"]["derived"][TankModel.derived_pressure]) / constants.Pressure.bar
      peakCCPressure = interpolateToFindPeak(t, models["combustion"]["state"][CombustionModel.states_pressure]) / constants.Pressure.bar
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

      if not mode in modeDatas:
        modeDatas[mode] = {
          "values": [],
          "baselines": [],
          "burnoutTimes": [],
          "leavingTowerTimes": [],
          "leavingTowerVelocity": [],
          "altitudes": [],
          "peakThrusts": [],
          "finalPortRadii": [],
          "impulse": []
        }
        legends = (*legends, mode)
        modes.append(mode)

      modeDatas[mode]["values"].append(currentValue)
      modeDatas[mode]["baselines"].append(baseline)
      modeDatas[mode]["burnoutTimes"].append(burnoutTime)
      modeDatas[mode]["leavingTowerTimes"].append(leavingTowerTime)
      modeDatas[mode]["leavingTowerVelocity"].append(leavingTowerVelocity)
      modeDatas[mode]["altitudes"].append(interpolateToFindPeak(t, altitude) / 1e3)
      modeDatas[mode]["peakThrusts"].append(peakThrust)
      modeDatas[mode]["finalPortRadii"].append(finalPortRadius)
      modeDatas[mode]["impulse"].append(totalImpulse)


  except Exception as e:
    print(e)

  for mode in modeDatas:
    for key in modeDatas[mode]:
      modeDatas[mode][key] = np.array(modeDatas[mode][key])
    modeDatas[mode]["deviation"] = (modeDatas[mode]["values"] - modeDatas[mode]["baselines"]) / modeDatas[mode]["baselines"] * 100
     
  plotBars = True
  plotGrouped = True

  plotBarPercentage = False

  index = 1
  width = 3 if plotBars else 6
  height = 2 if plotBars else 1
  
  def nextSubplot():
    if plotGrouped:
      return

    global index
    ax = plt.subplot(height, width, index)
    index = (index) % (width * height) + 1

  plt.figure(figsize=(12, 8))
  plotLabels = modes


  groupIndex = 0
  groupCount = 6
  groupLegend = []
  wantedWidth = 0.8 if plotGrouped else 0.8
  barWidth = wantedWidth/(groupCount) if plotGrouped else wantedWidth
  groupOffsets = -np.linspace(-wantedWidth/2 + barWidth/2, wantedWidth/2 - barWidth/2, groupCount)

  def plotSeries(series, unit, title, thresholds = None):
    global groupIndex

    slopes = []
    stds = []
    tickLabels = []
    numIndices = 0

    centerValues = []
    for key in plotLabels:
      # https://stackoverflow.com/a/37415568
      x, y = zip(*sorted(zip(modeDatas[key]["deviation"], modeDatas[key][series])))
      x = np.array(x)
      y = np.array(y)

      try:
        centerValues.append(interpolate.interp1d(x, y)(0))
      except Exception:
        pass

    centerValue = np.mean(centerValues)
      
    for key in plotLabels:
      # https://stackoverflow.com/a/37415568
      if not plotBars:
        x, y = zip(*sorted(zip(modeDatas[key]["deviation"], modeDatas[key][series])))
        x = np.array(x)
        y = np.array(y)
        indices = np.argwhere(np.abs(x) <= 10)
        x = x[indices]
        y = y[indices]


      else:
        if plotBarPercentage:
          x, y = zip(*sorted(zip(modeDatas[key]["deviation"], (modeDatas[key][series] - centerValue) / centerValue)))
          x = np.array(x)
          y = np.array(y) * 100 * 10
        else:
          x, y = zip(*sorted(zip(modeDatas[key]["deviation"], modeDatas[key][series] - centerValue)))
          x = np.array(x)
          y = np.array(y) * 10

        slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)

        # if not plotBars and slope > -5 and slope < 5:
        #   continue

        slopes.append(slope)
        stds.append(std_err)

      legendUnit = assumptions.availableToRandomize[key]["unit"] if "unit" in assumptions.availableToRandomize[key] else None
      if legendUnit is None:
        tickLabels.append(assumptions.availableToRandomize[key]["title"])
      else:
        scale = assumptions.availableToRandomize[key]["unitScale"] if "unitScale" in assumptions.availableToRandomize[key] else 1

        zeroValue = modeDatas[key]["baselines"][0] * scale
        baselineValue = zeroValue * 0.1
        if legendUnit == None:
          tenPercentValue = "{:.2f}".format(baselineValue)
          zeroPercentValue = "{:.2f}".format(zeroValue)
        elif legendUnit[0] == "°":
          tenPercentValue = "{:.2f}{:}".format(baselineValue, legendUnit)
          zeroPercentValue = "{:.2f}{:}".format(zeroValue, legendUnit)
        else:
          tenPercentValue = "{:.2f} {:}".format(baselineValue, legendUnit)
          zeroPercentValue = "{:.2f} {:}".format(zeroValue, legendUnit)

        tickLabels.append("{:} ({:})".format(
          assumptions.availableToRandomize[key]["title"],
          zeroPercentValue,
        ))
      numIndices += 1

      if not plotBars:
        # plt.plot(x, intercept + slope * x, '--')
        plt.plot(x, y, '.-', zorder=5,linewidth=1,markersize=2)
    
    if not plotBars:
      # plt.legend(legends)
      plt.ylabel("{:} ({:})".format(title, unit))
      plt.xlabel("% change")
      plt.legend(tickLabels)
      plt.xlim(-10, 10)
      plt.grid(zorder=-5)
      return

    indices = np.arange(numIndices)
    labels, x, err = zip(*sorted(zip(tickLabels, slopes, stds)))
    # x, err, labels = (slopes, stds, tickLabels)

    x = np.array(x)
    err = np.array(err)
    
  
    if not plotGrouped:
      plt.vlines(0, ymin=-0.75, ymax=numIndices - 0.25, color="#000000")
      plt.text(0,numIndices,"{:.2f} {:}".format(centerValue, unit), horizontalalignment="center", verticalalignment="bottom")

      if not thresholds is None:
        for thrsTitle in thresholds:
          if plotBarPercentage:
            xvalue = thresholds[thrsTitle] / centerValue * 100 - 100
          else:
            xvalue = thresholds[thrsTitle] - centerValue
          plt.vlines(xvalue, ymin=-0.75, ymax=numIndices - 0.25, color="#e377c2", linestyle="dashed",zorder=5)
          plt.text(xvalue,numIndices,thrsTitle, horizontalalignment="center", verticalalignment="bottom",zorder=5)

      plt.ylim(-1, numIndices + 0.5)
    
    barlocation = indices
    if plotGrouped:
      barlocation = barlocation + groupOffsets[groupIndex]
      groupIndex += 1
        
    plt.barh(barlocation, np.abs(x) if plotGrouped else x, xerr=err,zorder=5,height=barWidth)
    plt.yticks(indices, labels)

    if plotGrouped:
      plt.xlabel("relative (abs-value) change per 10%")
      groupLegend.append("{:} (10% = {:.2f} {:})".format(title, centerValue * 0.1, unit))
    else:
      if plotBarPercentage:
        plt.xlabel("{:} (% change per 10%)".format(title))
      else:
        plt.xlabel("{:} ({:} change per 10%)".format(title, unit))

    if plotGrouped and groupIndex == groupCount:
      left, right = plt.xlim()
      right *= 1.05
      plt.barh(indices, np.ones(len(barlocation)) * right, zorder=-5,height=wantedWidth, color="#efefef")
      plt.xscale("log")
      plt.xlim(left=0.1, right=right)
    else:
      plt.grid(axis="x",zorder=-5)

    return x, err, labels
    # plt.xlabel("Relative deviation (%)")
    # plt.ylabel("Max altitude (km)")
  
  nextSubplot()
  plotSeries("altitudes", "km", "Max altitude", { "10 km": 10, "14 km": 14 })
  nextSubplot()
  plotSeries("burnoutTimes", "s", "Burnout time")
  nextSubplot()
  plotSeries("peakThrusts", "kN", "Peak thrust")
  nextSubplot()
  plotSeries("leavingTowerVelocity", "m/s", "End of tower velocity")
  nextSubplot()
  plotSeries("finalPortRadii", "mm", "Final port radius", { "max tolerable": 70 })
  nextSubplot()
  plotSeries("impulse", "kNs", "Total impulse")

  if plotGrouped:
    plt.xlim(left=0)
    plt.legend(groupLegend)
    from matplotlib.ticker import PercentFormatter

    plt.gca().xaxis.set_major_formatter(PercentFormatter())
  
#
#  nextSubplot()
#  for key in plotLabels:
#    x, y = zip(*sorted(zip(modeDatas[key]["deviation"], modeDatas[key]["burnoutTimes"])))
#    plt.plot(x, y, '.-')
#  plt.xlabel("Relative deviation (%)")
#  plt.ylabel("Burnout time (s)")
#  plt.legend(legends)
#
#  nextSubplot()
#  for key in plotLabels:
#    x, y = zip(*sorted(zip(modeDatas[key]["deviation"], modeDatas[key]["peakThrusts"])))
#    plt.plot(x, y, '.-')
#  plt.xlabel("Relative deviation (%)")
#  plt.ylabel("Peak thrust (kN)")
#  plt.legend(legends)
#
#  nextSubplot()
#  for key in plotLabels:
#    x, y = zip(*sorted(zip(modeDatas[key]["deviation"], modeDatas[key]["leavingTowerVelocity"])))
#    plt.plot(x, y, '.-')
#  plt.xlabel("Relative deviation (%)")
#  plt.ylabel("Velocity leaving tower (m/s)")
#  plt.legend(legends)
#  
#  nextSubplot()
#  for key in plotLabels:
#    x, y = zip(*sorted(zip(modeDatas[key]["deviation"], modeDatas[key]["finalPortRadii"])))
#    plt.plot(x, y, '.-')
#  plt.xlabel("Relative deviation (%)")
#  plt.ylabel("Final port radius (mm)")
#  plt.legend(legends)
#  
#  nextSubplot()
#  for key in plotLabels:
#    x, y = zip(*sorted(zip(modeDatas[key]["deviation"], modeDatas[key]["impulse"])))
#    plt.plot(x, y, '.-')
#  plt.xlabel("Relative deviation (%)")
#  plt.ylabel("Total impulse (kNs)")
#  plt.legend(legends)

  if plotBars:
    if plotGrouped:
      plt.tight_layout()
      plt.subplots_adjust(
        left=0.25,
        bottom=0.05,
        top=1-0.05,
        right=1-0.025,
        wspace=0.2,
        hspace=0.15
      )
    else:
      plt.tight_layout()
      plt.subplots_adjust(
        left=0.115,
        bottom=0.05,
        top=1-0.025,
        right=1-0.025,
        wspace=0.55,
        hspace=0.15
      )
  else:
    plt.tight_layout()
    plt.subplots_adjust(
      left=0.05,
      bottom=0.05,
      top=1-0.05,
      right=1-0.025,
      wspace=0.2,
      hspace=0.15
    )
  plt.show()