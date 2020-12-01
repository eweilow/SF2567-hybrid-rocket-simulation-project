import matplotlib.pyplot as plt
import numpy as np
import scipy.interpolate

from models.tank import TankModel
from models.injector import InjectorModel
from models.combustion import CombustionModel
from models.nozzle import NozzleModel
from models.passiveVent import PassiveVentModel
from models.environment import EnvironmentModel
from models.flight import FlightModel
import utils.constants as constants

with open('./tmp/benchmark.npy', 'rb') as f:
  data = np.load(f)
  baselineRange = np.load(f)
  baselineTimes = np.load(f, allow_pickle=True)
  baselineOutput = np.load(f, allow_pickle=True)
  rest = np.load(f, allow_pickle=True)
  
  baselineOutput = baselineOutput.tolist()

  errors = {}


  for key in ["flight", "combustion"]:
    state = baselineOutput[key]["state"]
    baselineTimes, unique_indices = np.unique(baselineTimes, return_index=True)
    state = state[:,unique_indices]
    normalInterpolator = scipy.interpolate.interp1d(baselineTimes, state, axis=1, kind=3)

    [start, mid, end] = baselineRange

    endTime = end
    if key == "combustion":
      endTime = mid
      
    index = np.searchsorted(baselineTimes, endTime)

    time = baselineTimes[1:index]
    baseline = normalInterpolator(time)
    baselineNorm = np.linalg.norm(baseline, axis=0)

    errors[key] = {
      "time": time,
      "baseline": baseline,
      "baselineNorm": baselineNorm
    }

  rows = []
  groupings = {}

  tolKeys = []

  try:
    while True:
      tols, interpolationEnabled, enableDAESolver, rootFindingType = np.load(f, allow_pickle=True)

      timeRanges = np.load(f)
      times = np.load(f, allow_pickle=True)
      modelsOutput = np.load(f, allow_pickle=True)
      rest = np.load(f, allow_pickle=True)

      nfev, njev, nlu, firstSolverTime, nfev2, njev2, nlu2, secondSolverTime = rest

      tols = float(tols)
      enableDAESolver = False if enableDAESolver == "False" or enableDAESolver == False or enableDAESolver == "0.0" or enableDAESolver == 0.0 else True
      interpolationEnabled = False if interpolationEnabled == "False" or interpolationEnabled == False or interpolationEnabled == "0.0" or interpolationEnabled == 0.0 else True
      rootFindingType = False if rootFindingType == False or rootFindingType == "0.0" or rootFindingType == 0.0 else rootFindingType

      if interpolationEnabled:
        groupBy = "{:}{:}".format(
          "DAE (IDA)" if enableDAESolver else "ODE (LSODI)",
          "\n({:} iterative\nroot finding)".format(rootFindingType) if not rootFindingType == False else "\n(algebraic variables\nroot finding)"
        )

        if not groupBy in groupings:
          groupings[groupBy] = {}

        tolKey = "{:.0e}".format(tols)
        if not tolKey in groupings[groupBy]:
          groupings[groupBy][tolKey] = {
            "time": float(firstSolverTime),
            "nfev": int(nfev)
          }

        if not tolKey in tolKeys:
          tolKeys.append(tolKey)
      
      rows.append((
        tols,
        enableDAESolver,
        interpolationEnabled,
        rootFindingType,
        timeRanges,
        times,
        modelsOutput,
        rest
      ))

  except IOError:
    pass
  
  plotBars = False

  if plotBars:
    groupCount = len(tolKeys)
    wantedWidth = 0.8
    barWidth = wantedWidth/(groupCount)
    groupOffsets = -np.linspace(-wantedWidth/2 + barWidth/2, wantedWidth/2 - barWidth/2, groupCount)
    
    plt.subplot(1,2,1)
    i = 0
    for tolKey in tolKeys:
      indices = np.arange(1, len(groupings.keys()) + 1, 1)
      values = []
      labels = []

      for groupBy in groupings:
        values.append(groupings[groupBy][tolKey]["nfev"])
        labels.append(groupBy)
      plt.barh(indices + groupOffsets[i], values, zorder=5, height=barWidth)
      if i == 0:
        plt.yticks(indices, labels)
      i += 1

    from matplotlib.ticker import ScalarFormatter
    
    plt.xscale("log")
    plt.gca().xaxis.set_major_formatter(ScalarFormatter())
    plt.xlim(left=900)
    plt.grid(axis="x", which="minor")
    plt.grid(axis="x", which="major")
    plt.xlabel("Number of function evaluations\nper simulation")
    plt.legend(["rtol = {:}".format(val) for val in tolKeys], bbox_to_anchor=(0.5, 1.2), loc="upper center", ncol=3)

    plt.subplot(1,2,2)
    i = 0
    for tolKey in tolKeys:
      indices = np.arange(1, len(groupings.keys()) + 1, 1)
      values = []
      labels = []

      for groupBy in groupings:
        values.append(groupings[groupBy][tolKey]["time"])
        labels.append(groupBy)
      plt.barh(indices + groupOffsets[i], values, zorder=5, height=barWidth)
      if i == 0:
        plt.yticks(indices, labels)
      i += 1

    from matplotlib.ticker import ScalarFormatter
    
    plt.xscale("log")
    plt.gca().xaxis.set_major_formatter(ScalarFormatter())
    plt.xlim(left=9)
    plt.grid(axis="x", which="minor")
    plt.grid(axis="x", which="major")
    plt.xlabel("Simulation time (s)")
    plt.legend(["rtol = {:}".format(val) for val in tolKeys], bbox_to_anchor=(0.5, 1.2), loc="upper center", ncol=3)

    plt.tight_layout()
    plt.subplots_adjust(
      left=0.25,
      bottom=0.15,
      top=0.386,
      right=0.536,
      wspace=0.93,
      hspace=0.15
    )
    plt.show()
  else:


    def plotWithFilter(title, errorKey, filterFn):
      for row in rows:
        tols, enableDAESolver, interpolationEnabled, rootFindingType, timeRanges, times, modelsOutput, rest = row

        if not filterFn(tols, enableDAESolver, interpolationEnabled, rootFindingType):
          continue
        
        name = "tol={:.0e}".format(tols)

        modelsOutput = modelsOutput.tolist()
        state = modelsOutput[errorKey]["state"]
        times, unique_indices = np.unique(times, return_index=True)
        state = state[:,unique_indices]
        interpolator = scipy.interpolate.interp1d(times, state, axis=1, kind=3)

        time = errors[errorKey]["time"]
        baseline = errors[errorKey]["baseline"]
        baselineNorm = errors[errorKey]["baselineNorm"]

        values = interpolator(time)

        error = (values - baseline) / baselineNorm
        errorNorm = np.linalg.norm(error, axis=0)

        plt.title(title)
        plt.semilogy(time, errorNorm, linewidth=1, label=name)
        plt.ylabel("relative error")
        plt.xlabel("time (s)")
        plt.grid(axis="y", zorder=-5)
        plt.ylim(top=1, bottom=1e-8)
        plt.legend()     

    plt.subplot(2,2,1)
    plotWithFilter("DAE + ODE (combustion phase)", "combustion", lambda tols, enableDAESolver, interpolationEnabled, rootFindingType: enableDAESolver and interpolationEnabled and rootFindingType == False)
    plt.subplot(2,2,2)
    plotWithFilter("DAE + ODE (flight phase)", "flight", lambda tols, enableDAESolver, interpolationEnabled, rootFindingType: enableDAESolver and interpolationEnabled and rootFindingType == "brentq")
    plt.subplot(2,2,3)
    plotWithFilter("ODE + ODE (combustion phase)", "combustion", lambda tols, enableDAESolver, interpolationEnabled, rootFindingType: not enableDAESolver and interpolationEnabled and rootFindingType == "brentq")
    plt.subplot(2,2,4)
    plotWithFilter("ODE + ODE (flight phase)", "flight", lambda tols, enableDAESolver, interpolationEnabled, rootFindingType: not enableDAESolver and not interpolationEnabled and rootFindingType == "brentq")
  
    plt.show()