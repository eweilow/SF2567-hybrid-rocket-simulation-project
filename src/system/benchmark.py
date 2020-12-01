from runSingle import runSystem
import options
import os, sys
import numpy as np

import utils.cea

hiddenPrints = utils.cea.HiddenPrints()

options.enableCeaLookup = True
options.enableTankTemperatureInterpolation = True
options.enableTankWallHeatTransfer = True
options.printFailedInjector = False
options.printTime = False

options.solveTankVolumeContraintWithDAE = True
options.rootFindingType = "brentq"
options.enableDAESolver = True

options.combustion_rtol = 1e-6
options.combustion_atol = 1e-4
options.flight_rtol = 1e-6
options.flight_atol = 1e-4

a, baselineRange, baselineTimes, baselineOutput, *rest = runSystem()
with open("/data/benchmark.npy", 'wb') as f:
  np.save(f, (1e-6, True, True, False))
  np.save(f, baselineRange)
  np.save(f, baselineTimes)
  np.save(f, baselineOutput)
  np.save(f, rest)

for repetitions in range(1):
  for tols in [1e-3, 1e-4, 1e-5]:
    options.combustion_rtol = tols
    options.combustion_atol = 1e-4
    options.flight_rtol = tols
    options.flight_atol = 1e-4

    for interpolationEnabled in [True, False]:
      options.enableCeaLookup = interpolationEnabled
      options.enableTankTemperatureInterpolation = interpolationEnabled

      for enableDAESolver in [True, False]:
        options.enableDAESolver = enableDAESolver
        for rootFindingType in ["brentq", "bisect", False]:
          if rootFindingType is False:
            if not enableDAESolver:
              continue # unsupported

            options.solveTankVolumeContraintWithDAE = True
          else:
            options.solveTankVolumeContraintWithDAE = False
            options.rootFindingType = rootFindingType

          utils.cea.NasaCEA.clear()

          with hiddenPrints:
            cpuTime, timeRanges, times, modelsOutput, *rest = runSystem()
            nfev, njev, nlu, firstSolverTime, nfev2, njev2, nlu2, secondSolverTime = rest

          print(
            "tol={:.2e}, interpolationEnabled={:}, enableDAESolver={:}, rootFindingType={:}, nfev = {:}, njev = {:}, nlu = {:}, firstSolverTime = {:.2f}, nfev2 = {:}, njev2 = {:}, nlu2 = {:}, secondSolverTime = {:.2f}".format(
              tols,
              interpolationEnabled,
              enableDAESolver,
              rootFindingType, nfev, njev, nlu, firstSolverTime, nfev2, njev2, nlu2, secondSolverTime
            )
          )
          with open("/data/benchmark.npy", 'ab') as f:
            np.save(f, (tols, interpolationEnabled, enableDAESolver, rootFindingType))
            np.save(f, timeRanges)
            np.save(f, times)
            np.save(f, modelsOutput)
            np.save(f, rest)
