from runSingle import runSystem
import options
import os, sys

import utils.cea

hiddenPrints = utils.cea.HiddenPrints()


options.enableCeaLookup = False
options.enableTankTemperatureInterpolation = False
options.enableTankWallHeatTransfer = True
options.printFailedInjector = False
options.printTime = False

options.solveTankVolumeContraintWithDAE = True
options.rootFindingType = "brentq"
options.enableDAESolver = True

options.combustion_rtol = 1e-4
options.combustion_atol = 1e-4
options.flight_rtol = 1e-4
options.flight_atol = 1e-4

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
          cpuTime, timeRanges, times, modelsOutput, nfev, njev, nlu, firstSolverTime, nfev2, njev2, nlu2, secondSolverTime = runSystem()

        print(
          "tol={:.2e}, interpolationEnabled={:}, enableDAESolver={:}, rootFindingType={:}, nfev = {:}, njev = {:}, nlu = {:}, firstSolverTime = {:.2f}, nfev2 = {:}, njev2 = {:}, nlu2 = {:}, secondSolverTime = {:.2f}".format(
            tols,
            interpolationEnabled,
            enableDAESolver,
            rootFindingType, nfev, njev, nlu, firstSolverTime, nfev2, njev2, nlu2, secondSolverTime
          )
        )