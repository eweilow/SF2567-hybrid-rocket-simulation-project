import numpy as np
import time
import options
import os

from utils.cea import NasaCEA
from solve.ode import makeODE
from scipy.integrate import solve_ivp

import assumptions

from solve.dependencies import recurseModelDependencies, applyDerivedVariablesToResult
from solve.states import applyModelStates, collectModelStates

from multiprocessing import Process, Lock

from runSingle import solver, runSystem

def stochastic(fileLock, N, output):
  print('module name:', __name__)
  print('parent process:', os.getppid())
  print('process id:', os.getpid())

  for i in range(N):
    start = time.time()
    # Reset all variables to defaults
    assumptions.Variable.reset()
    assumptions.tankFillingGrade.randomizeInRange(0.85, 0.95)
    assumptions.tankFilledTemperature.randomizeInRange(273, 298)

    assumptions.dragLevelAtPeak.randomize(0.05)
    assumptions.dragLevelAtZero.randomize(0.05)
    assumptions.dragLevelMachAsymptote.randomize(0.05)
    assumptions.dragDropoffConstant.randomize(0.05)
    assumptions.dragPeakMachNumber.randomize(0.05)
    assumptions.dragPeakSmoothingRadius.randomize(0.05)
    
    assumptions.combustionEfficiency.randomize(0.05)
    assumptions.fuelDensitySolid.randomize(0.05)
    assumptions.fuelDensityLiquid.randomize(0.05)
    assumptions.fuelEnthalpyOfFormation.randomize(0.01)
    assumptions.carbonBlackFraction.randomize(0.05)
    assumptions.fuelGrainAConstant.randomize(0.001)
    assumptions.fuelGrainNConstant.randomize(0.001)
    assumptions.fuelPortInitialRadius.randomize(0.001)
    assumptions.fuelPortLength.randomize(0.001)
    assumptions.fuelPortMaximumRadius.randomize(0.001)
    assumptions.initialAtmosphericPressure.randomize(0.01)
    assumptions.injectorHoleDiameter.randomize(0.01)
    assumptions.injectorHoleDischargeCoefficient.randomize(0.05)
    assumptions.launchTowerDirectionAngle.randomize(0.05)
    assumptions.launchSeaLevelAltitude.randomize(0.05)
    assumptions.launchTowerLength.randomize(0.05)
    assumptions.launchTowerDirectionAngle.randomize(0.05)

    assumptions.nozzleEfficiency.randomize(0.05)
    assumptions.nozzleErosionConstant.randomize(0.05)
    assumptions.nozzleErosionStart.randomize(0.05)
    assumptions.nozzleErosionStartRadius.randomize(0.05)

    assumptions.nozzleExhaustRadius.randomize(0.01)
    assumptions.nozzleThroatRadius.randomize(0.01)

    assumptions.preCombustionChamberVolume.randomize(0.1)
    assumptions.postCombustionChamberVolume.randomize(0.1)

    assumptions.tankPassiveVentDischargeCoefficient.randomize(0.05)
    assumptions.tankPassiveVentDiameter.randomize(0.01)

    assumptions.rocketBodyDiameter.randomize(0.01)

    assumptions.rocketOnBoardRecoverySystemMass.randomize(0.05)
    assumptions.rocketOnBoardElectronicsMass.randomize(0.05)
    assumptions.rocketPayloadMass.randomize(0.05)
    assumptions.rocketBodyMass.randomize(0.05)
    assumptions.rocketOxidizerTankMass.randomize(0.05)
    assumptions.rocketFuelCasingMass.randomize(0.05)
    assumptions.kastrullMass.randomize(0.05)
    assumptions.fastenersMass.randomize(0.05)
    assumptions.injectorMass.randomize(0.05)
    assumptions.chamberWallsMass.randomize(0.05)
    assumptions.nozzleMass.randomize(0.05)

    assumptions.combustionEfficiencyStartupTransientTime.randomize(0.1)
    assumptions.injectorStartupTransientTime.randomize(0.1)
    
    options.printTime = False

    cpuTime, timeRanges, times, modelsOutput = runSystem()
    print("Running simulation took {:}".format(cpuTime))
    
    with fileLock:
      with open(output, 'ab') as f:
        np.save(f, timeRanges)
        np.save(f, times)
        for key in modelsOutput:
          np.save(f, modelsOutput[key]["state"])
          np.save(f, modelsOutput[key]["derivedResult"])

        np.save(f, cpuTime)
      # with open('/data/simulation.npy', 'wb') as f:
      #   np.save(f, sol.t)
      #   
      #   for key in models:
      #     np.save(f, models[key]["state"])
      #     np.save(f, models[key]["derivedResult"])



if __name__ == '__main__':
  P = 6
  N = 200

  with open('/data/montecarlo.npy', "wb") as f:
    np.save(f, [P * N])
  
  fileLock = Lock()

  print("If each takes ~16 seconds, time to completion is ~{:.1f} minutes".format(N*16/60))

  processes = []
  for num in range(P):
    p = Process(target=stochastic, args=(fileLock,N,'/data/montecarlo.npy'))
    p.start()
    processes.append(p)

  for p in processes:
    p.join()