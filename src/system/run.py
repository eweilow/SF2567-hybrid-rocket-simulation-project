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

def stochastic(fileLock, N, output):
  print('module name:', __name__)
  print('parent process:', os.getppid())
  print('process id:', os.getpid())

  for i in range(N):
    start = time.time()
    # Reset all variables to defaults
    assumptions.Variable.reset()
    assumptions.tankFillingGrade.randomizeInRange(0.75, 0.95)
    assumptions.tankFilledTemperature.randomizeInRange(273, 293)

    assumptions.dragPeak.randomize(0.05)
    assumptions.dragBaseLevel.randomize(0.05)
    assumptions.dragPeakAroundMach.randomize(0.05)
    assumptions.dragPeakSmoothingRadius.randomize(0.05)
    
    assumptions.combustionEfficiency.randomize(0.05)
    assumptions.fuelDensity.randomize(0.05)
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

    assumptions.tankVolume.randomize(0.01)

    assumptions.tankPassiveVentDischargeCoefficient.randomize(0.05)
    assumptions.tankPassiveVentDiameter.randomize(0.01)

    assumptions.rocketBodyDiameter.randomize(0.01)

    assumptions.rocketOnBoardRecoverySystemMass.randomize(0.05)
    assumptions.rocketOnBoardElectronicsMass.randomize(0.05)
    assumptions.rocketPayloadMass.randomize(0.05)
    assumptions.rocketBodyMass.randomize(0.05)
    assumptions.rocketOxidizerTankMass.randomize(0.05)
    assumptions.rocketFuelCasingMass.randomize(0.05)
    assumptions.rocketEngineMass.randomize(0.05)
    
    models, system, events, initialState = makeODE()
    maximumSolveTime = 40
    options.printTime = False

    sol = solve_ivp(
      system, 
      [0, maximumSolveTime], 
      initialState, 
      'LSODA', 
      t_eval=np.linspace(0, maximumSolveTime, 500), 
      dense_output=False, 
      events=events,
      #max_step=1e-1
    )

    applyDerivedVariablesToResult(sol.t, sol.y, models)

    # Recover state for different parts
    for key in models:
      models[key]["state"] = np.dot(models[key]["invMatrix"], sol.y)

    end = time.time()
    print("Running simulation {:} took {:}".format(i, end - start))
    
    with fileLock:
      with open(output, 'ab') as f:
        np.save(f, sol.t)
        for key in models:
          np.save(f, models[key]["state"])
          np.save(f, models[key]["derivedResult"])

        np.save(f, [end - start])
      # with open('/data/simulation.npy', 'wb') as f:
      #   np.save(f, sol.t)
      #   
      #   for key in models:
      #     np.save(f, models[key]["state"])
      #     np.save(f, models[key]["derivedResult"])



if __name__ == '__main__':
  P = 4
  N = 250

  with open('/data/montecarlo.npy', "wb") as f:
    np.save(f, [P * N])
  
  fileLock = Lock()

  print("If each takes ~11 seconds, time to completion is ~{:.1f} minutes".format(N*11/60))

  processes = []
  for num in range(P):
    p = Process(target=stochastic, args=(fileLock,N,'/data/montecarlo.npy'))
    p.start()
    processes.append(p)

  for p in processes:
    p.join()