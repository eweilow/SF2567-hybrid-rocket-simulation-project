import numpy as np
import time
import options
import os
import random

from utils.cea import NasaCEA
from solve.ode import makeODE
from scipy.integrate import solve_ivp

import assumptions

from solve.dependencies import recurseModelDependencies, applyDerivedVariablesToResult
from solve.states import applyModelStates, collectModelStates

from multiprocessing import Process, Lock

from runSingle import solver, runSystem

def stochastic(fileLock, N, output):
  random.seed(int(os.getpid() + time.time()))
  print('module name:', __name__)
  print('parent process:', os.getppid())
  print('process id:', os.getpid())

  mode = None
  for i in range(N):
    if i % 5 == 0:
      mode = None

    mode, baseline, currentValue = assumptions.randomizeOne(mode)
    
    options.printTime = False

    try:
      cpuTime, timeRanges, times, modelsOutput = runSystem()
      print("Running simulation took {:}".format(cpuTime))
      
      with fileLock:
        with open(output, 'ab') as f:
          np.save(f, mode)
          np.save(f, baseline)
          np.save(f, currentValue)
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
    except:
      pass


if __name__ == '__main__':
  P = 4
  N = 5000

  # with open('/data/sensitivity.npy', "wb") as f:
  #   np.save(f, [P * N])
  
  fileLock = Lock()

  print("If each takes ~50 seconds, time to completion is ~{:.1f} minutes".format(N*50/60))

  processes = []
  for num in range(P):
    p = Process(target=stochastic, args=(fileLock,N,'/data/sensitivity.npy'))
    p.start()
    processes.append(p)

  for p in processes:
    p.join()