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

def solver(
  sys,
  y0,
  fromTime,
  toTime,
  outputTimestep,
  events = None
):
  sol = solve_ivp(
    sys, 
    [fromTime, toTime], 
    y0, 
    'LSODA', 
    t_eval=np.arange(fromTime, toTime, outputTimestep), 
    dense_output=False, 
    events=events
  )

  return sol.t, sol.y, sol.t[-1]

def recoverModelState(t, y, models):
  applyDerivedVariablesToResult(t, y, models)
  # Recover state for different parts
  for key in models:
    models[key]["state"] = np.dot(models[key]["invMatrix"], y)

    if not models[key]["simplifiedInit"] is None:
      mask, args = models[key]["simplifiedInit"]
      for i in range(len(mask)):
        if mask[i] is None:
          models[key]["state"][i,:] = None

start = time.time()

maximumSolveTime = 40
options.printTime = False

models, system, events, initialState, hit_ground_event = makeODE()
t, y, tend = solver(system, initialState, 0, maximumSolveTime, 0.05, events)
recoverModelState(t, y, models)

models2, system2, events2, initialState2, hit_ground_event2 = makeODE(simplified = True, timeHistory=t, previousModels=models)
t2, y2, tend2 = solver(system2, y[:,-1], tend, tend + 240, 0.1, events=(hit_ground_event2))
applyDerivedVariablesToResult(t2, y2, models2)
recoverModelState(t2, y2, models2)

end = time.time()
print("Running simulation took {:}".format(end - start))

with open("/data/simulation.npy", 'wb') as f:
  np.save(f, [0, t[-1], t2[-2]])
  np.save(f, np.concatenate((t, t2)))
  for key in models:
    np.save(f, np.concatenate((models[key]["state"], models2[key]["state"]), axis=1))
    np.save(f, np.concatenate((models[key]["derivedResult"], models2[key]["derivedResult"]), axis=1))

  np.save(f, [end - start])