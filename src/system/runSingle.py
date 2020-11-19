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


start = time.time()

models, system, events, initialState, hit_ground_event = makeODE()
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

models2, system2, events2, initialState2, hit_ground_event2 = makeODE(simplified = True, timeHistory=sol.t, previousModels=models)

sol2 = solve_ivp(
  system2, 
  [sol.t[-1], sol.t[-1] + 120], 
  sol.y[:,-1], 
  'LSODA', 
  t_eval=np.linspace(sol.t[-1], sol.t[-1] + 120, 100), 
  dense_output=False, 
  events=(hit_ground_event2)
  #max_step=1e-1
)
applyDerivedVariablesToResult(sol2.t, sol2.y, models2)

# Recover state for different parts
for key in models2:
  models2[key]["state"] = np.dot(models2[key]["invMatrix"], sol2.y)

  if not models[key]["simplifiedInit"] is None:
    mask, args = models[key]["simplifiedInit"]
    for i in range(len(mask)):
      if mask[i] is None:
        models2[key]["state"][i,:] = None


end = time.time()
print("Running simulation took {:}".format(end - start))

with open("/data/simulation.npy", 'wb') as f:
  np.save(f, [0, sol.t[-1], sol2.t[-2]])
  np.save(f, np.concatenate((sol.t, sol2.t)))
  for key in models:
    np.save(f, np.concatenate((models[key]["state"], models2[key]["state"]), axis=1))
    np.save(f, np.concatenate((models[key]["derivedResult"], models2[key]["derivedResult"]), axis=1))

  np.save(f, [end - start])