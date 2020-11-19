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

sol2 = solve_ivp(
  system, 
  [sol.t[-1], sol.t[-1] + 0.1], 
  sol.y[:,-1], 
  'LSODA', 
  t_eval=np.linspace(sol.t[-1], sol.t[-1] + 0.1, 25), 
  dense_output=False, 
  #max_step=1e-1
)
applyDerivedVariablesToResult(sol.t, sol.y, models)

# Recover state for different parts
for key in models:
  models[key]["state"] = np.dot(models[key]["invMatrix"], sol.y)

end = time.time()
print("Running simulation took {:}".format(end - start))

with open("/data/simulation.npy", 'wb') as f:
  np.save(f, sol.t)
  for key in models:
    np.save(f, models[key]["state"])
    np.save(f, models[key]["derivedResult"])

  np.save(f, [end - start])