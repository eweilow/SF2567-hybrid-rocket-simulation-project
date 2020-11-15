import numpy as np

from solve.ode import makeODE
from scipy.integrate import solve_ivp

from solve.dependencies import recurseModelDependencies, applyDerivedVariablesToResult
from solve.states import applyModelStates, collectModelStates

models, system, events, initialState = makeODE()

maximumSolveTime = 40

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

with open('/data/simulation.npy', 'wb') as f:
  np.save(f, sol.t)
  
  for key in models:
    np.save(f, models[key]["state"])
    np.save(f, models[key]["derivedResult"])
