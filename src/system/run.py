import numpy as np

from solve.ode import makeODE
from scipy.integrate import solve_ivp

from solve.dependencies import recurseModelDependencies
from solve.states import applyModelStates, collectModelStates

models, system, events, initialState = makeODE()

maximumSolveTime = 20

sol = solve_ivp(
  system, 
  [0, maximumSolveTime], 
  initialState, 
  'LSODA', 
  t_eval=np.linspace(0, maximumSolveTime, 5000), 
  dense_output=False, 
  events=events,
  max_step=1e-1
)

sol.t = sol.t[1:]
sol.y = sol.y[:,1:]
states = sol.y

for key in models:
  models[key]["derivedResult"] = np.zeros((len(models[key]["derived"]), len(sol.t)))

# Recover derived variables
for i in range(len(sol.t)):
  y = sol.y[:,i]
  print(sol.t[i])
  applyModelStates(models, y)

  # Set up derived variables
  visited = []
  for key in models:
    visited = recurseModelDependencies(sol.t[i], models[key], 0, models, visited)

  for key in models:
    models[key]["derivedResult"][:, i] = models[key]["derived"]


applyModelStates(models, sol.y)
# Recover state for different parts
for key in models:
  models[key]["state"] = np.dot(models[key]["invMatrix"], sol.y)


with open('/data/simulation.npy', 'wb') as f:
  np.save(f, sol.t)
  
  np.save(f, models["tank"]["state"])
  np.save(f, models["tank"]["derivedResult"])
  np.save(f, models["injector"]["state"])
  np.save(f, models["injector"]["derivedResult"])
  np.save(f, models["combustion"]["state"])
  np.save(f, models["combustion"]["derivedResult"])
  np.save(f, models["nozzle"]["state"])
  np.save(f, models["nozzle"]["derivedResult"])

