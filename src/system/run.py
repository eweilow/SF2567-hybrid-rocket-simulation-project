import numpy as np
import time
import options

from utils.cea import NasaCEA
from solve.ode import makeODE
from scipy.integrate import solve_ivp

import assumptions

from solve.dependencies import recurseModelDependencies, applyDerivedVariablesToResult
from solve.states import applyModelStates, collectModelStates

N = 200
with open('/data/montecarlo.npy', 'wb') as f:
  np.save(f, [N])
  times = []
  for i in range(N):
    start = time.time()
    # Reset all variables to defaults
    assumptions.Variable.reset()
    assumptions.tankFillingGrade.randomizeInRange(0.75, 0.95)
    assumptions.tankFilledTemperature.randomizeInRange(273, 293)
    
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

    times.append(end - start)
    np.save(f, sol.t)
    for key in models:
      np.save(f, models[key]["state"])
      np.save(f, models[key]["derivedResult"])

  np.save(f, times)
    # with open('/data/simulation.npy', 'wb') as f:
    #   np.save(f, sol.t)
    #   
    #   for key in models:
    #     np.save(f, models[key]["state"])
    #     np.save(f, models[key]["derivedResult"])
