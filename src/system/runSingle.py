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
  events = None,
  useDae = False
):
  times = np.arange(fromTime, toTime, outputTimestep)

  # sys is of the form def system(t, y)
  if useDae:
    from scikits.odes import dae
    def addafunc(t, y, ml, mu, p, nrowp):
      return p - np.eye(nrowp)
      
    def residual(t, y, ydot, result):
      try:
        derivative = sys(t, y)
      except Exception:
        return 1

      N = np.size(derivative)
      for i in range(N):
        result[i] = ydot[i] - derivative[0, i]
      
      # print(t, np.linalg.norm(result))
      # print(", ".join(["{:9.2e}".format(result[i]) for i in range(N)]))
      # print(", ".join(["{:9.2e}".format(ydot[i]) for i in range(N)]))
      # print(", ".join(["{:9.2e}".format(derivative[0, i]) for i in range(N)]))
      
      return 0
    
    daeSolver = "lsodi"

    rtol = 1e-3
    atol = 1e-3

    if daeSolver == "ida":
      solver = dae('ida', residual, 
        exclude_algvar_from_error=True,
        adda_func=addafunc,
        algebraic_vars_idx=None,
        atol=atol,
        rtol=rtol,
        max_steps=500000
      )

      yp0 = sys(fromTime, y0)
      solution = solver.solve(times, y0, yp0[0,:])

      if not solution.flag == 1:
        raise Exception("Flag is not 0 (it is " + str(solution.flag) + ")")

      t = np.array(solution.values.t)
      v = np.array(solution.values.y)
      return np.transpose(t), np.transpose(v), t[-1]
    elif daeSolver == "lsodi":
      solver = dae('lsodi', residual, 
        method="adams",
        order=12,
        exclude_algvar_from_error=True,
        adda_func=addafunc,
        algebraic_vars_idx=None,
        atol=atol,
        rtol=rtol,
        max_steps=500000
      )

      yp0 = sys(fromTime, y0)
      solution = solver.solve(times, y0, yp0[0,:])

      flag, t, v, *rest = solution
      return np.transpose(t), np.transpose(v), t[-1]
    else:
      raise Exception("Unknown solver")
    # t, y = values
    #return t, y, t[-1]

  sol = solve_ivp(
    sys, 
    [fromTime, toTime], 
    y0, 
    'LSODA', 
    t_eval=times, 
    dense_output=False, 
    events=events
  )
  
  print(np.shape(sol.t))
  print(np.shape(sol.y))
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

maximumSolveTime = 10
options.printTime = False

models, system, events, initialState, hit_ground_event = makeODE()
t, y, tend = solver(system, initialState, 0, maximumSolveTime, 0.05, events, useDae=True)
print(tend)
recoverModelState(t, y, models)

models2, system2, events2, initialState2, hit_ground_event2 = makeODE(simplified = True, timeHistory=t, previousModels=models)
t2, y2, tend2 = solver(system2, y[:,-1], tend, tend + 240, 0.1, events=(hit_ground_event2), useDae=False)
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