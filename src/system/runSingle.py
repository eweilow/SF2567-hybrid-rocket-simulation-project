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

  if useDae:
    print("Running solution from {:.2f} to {:.2f} with DAE solver".format(fromTime, toTime))
  else:
    print("Running solution from {:.2f} to {:.2f} with ODE solver".format(fromTime, toTime))

  options.solvingWithDAE = useDae
  # sys is of the form def system(t, y)
  if useDae:
    from scikits.odes import dae

    yp0 = sys(fromTime, y0)
    yp0 = yp0[0,:]

    def addafunc(t, y, ml, mu, p, nrowp):
      return p - np.eye(nrowp)
      
    N = np.size(y0)
    def residual(t, y, ydot, result):
      try:
        derivative = sys(t, y[0:N])
      except Exception:
        return 1

      gamma = derivative[0, 10]
      volume = derivative[0, 11]
      dgamma_dt = ydot[10]
      dvolume_dt = ydot[11]

      pressure = y[7]

      derivative[0, 7] += pressure / (volume - 1) * dgamma_dt
      derivative[0, 7] += pressure * gamma / volume * dvolume_dt

      tankVolume = derivative[0, 4]

      for i in range(N):
        result[i] = ydot[i] - derivative[0, i]
      

      # DAE part... not pretty but yes, for now, it works
      result[4] = tankVolume
      result[10] = ydot[10]
      result[11] = ydot[11]
      result[N] = y[N] - derivative[0, 10]
      result[N+1] = y[N+1] - derivative[0, 11]

      return 0

    rtol = 1e-4
    atol = 1e-4

    extra_vars = 2
    if extra_vars > 0:
      y0 = np.concatenate((y0, np.zeros(extra_vars)))
      yp0 = np.concatenate((yp0, np.zeros(extra_vars)))

    y0[10] = 0
    yp0[10] = 0
    y0[11] = 0
    yp0[11] = 0
    y0[N] = yp0[10]
    yp0[N] = 0
    y0[N+1] = yp0[11]
    yp0[N+1] = 0
    algebraic_vars_idx = [4, 7, 10, 11, N, N+1]

    def root_fn(t, y, yp, result):
      for i in range(len(events)):
        result[i] = events[i](t, y[0:N])

      # print(t, result)

    solver = dae('ida', residual, 
      compute_initcond="yp0",
      exclude_algvar_from_error=False,
      algebraic_vars_idx=algebraic_vars_idx,
      atol=atol,
      rtol=rtol,
      max_steps=500000,
      rootfn=root_fn,
      nr_rootfns=len(events)
    )


    solution = solver.solve(times, y0, yp0)

    # 0 is done
    # 2 is that it hit root
    if not (solution.flag == 0 or solution.flag == 2):
      raise Exception("Flag is not 0 or 2 (it is " + str(solution.flag) + ")")

    t = np.array(solution.values.t)
    v = np.array(solution.values.y)
    return np.transpose(t), np.transpose(v[:,0:N]), t[-1]

  sol = solve_ivp(
    sys, 
    [fromTime, toTime], 
    y0, 
    'LSODA', 
    t_eval=times, 
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

def runSystem():
  start = time.time()

  maximumSolveTime = 100
  # options.printTime = False

  solveWithDAE = True

  models, system, events, initialState, hit_ground_event = makeODE()
  t, y, tend = solver(system, initialState, 0, maximumSolveTime, 0.01, events, useDae=solveWithDAE)
  recoverModelState(t, y, models)

  models2, system2, events2, initialState2, hit_ground_event2 = makeODE(simplified = True, timeHistory=t, previousModels=models)
  t2, y2, tend2 = solver(system2, y[:,-1], tend, tend + 240, 0.25, events=(hit_ground_event2), useDae=False)
  applyDerivedVariablesToResult(t2, y2, models2)
  recoverModelState(t2, y2, models2)

  end = time.time()

  cpuTime = end - start

  timeRanges = [0, t[-1], t2[-2]]
  times = np.concatenate((t, t2))
  
  modelsOutput = {}
  for key in models:
    modelsOutput[key] = {
      "state": np.concatenate((models[key]["state"], models2[key]["state"]), axis=1),
      "derivedResult": np.concatenate((models[key]["derivedResult"], models2[key]["derivedResult"]), axis=1),
    }

  return cpuTime, timeRanges, times, modelsOutput

if __name__ == '__main__':
  options.printTime = True

  cpuTime, timeRanges, times, modelsOutput = runSystem()
  print("Running simulation took {:}".format(cpuTime))

  with open("/data/simulation.npy", 'wb') as f:
    np.save(f, timeRanges)
    np.save(f, times)
    for key in modelsOutput:
      np.save(f, modelsOutput[key]["state"])
      np.save(f, modelsOutput[key]["derivedResult"])

    np.save(f, cpuTime)
