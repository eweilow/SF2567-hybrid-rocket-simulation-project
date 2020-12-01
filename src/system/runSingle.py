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
  useDae = False,
  rtol = None,
  atol = None
):
  times = np.arange(fromTime, toTime, outputTimestep)

  options.currentlySolvingWithDAE = useDae
  # sys is of the form def system(t, y)
  if useDae:
    from scikits.odes import dae

    yp0 = sys(fromTime, y0)
    yp0 = yp0[0,:]

    def addafunc(t, y, ml, mu, p, nrowp):
      return p - np.eye(nrowp)
    
    nfev = 0
    njev = None
    nlu = None

    N = np.size(y0)
    def residual(t, y, ydot, result):
      nonlocal nfev
      nfev += 1
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

      for i in range(N):
        result[i] = ydot[i] - derivative[0, i]
      

      # DAE part... not pretty but yes, for now, it works
      if options.solveTankVolumeContraintWithDAE:
        tankVolume = derivative[0, 4]
        result[4] = tankVolume
      else:
        result[4] = y[4]

      result[10] = ydot[10]
      result[11] = ydot[11]
      result[N] = y[N] - derivative[0, 10]
      result[N+1] = y[N+1] - derivative[0, 11]

      return 0

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
        if hasattr(events[i], "disabled") and events[i].disabled:
          result[i] = 1
          continue

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


    allTimes = np.zeros((0))
    allY = np.zeros((N,0))
    currentTime = 0
    
    currentY0 = y0
    currentYP0 = yp0


    while True:
      computeTimes = np.concatenate(([currentTime], times[times > currentTime]))

      if len(computeTimes) == 0:
        break

      print("Running solution from {:.4f} to {:.4f} with DAE (IDA) solver".format(np.min(computeTimes), np.max(computeTimes)))
      solution = solver.solve(computeTimes, currentY0, currentYP0)

      # 0 is done
      # 2 is that it hit root
      if not (solution.flag == 0 or solution.flag == 2):
        raise Exception("Flag is not 0 or 2 (it is " + str(solution.flag) + ")")

      t = np.array(solution.values.t)
      v = np.array(solution.values.y)
      allTimes = np.concatenate((allTimes, np.transpose(t)))
      allY = np.concatenate((allY, np.transpose(v[:,0:N])), axis=1)

      if solution.flag == 2:
        root = solution.roots
        result = np.zeros(len(events))
        root_fn(root.t[0], root.y[0], root.ydot[0], result)
        result = np.abs(result)
        argmin = np.argmin(result)
        event = events[argmin]
        event.disabled = True

        currentTime = root.t[0]
        currentY0 = root.y[0]
        currentYP0 = root.ydot[0]
        print("Root '{:}' hit at {:.4f}".format(event.__name__, currentTime))
        if hasattr(event, "terminal"):
          if event.terminal:
            break # stop!
          

    return allTimes, allY, allTimes[-1], nfev, njev, nlu
    
  print("Running solution from {:.2f} to {:.2f} with ODE (LSODA) solver".format(fromTime, toTime))
  sol = solve_ivp(
    sys, 
    [fromTime, toTime], 
    y0, 
    'LSODA', 
    t_eval=times, 
    dense_output=False, 
    events=events,
    atol=atol,
    rtol=rtol
  )

  for i in range(len(sol.t_events)):
    if(len(sol.t_events[i]) == 0):
      continue
    sol.t = np.concatenate((sol.t, sol.t_events[i]), axis=0)
    sol.y = np.concatenate((sol.y, np.transpose(sol.y_events[i])), axis=1)
  #   sol.t = np.concatenate(sol.t, sol.t_events[i])
  #   sol.y = np.concatenate(sol.y, sol.y_events[i])
  # 
  # sort_indices = np.argsort(sol.t, axis=0)
  # sol.t = sol.t[sort_indices]
  # sol.y = sol.y[sort_indices]

  indices = np.argsort(sol.t)
  sol.t = sol.t[indices]
  sol.y = sol.y[:,indices]
  return sol.t, sol.y, sol.t[-1], sol.nfev, sol.njev, sol.nlu

def recoverModelState(t, y, models, canApplySimplified = False):
  applyDerivedVariablesToResult(t, y, models)
  # Recover state for different parts
  for key in models:
    models[key]["state"] = np.dot(models[key]["invMatrix"], y)

    if canApplySimplified and not models[key]["simplifiedInit"] is None:
      mask, args = models[key]["simplifiedInit"]
      for i in range(len(mask)):
        if mask[i] is None:
          models[key]["state"][i,:] = None

def runSystem():
  start = time.time()

  maximumSolveTime = 25
  # options.printTime = False

  models, system, events, initialState, hit_ground_event, hit_apoapsis = makeODE()
  solve1 = time.time()
  t, y, tend, nfev, njev, nlu = solver(system, initialState, 0, maximumSolveTime, 0.01, events, useDae=options.enableDAESolver, rtol=options.combustion_rtol, atol=options.combustion_atol)
  recoverModelState(t, y, models, canApplySimplified = False)

  solve2 = time.time()
  models2, system2, events2, initialState2, hit_ground_event2, hit_apoapsis2 = makeODE(simplified = True, timeHistory=t, previousModels=models)
  t2, y2, tend2, nfev2, njev2, nlu2 = solver(system2, y[:,-1], tend, tend + 240, 0.25, events=(hit_ground_event2, hit_apoapsis2), useDae=False, rtol=options.flight_rtol, atol=options.flight_atol)
  recoverModelState(t2, y2, models2, canApplySimplified = True)
  solve3 = time.time()

  firstSolverTime = solve2 - solve1
  secondSolverTime = solve3 - solve2

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

  return cpuTime, timeRanges, times, modelsOutput, nfev, njev, nlu, firstSolverTime, nfev2, njev2, nlu2, secondSolverTime

if __name__ == '__main__':
  # options.printTime = True

  cpuTime, timeRanges, times, modelsOutput, *rest = runSystem()
  print("Running simulation took {:}".format(cpuTime))

  with open("/data/simulation.npy", 'wb') as f:
    np.save(f, timeRanges)
    np.save(f, times)
    for key in modelsOutput:
      np.save(f, modelsOutput[key]["state"])
      np.save(f, modelsOutput[key]["derivedResult"])

    np.save(f, cpuTime)
