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

from multiprocessing import Process, Lock

from runSingle import solver, runSystem

def stochastic(fileLock, N, output):
  print('module name:', __name__)
  print('parent process:', os.getppid())
  print('process id:', os.getpid())

  for i in range(N):
    # Reset all variables to defaults
    assumptions.Variable.reset()

    knownLengthTolerance = 0.01 # +- 1%
    knownVolumeTolerance = 0.05 # +- 5%
    knownMassTolerance = 0.025 # +- 2.5%
    knownConstantTolerance = 0.025
    efficiencyTolerance = 0.025  # +- 2.5%

    assumptions.tankOutletLiquidGasTransferRadius.randomize(knownConstantTolerance)
    assumptions.tankTopVentLiquidGasTransferRadius.randomize(knownConstantTolerance)
    assumptions.initialAtmosphericPressure.randomize(knownConstantTolerance)
    assumptions.initialAtmosphericTemperature.randomize(knownConstantTolerance)
    assumptions.tankFillingGrade.randomizeInRange(0.85, 0.99)
    assumptions.tankFilledTemperature.randomizeInRange(268, 298)
    assumptions.tankVolume.randomize(0.025)
    assumptions.tankThickness.randomize(knownLengthTolerance)
    assumptions.tankWallDensity.randomize(knownConstantTolerance)
    assumptions.tankWallThermalConductivity.randomize(knownConstantTolerance)
    assumptions.tankWallSpecificHeatCapacity.randomize(knownConstantTolerance)
    assumptions.tankPassiveVentDischargeCoefficient.randomize(knownConstantTolerance)
    assumptions.tankPassiveVentDiameter.randomize(knownLengthTolerance)
    assumptions.fuelDensitySolid.randomize(knownConstantTolerance)
    assumptions.fuelDensityLiquid.randomize(knownConstantTolerance)
    assumptions.fuelEnthalpyOfFormation.randomize(knownConstantTolerance)
    assumptions.fuelPortLength.randomize(knownLengthTolerance)
    assumptions.fuelPortInitialRadius.randomize(knownLengthTolerance)
    assumptions.carbonBlackFraction.randomize(knownConstantTolerance)
    assumptions.preCombustionChamberVolume.randomize(knownVolumeTolerance)
    assumptions.postCombustionChamberVolume.randomize(knownVolumeTolerance)
    assumptions.combustionEfficiency.randomize(efficiencyTolerance)
    assumptions.fuelGrainAConstant.randomizeInRange(0.132e-3 - 0.015e-3, 0.155e-3 + 0.015e-3)
    assumptions.fuelGrainNConstant.randomizeInRange(0.5-0.015, 0.555+0.015)
    assumptions.injectorHoleCount.randomize(knownConstantTolerance)
    assumptions.injectorHoleDischargeCoefficient.randomize(knownConstantTolerance)
    assumptions.injectorHoleDiameter.randomize(knownLengthTolerance)
    assumptions.nozzleEfficiency.randomize(efficiencyTolerance)
    assumptions.nozzleExhaustRadius.randomize(knownLengthTolerance)
    assumptions.nozzleThroatRadius.randomize(knownLengthTolerance)
    
    assumptions.launchLatitudeDegrees.randomize(knownConstantTolerance)
    assumptions.launchLongitudeDegrees.randomize(knownConstantTolerance)
    assumptions.launchSeaLevelAltitude.randomize(knownConstantTolerance)
    assumptions.launchTowerLength.randomizeInRange(9, 11)
    assumptions.launchTowerVerticalAngle.randomizeInRange(78, 82)
    assumptions.launchTowerDirectionAngle.randomizeInRange(0, 360)
    assumptions.rocketBodyDiameter.randomize(knownLengthTolerance)
    assumptions.rocketOnBoardRecoverySystemMass.randomize(knownMassTolerance)
    assumptions.rocketOnBoardElectronicsMass.randomize(knownMassTolerance)
    assumptions.rocketPayloadMass.randomize(knownMassTolerance)
    assumptions.rocketBodyMass.randomize(knownMassTolerance)
    assumptions.rocketOxidizerTankMass.randomize(knownMassTolerance)
    assumptions.rocketFuelCasingMass.randomize(knownMassTolerance)
    assumptions.kastrullMass.randomize(knownMassTolerance)
    assumptions.fastenersMass.randomize(knownMassTolerance)
    assumptions.injectorMass.randomize(knownMassTolerance)
    assumptions.chamberWallsMass.randomize(knownMassTolerance)
    assumptions.nozzleMass.randomize(knownMassTolerance)
    assumptions.rocketEngineMass.randomize(knownMassTolerance)

    assumptions.dragLevelAtZero.randomize(knownConstantTolerance)
    assumptions.dragLevelAtPeak.randomize(knownConstantTolerance)
    assumptions.dragPeakMachNumber.randomize(knownConstantTolerance)
    assumptions.dragLevelMachAsymptote.randomize(knownConstantTolerance)
    assumptions.dragDropoffConstant.randomize(knownConstantTolerance)
    assumptions.dragPeakSmoothingRadius.randomize(knownConstantTolerance)
    assumptions.combustionEfficiencyStartupTransientTime.randomize(knownConstantTolerance)
    assumptions.combustionEfficiencyStartupDelay.randomize(knownConstantTolerance)
    assumptions.injectorStartupTransientTime.randomize(knownConstantTolerance)
    assumptions.injectorStartupDelay.randomize(knownConstantTolerance)
    assumptions.maximumRegressionRateAt.randomize(knownConstantTolerance)
    
    options.printTime = False
    try:
      cpuTime, timeRanges, times, modelsOutput, *rest = runSystem()
      print("Running simulation took {:}".format(cpuTime))
      
      with fileLock:
        with open(output, 'ab') as f:
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
  P = 6
  N = 5000

  fileLock = Lock()

  print("If each takes ~50 seconds, time to completion is ~{:.1f} minutes".format(N*40/60))

  processes = []
  for num in range(P):
    p = Process(target=stochastic, args=(fileLock,N,'/data/montecarlo.npy'))
    p.start()
    processes.append(p)

  for p in processes:
    p.join()