from models.tank import TankModel
from models.injector import InjectorModel
from models.combustion import CombustionModel
from models.nozzle import NozzleModel
import matplotlib.pyplot as plt
import numpy as np
import utils.constants as constants
from scipy.integrate import solve_ivp

tank = TankModel()
injector = InjectorModel()
combustion = CombustionModel()
nozzle = NozzleModel()

modelsList = [tank, injector, combustion, nozzle]

visited = []
def recurse(t, obj, currentDepth):
  if obj in visited:
    # print("passing", obj["name"], "at depth", currentDepth)
    return

  # print("computing", obj["name"], "at depth", currentDepth)
  
  visited.append(obj)
  dependsOn = obj["model"].derivedVariablesDependsOn(models)

  for dependency in dependsOn:
    recurse(t, dependency, currentDepth + 1)
  
  obj["derived"] = obj["model"].computeDerivedVariables(t, obj["state"], models)

models = {
  "tank": {
    "name": "tank",
    "model": tank,
    "state": np.array(tank.initializeState()),
    "derived": [],
  },
  "injector": {
    "name": "injector",
    "model": injector,
    "state": np.array(injector.initializeState()),
    "derived": [],
  },
  "combustion": {
    "name": "combustion",
    "model": combustion,
    "state": np.array(combustion.initializeState()),
    "derived": [],
  },
  "nozzle": {
    "name": "nozzle",
    "model": nozzle,
    "state": np.array(nozzle.initializeState()),
    "derived": [],
  },
}

tankStateLength = np.size(models["tank"]["state"])
injectorStateLength = np.size(models["injector"]["state"])
combustionStateLength = np.size(models["combustion"]["state"])
nozzleStateLength = np.size(models["nozzle"]["state"])
fullSystemLength = tankStateLength + injectorStateLength + combustionStateLength + nozzleStateLength

nextIndex = 0
tankStateMatrix = np.zeros((fullSystemLength, tankStateLength))
for i in range(tankStateLength):
  tankStateMatrix[nextIndex][i] = 1
  nextIndex = nextIndex + 1
tankStateMatrixInverse = np.transpose(tankStateMatrix)

injectorStateMatrix = np.zeros((fullSystemLength, injectorStateLength))
for i in range(injectorStateLength):
  injectorStateMatrix[nextIndex][i] = 1
  nextIndex = nextIndex + 1
injectorStateMatrixInverse = np.transpose(injectorStateMatrix)

combustionStateMatrix = np.zeros((fullSystemLength, combustionStateLength))
for i in range(combustionStateLength):
  combustionStateMatrix[nextIndex][i] = 1
  nextIndex = nextIndex + 1
combustionStateMatrixInverse = np.transpose(combustionStateMatrix)

nozzleStateMatrix = np.zeros((fullSystemLength, nozzleStateLength))
for i in range(nozzleStateLength):
  nozzleStateMatrix[nextIndex][i] = 1
  nextIndex = nextIndex + 1
nozzleStateMatrixInverse = np.transpose(nozzleStateMatrix)

def f(t, y):
  print(t)
  models["tank"]["state"] = np.dot(tankStateMatrixInverse, y)
  models["injector"]["state"] = np.dot(injectorStateMatrixInverse, y)
  models["combustion"]["state"] = np.dot(combustionStateMatrixInverse, y)
  models["nozzle"]["state"] = np.dot(nozzleStateMatrixInverse, y)

  # Set up derived variables
  recurse(t, models["tank"], 0)
  recurse(t, models["injector"], 0)
  recurse(t, models["combustion"], 0)
  recurse(t, models["nozzle"], 0)
  visited.clear()

  tankDerivatives = np.array(models["tank"]["model"].computeDerivatives(t, models["tank"]["state"], models["tank"]["derived"], models))
  injectorDerivatives = np.array(models["injector"]["model"].computeDerivatives(t, models["injector"]["state"], models["injector"]["derived"], models))
  combustionDerivatives = np.array(models["combustion"]["model"].computeDerivatives(t, models["combustion"]["state"], models["combustion"]["derived"], models))
  nozzleDerivatives = np.array(models["nozzle"]["model"].computeDerivatives(t, models["nozzle"]["state"], models["nozzle"]["derived"], models))

  # Merge compartmentalized derivatives into a "global" derivative
  return np.dot(tankStateMatrix, tankDerivatives) + np.dot(injectorStateMatrix, injectorDerivatives) + np.dot(combustionStateMatrix, combustionDerivatives) + np.dot(nozzleStateMatrix, nozzleDerivatives)

def no_oxidizer_mass(t, y): 
  models["tank"]["state"] = np.dot(tankStateMatrixInverse, y)
  return models["tank"]["state"][TankModel.states_oxidizerMass] - 0.25
no_oxidizer_mass.terminal = True

def no_oxidizer_energy(t, y): 
  models["tank"]["state"] = np.dot(tankStateMatrixInverse, y)
  return models["tank"]["state"][TankModel.states_totalEnergy]
no_oxidizer_energy.terminal = True

def no_chamber_pressure(t, y): 
  models["combustion"]["state"] = np.dot(combustionStateMatrixInverse, y)
  return models["combustion"]["state"][CombustionModel.states_pressure] - 2 * constants.Pressure.bar
no_chamber_pressure.terminal = True
no_chamber_pressure.direction = -1


initialState = np.dot(tankStateMatrix, models["tank"]["state"]) + np.dot(injectorStateMatrix, models["injector"]["state"]) + np.dot(combustionStateMatrix, models["combustion"]["state"]) + np.dot(nozzleStateMatrix, models["nozzle"]["state"])

T = 20

sol = solve_ivp(
  f, 
  [0, T], 
  initialState, 
  'LSODA', 
  t_eval=np.linspace(0, T, 250), 
  dense_output=False, 
  events=(no_oxidizer_mass, no_oxidizer_energy, no_chamber_pressure),
  rtol=1e-2,
  atol=1e-2
)
sol.t = sol.t[1:]
sol.y = sol.y[:,1:]
states = sol.y



tankDerivedVariables = np.zeros((len(models["tank"]["derived"]), len(sol.t)))
injectorDerivedVariables = np.zeros((len(models["injector"]["derived"]), len(sol.t)))
combustionDerivedVariables = np.zeros((len(models["combustion"]["derived"]), len(sol.t)))
nozzleDerivedVariables = np.zeros((len(models["nozzle"]["derived"]), len(sol.t)))

# Recover derived variables
for i in range(len(sol.t)):
  y = sol.y[:,i]
  print(sol.t[i])
  models["tank"]["state"] = np.dot(tankStateMatrixInverse, y)
  models["injector"]["state"] = np.dot(injectorStateMatrixInverse, y)
  models["combustion"]["state"] = np.dot(combustionStateMatrixInverse, y)
  models["nozzle"]["state"] = np.dot(nozzleStateMatrixInverse, y)

  # Set up derived variables
  recurse(sol.t[i], models["tank"], 0)
  recurse(sol.t[i], models["injector"], 0)
  recurse(sol.t[i], models["combustion"], 0)
  recurse(sol.t[i], models["nozzle"], 0)
  visited.clear()

  tankDerivedVariables[:, i] = models["tank"]["derived"]
  injectorDerivedVariables[:, i] = models["injector"]["derived"]
  combustionDerivedVariables[:, i] = models["combustion"]["derived"]
  nozzleDerivedVariables[:, i] = models["nozzle"]["derived"]

# Recover state for different parts
models["tank"]["state"] = np.dot(tankStateMatrixInverse, sol.y)
models["injector"]["state"] = np.dot(injectorStateMatrixInverse, sol.y)
models["combustion"]["state"] = np.dot(combustionStateMatrixInverse, sol.y)
models["nozzle"]["state"] = np.dot(nozzleStateMatrixInverse, sol.y)

models["tank"]["derived"] = tankDerivedVariables
models["injector"]["derived"] = injectorDerivedVariables
models["combustion"]["derived"] = combustionDerivedVariables
models["nozzle"]["derived"] = nozzleDerivedVariables

index = 1
width = 5
height = 4



plt.plot(models["combustion"]["derived"][CombustionModel.derived_ofRatio], models["nozzle"]["derived"][NozzleModel.derived_specificImpulse])
plt.xlabel("Oxidizer-Fuel Ratio")
plt.ylabel("Specific Impulse [s]")
plt.grid()
plt.show()

def nextSubplot():
  global index
  plt.subplot(height, width, index)
  plt.grid()
  index = (index) % (width * height) + 1

nextSubplot()
plt.plot(sol.t, models["tank"]["state"][TankModel.states_oxidizerMass], '--')
plt.plot(sol.t, models["tank"]["derived"][TankModel.derived_liquidMass], '-')
plt.plot(sol.t, models["tank"]["derived"][TankModel.derived_gasMass], '-')
plt.legend(("Total", "Liquid", "Gas"))
plt.xlabel("Time (s)")
plt.ylabel("Oxidizer mass (kg)")
plt.title("Oxidizer")

# nextSubplot()
# plt.plot(sol.t, models["combustion"]["state"][CombustionModel.states_burntFuel], '-')
# plt.xlabel("Time (s)")
# plt.ylabel("Burnt mass (kg)")
# plt.title("Fuel")

nextSubplot()
plt.plot(sol.t, models["injector"]["derived"][InjectorModel.derived_massFlow], '-')
plt.plot(sol.t, models["combustion"]["derived"][CombustionModel.derived_fuelFlow], '-')
plt.plot(sol.t, models["nozzle"]["derived"][NozzleModel.derived_massFlow], '--')
plt.legend(("Oxidizer", "Fuel", "Nozzle"))
plt.xlabel("Time (s)")
plt.ylabel("Mass flow (kg/s)")
plt.title("Propellant flow")

nextSubplot()
plt.plot(sol.t, models["tank"]["derived"][TankModel.derived_massFlowTop]*1e3, '-')
plt.xlabel("Time (s)")
plt.ylabel("Mass flow (g/s)")
plt.title("Passive vent")

nextSubplot()
plt.plot(sol.t, models["tank"]["derived"][TankModel.derived_pressure] / constants.Pressure.bar, '-')
plt.plot(sol.t, models["combustion"]["state"][CombustionModel.states_pressure] / constants.Pressure.bar, '-')
plt.plot(sol.t, models["combustion"]["derived"][CombustionModel.derived_exhaustPressure] / constants.Pressure.bar, '-')
plt.legend(("Tank", "Chamber", "Exhaust"))
plt.xlabel("Time (s)")
plt.ylabel("Pressure (bar)")
plt.title("Pressures")

nextSubplot()
plt.plot(sol.t, models["tank"]["derived"][TankModel.derived_outletDensity], '-')
plt.plot(sol.t, models["tank"]["derived"][TankModel.derived_outletTopDensity], '-')
plt.legend(("Main outlet", "Passive vent"))
plt.xlabel("Time (s)")
plt.ylabel("Density (kg/m^3)")
plt.title("Tank outlets")

nextSubplot()
plt.plot(sol.t, models["tank"]["derived"][TankModel.derived_temperature], '-')
plt.xlabel("Time (s)")
plt.ylabel("Temperature (K)")
plt.title("Tank temperature")

nextSubplot()
plt.plot(sol.t, models["tank"]["derived"][TankModel.derived_liquidLevel] * 100, '-')
plt.xlabel("Time (s)")
plt.ylabel("Liquid level (%)")
plt.title("Tank liquid level")

nextSubplot()
plt.plot(sol.t, models["tank"]["derived"][TankModel.derived_vaporQuality] * 100, '-')
plt.xlabel("Time (s)")
plt.ylabel("Vapor quality (%)")
plt.title("Tank vapor quality")

nextSubplot()
plt.plot(sol.t, models["tank"]["derived"][TankModel.derived_liquidDensity], '-')
plt.plot(sol.t, models["tank"]["derived"][TankModel.derived_gasDensity], '-')
plt.legend(('Liquid', 'Gas'))
plt.xlabel("Time (s)")
plt.ylabel("Density (kg/m^3)")
plt.title("Oxidizer density")


nextSubplot()
plt.plot(sol.t, models["combustion"]["state"][CombustionModel.states_portRadius] / constants.Lengths.mm, '-')
plt.xlabel("Time (s)")
plt.ylabel("Radius (mm)")
plt.title("Port radius")


nextSubplot()
plt.plot(sol.t, models["combustion"]["derived"][CombustionModel.derived_volume] / constants.Volume.liter, '-')
plt.xlabel("Time (s)")
plt.ylabel("Volume (liter)")
plt.title("Chamber volume")


nextSubplot()
plt.plot(sol.t, models["combustion"]["derived"][CombustionModel.derived_temperature], '-')
plt.xlabel("Time (s)")
plt.ylabel("Temperature (K)")
plt.title("Adiabatic chamber flame temperature")

nextSubplot()
plt.plot(sol.t, models["combustion"]["derived"][CombustionModel.derived_gamma], '-')
plt.xlabel("Time (s)")
plt.ylabel("Gamma")
plt.title("Combustion gamma")

nextSubplot()
plt.plot(sol.t, models["combustion"]["derived"][CombustionModel.derived_CpT], '-')
plt.xlabel("Time (s)")
plt.ylabel("CpT")
plt.title("Combustion CpT")

nextSubplot()
plt.plot(sol.t, models["combustion"]["derived"][CombustionModel.derived_rDot] / constants.Lengths.mm, '-')
plt.xlabel("Time (s)")
plt.ylabel("Rate (mm/s)")
plt.title("Regression rate")

nextSubplot()
plt.plot(sol.t, models["combustion"]["derived"][CombustionModel.derived_ofRatio], '-')
plt.xlabel("Time (s)")
plt.ylabel("OF")
plt.title("Oxidizer-fuel ratio")

nextSubplot()
plt.plot(sol.t, models["combustion"]["derived"][CombustionModel.derived_cStar], '-')
plt.xlabel("Time (s)")
plt.ylabel("C*")
plt.title("Characteristic velocity")

nextSubplot()
plt.plot(sol.t, models["combustion"]["derived"][CombustionModel.derived_thrustCoefficient], '-')
plt.xlabel("Time (s)")
plt.ylabel("Coefficient")
plt.title("Thrust coefficient")

nextSubplot()
plt.plot(sol.t, models["nozzle"]["derived"][NozzleModel.derived_thrust], '-')
plt.xlabel("Time (s)")
plt.ylabel("Thrust (N)")
plt.title("Thrust")

nextSubplot()
plt.plot(sol.t, models["nozzle"]["derived"][NozzleModel.derived_specificImpulse], '-')
plt.xlabel("Time (s)")
plt.ylabel("Specific impulse (s)")
plt.title("Specific impulse")

plt.tight_layout()
plt.subplots_adjust(
  left=0.05,
  bottom=0.05,
  top=1-0.05,
  right=1-0.05,
  wspace=0.25,
  hspace=0.35
)
plt.show()