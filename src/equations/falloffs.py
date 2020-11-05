import matplotlib.pyplot as plt
import numpy as np

def sigmoid(x, center, radius):
  return 1 / (1 + np.exp(-((x - center)/(radius/np.pi/2))))

def outletPhaseFalloff(liquidLevel):
  bottomFalloff = 0.01
  return 1.0 - sigmoid(liquidLevel, bottomFalloff, bottomFalloff)

def inletPhaseFalloff(liquidLevel):
  topFalloff = 0.02
  return 1.0 - sigmoid(liquidLevel, 1 - topFalloff, topFalloff)

def combustionEfficiencyTransient(time):
  center = 0.25

  minimumLevel = 1e-2
  return minimumLevel + (1 - minimumLevel) * sigmoid(time, center, center)

def injectorTransientFalloff(time):
  center = 0.1

  minimumLevel = 1e-2
  return minimumLevel + (1 - minimumLevel) * sigmoid(time, center, center)



if __name__ == '__main__':
  x = np.linspace(0, 0.05, 100)
  y = []
  for i in range(np.size(x)):
    y.append(outletPhaseFalloff(x[i]))

  plt.subplot(2,2,1)
  plt.plot(x, y)


  x = np.linspace(0.95, 1, 100)
  y = []
  for i in range(np.size(x)):
    y.append(inletPhaseFalloff(x[i]))

  plt.subplot(2,2,2)
  plt.plot(x, y)


  x = np.linspace(0, 3, 100)
  y = []
  for i in range(np.size(x)):
    y.append(combustionEfficiencyTransient(x[i]))

  plt.subplot(2,2,3)
  plt.plot(x, y)

  plt.show()