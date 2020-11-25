import numpy as np
import matplotlib.pyplot as plt

from equations.falloffs import sigmoid

import assumptions

def getDragCoefficient(mach):
  if mach > 10:
    mach = 10

  a = 0.9
  b = -0.6
  c = assumptions.dragPeakAroundMach.get()

  C = assumptions.dragBaseLevel.get()
  A = assumptions.dragPeak.get() - C
  B = A
  
  baseEdge = sigmoid(mach, c, assumptions.dragPeakSmoothingRadius.get())
  return C  + (1.0 - baseEdge) * A*np.exp(a*(mach - c)) + baseEdge * B*np.exp(b*(mach - c))
