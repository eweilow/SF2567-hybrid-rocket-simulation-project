import numpy as np

from equations.falloffs import sigmoid

import assumptions

def getDragCoefficient(mach):
  # José Carlos Santos (https://math.stackexchange.com/users/446262/jos%c3%a9-carlos-santos), How to find the exponential curve between two generic points?, URL (version: 2020-07-21): https://math.stackexchange.com/q/3276964
  c = assumptions.dragPeakMachNumber.get()

  Cleft = assumptions.dragLevelAtZero.get()
  Cright = assumptions.dragLevelMachAsymptote.get()
  Cpeak = assumptions.dragLevelAtPeak.get()

  leftRateConstant = 2
  rightRateConstant = assumptions.dragDropoffConstant.get()

  transformedLeftMach = leftRateConstant*(mach - c) / c
  transformedRightMach = -rightRateConstant*(mach - c) / c

  yValueAtLeft = np.exp(-leftRateConstant)
  yValueAtPeak = 1
  transformedLeftY = (np.exp(transformedLeftMach) - yValueAtLeft) / (yValueAtPeak - yValueAtLeft)
  leftDragCoefficient = Cleft + transformedLeftY * (Cpeak - Cleft)

  transformedRightY = np.exp(transformedRightMach) / yValueAtPeak
  rightDragCoefficient = Cright + transformedRightY * (Cpeak - Cright)

  leftRightFactor = sigmoid(mach, c, assumptions.dragPeakSmoothingRadius.get())

  return (1 - leftRightFactor) * leftDragCoefficient + (leftRightFactor) * rightDragCoefficient