from equations import drag

import matplotlib.pyplot as plt
import numpy as np
import assumptions

machValues = np.linspace(0, 20, 200)

dragValues = drag.getDragCoefficient(machValues)

plt.plot(machValues, dragValues)
plt.vlines(assumptions.dragPeakMachNumber.get(), np.min(dragValues), np.max(dragValues), linestyles='dashed')
plt.hlines(assumptions.dragLevelAtZero.get(), np.min(machValues), np.max(machValues), linestyles='dashed')
plt.hlines(assumptions.dragLevelMachAsymptote.get(), np.min(machValues), np.max(machValues), linestyles='dashed')
plt.hlines(assumptions.dragLevelAtPeak.get(), np.min(machValues), np.max(machValues), linestyles='dashed')
plt.show()