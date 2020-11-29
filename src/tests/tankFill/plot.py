
import matplotlib.pyplot as plt
import numpy as np
import constants

with open('./tmp/tankFilling.npy', 'rb') as f:
  t = np.load(f)
  y = np.load(f)

plt.figure(figsize=(12, 8))
plt.subplot(2,3,1)
plt.plot(t, y[:,0])
plt.plot(t, y[:,8])
plt.plot(t, y[:,7])
#plt.plot(t, y[:,1])
#plt.plot(t, y[:,0] + y[:,1])
plt.legend(("Total", "Liquid", "Gas"))
plt.xlabel("Time (s)")
plt.ylabel("Mass (kg)")
plt.title("Mass in tank")

plt.subplot(2,3,2)
plt.plot(t, y[:,2] - 273.15)
#plt.plot(t, y[:,3] - 273.15)
#plt.legend(("Liquid", "Gas"))
plt.xlabel("Time (s)")
plt.ylabel("Temperature (°C)")
plt.title("Temperatures in tank")

plt.subplot(2,3,3)
plt.plot(t, y[:,1] / 1e3)
#plt.plot(t, y[:,5] / constants.Volume.liter)
#plt.plot(t, (y[:,4] + y[:,5]) / constants.Volume.liter)
#plt.legend(("Liquid", "Gas", "Total"))
plt.xlabel("Time (s)")
plt.ylabel("Energy (kJ)")
plt.title("Energy in tank")

plt.subplot(2,3,4)
plt.plot(t, y[:,3] / constants.Pressure.bar)
plt.xlabel("Time (s)")
plt.ylabel("Pressure (bar)")
plt.title("Pressures in tank")

plt.subplot(2,3,5)
plt.plot(t, y[:,4] * 100)
plt.xlabel("Time (s)")
plt.ylabel("Liquid level (%)")
plt.title("Liquid level in tank")
# plt.subplot(2,3,3)
# plt.plot(t, y[:,4])
# plt.xlabel("Time (s)")
# plt.ylabel("Temperature (K)")
# plt.title("Temperature in tank")
# 
# plt.subplot(2,3,4)
# plt.plot(t, y[:,5] / constants.Pressure.bar)
# plt.xlabel("Time (s)")
# plt.ylabel("Pressure (bar)")
# plt.title("Pressure in tank")

plt.subplot(2,3,6)
plt.plot(t, y[:,5])
plt.plot(t, y[:,6])
plt.legend(("Liquid entered", "Gas exited"))
#plt.plot(t, y[:,1])
#plt.plot(t, y[:,0] + y[:,1])
#plt.legend(("Liquid", "Gas", "Total"))
plt.xlabel("Time (s)")
plt.ylabel("Mass (kg)")
plt.title("Inlet & outlet masses")

plt.tight_layout()
plt.subplots_adjust(
  left=0.075,
  bottom=0.075,
  top=1-0.075,
  right=1-0.075,
  wspace=0.25,
  hspace=0.35
)
plt.show()