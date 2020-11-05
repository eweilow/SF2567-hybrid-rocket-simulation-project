
import matplotlib.pyplot as plt
import numpy as np

fig = plt.figure()

legendEntries = []
with open('./tmp/combustionProperties.npy', 'rb') as f:
  mixtureRatios = np.load(f)
  temperatures = np.load(f) - 273.15
  chamberPressure = np.load(f)
  ambientPressure = np.load(f)

  for temperature in temperatures:
    Isp = np.load(f)
    CpAve = np.load(f)
    MolWt = np.load(f)
    Cstar = np.load(f)
    Tc = np.load(f)
    gamma = np.load(f)
    
    plt.subplot(2, 2, 1)
    plt.plot(mixtureRatios, Isp)

    plt.subplot(2, 2, 2)
    plt.plot(mixtureRatios, Cstar)

    plt.subplot(2, 2, 3)
    plt.plot(mixtureRatios, Tc)

    plt.subplot(2, 2, 4)
    plt.plot(mixtureRatios, gamma)

    legendEntries.append("{:.2f} °C".format(temperature))


fig.suptitle("Nitrous Oxide / SASOLWAX 907 + 2% Carbon Black, Pc = {:.2f} bar, Pa = {:.2f} bar".format(chamberPressure / 1e5, ambientPressure / 1e5))

plt.subplot(2, 2, 1)
plt.title("Specific Impulse")
plt.xlabel("Oxidizer-Fuel Ratio")
plt.ylabel("Specific Impulse [s]")
plt.grid()
plt.legend(legendEntries)

plt.subplot(2, 2, 2)
plt.title("Characteristic Velocity")
plt.xlabel("Oxidizer-Fuel Ratio")
plt.ylabel("Cstar [m/s]")
plt.grid()
plt.legend(legendEntries)

plt.subplot(2, 2, 3)
plt.title("Combustion Temperature")
plt.xlabel("Oxidizer-Fuel Ratio")
plt.ylabel("Temperature [K]")
plt.grid()
plt.legend(legendEntries)

plt.subplot(2, 2, 4)
plt.title("Gamma")
plt.xlabel("Oxidizer-Fuel Ratio")
plt.ylabel("Gamma")
plt.grid()
plt.legend(legendEntries)

plt.show()