import utils.constants as constants
import subprocess
import hashlib
import os
import json

import numpy as np

inMemoryCache = {}
class NasaCEA:
  def compute(self, nozzleAreaRatio, carbonBlackFraction, oxidizerFuelRatio, pressure, fuelTemperature, oxidizerTemperature, oxidizerDensity):
    inputStr = """
    problem  rocket  equilibrium 
      o/f={0:4f} 
      p,bar={1:4f} 
      pi/p={1:4f} 
      sup-ae/at={2:4f}
    reactants   
      oxid=N2O wt=1 t,K={5:4f} rho,g/cc = {7:4f}
      fuel=C(gr) wt={3:4f} t,K={5:4f}
      fuel=SASOL907 wt={4:4f} t,k={6:4f} h,kj/mol=-1438.200 rho,g/cc=0.720 C 50 H 102
    output
      short
      siunits
    end
    """.format(
        oxidizerFuelRatio, 
        pressure / constants.Pressure.bar, 
        nozzleAreaRatio, 
        carbonBlackFraction, 
        1 - carbonBlackFraction,
        fuelTemperature,
        oxidizerTemperature,
        oxidizerDensity / 1000
      )

    strHash = hashlib.md5(inputStr.encode()).hexdigest()

    cacheFile = "G:/ceaCache/{0}".format(strHash)

    try:
      if cacheFile in inMemoryCache:
        return inMemoryCache[cacheFile]

      if os.path.exists(cacheFile):
        with open(cacheFile, 'r') as file:
          inMemoryCache[cacheFile] = json.loads(file.read())
          return inMemoryCache[cacheFile]
    except:
      os.remove(cacheFile)
      
    f = open("./CEA/temp.inp", "w")
    f.write(inputStr)
    f.close()

    p = subprocess.Popen("FCEA2m.exe", shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, cwd="./CEA")
    try:
      p.communicate("temp\n".encode(), timeout=1)
    except:
      p.kill()

    def parseEntry(entry):
      if entry == "":
        return None
        
      if(entry[-2] == " " or entry[-2] == "+"):
        exponent = int(entry[-1])
        return float(entry[0:-2]) * pow(10, exponent)
      elif(entry[-2] == "-"):
        exponent = -int(entry[-1])
        return float(entry[0:-2]) * pow(10, exponent)
      else:
        return float(entry)

    def parseLine(line):
      a = line[17:25].strip()
      b = line[26:34].strip()
      c = line[35:43].strip()
      d = line[44:52].strip()
      return [parseEntry(a), parseEntry(b), parseEntry(c), parseEntry(d)]
    
    def mapEntries(array):
      return {
        "pressure": array[1] * constants.Pressure.bar,
        "temperature": array[2],
        "density": array[3],
        "H": array[4]*1000,
        "U": array[5]*1000,
        "G": array[6]*1000,
        "S": array[7]*1000,
        "M": array[8],
        "Cp": array[9]*1000,
        "gamma": array[10],
        "speedOfSound": array[11],
        "mach": array[13],
        "cStar": array[14],
        "CF": array[15],
        "Ivac": array[16],
        "Isp": array[17]
      }

    lines = []
    with open("./CEA/temp.out", 'r') as f:
      for line in f.readlines():
        if line.startswith(" Pinf/P"):
          lines.append(parseLine(line))
        if line.startswith(" P, BAR"):
          lines.append(parseLine(line))
        if line.startswith(" T, K"):
          lines.append(parseLine(line))
        if line.startswith(" RHO, KG/CU M"):
          lines.append(parseLine(line))
        if line.startswith(" H, KJ/KG"):
          lines.append(parseLine(line))
        if line.startswith(" U, KJ/KG"):
          lines.append(parseLine(line))
        if line.startswith(" G, KJ/KG"):
          lines.append(parseLine(line))
        if line.startswith(" S, KJ/(KG)(K)"):
          lines.append(parseLine(line))
        if line.startswith(" M, (1/n)"):
          lines.append(parseLine(line))
        if line.startswith(" Cp, KJ/(KG)(K)"):
          lines.append(parseLine(line))
        if line.startswith(" GAMMAs"):
          lines.append(parseLine(line))
        if line.startswith(" SON VEL,M/SEC"):
          lines.append(parseLine(line))
        if line.startswith(" MACH NUMBER"):
          lines.append(parseLine(line))
        if line.startswith(" Ae/At"):
          lines.append(parseLine(line))
        if line.startswith(" CSTAR, M/SEC"):
          lines.append(parseLine(line))
        if line.startswith(" CF"):
          lines.append(parseLine(line))
        if line.startswith(" Ivac, M/SEC"):
          lines.append(parseLine(line))
        if line.startswith(" Isp, M/SEC"):
          lines.append(parseLine(line))
      result = np.array(lines)
      
      pressureDistance = np.abs(result[1,2:] - 1*constants.Pressure.bar)
      [atmosphericConditionsIndex] = np.where(pressureDistance == np.amin(pressureDistance))
      
      areaRatioDistance = np.abs(result[13,2:] - nozzleAreaRatio)
      [areaRatioConditionsIndex] = np.where(areaRatioDistance == np.amin(areaRatioDistance))

      output = {
        "chamber": mapEntries(result[:,0]),
        "throat": mapEntries(result[:,1]),
        "atmospheric": mapEntries(result[:,atmosphericConditionsIndex[0] + 2]),
        "exhaust": mapEntries(result[:,areaRatioConditionsIndex[0] + 2]),
      }
      
      with open(cacheFile, 'w') as file:
        inMemoryCache[cacheFile] = output
        file.write(json.dumps(output))

      return output

# cea = NasaCEA()
# print(cea.compute(4, 0.02, 8, 5e5, 293, 273))