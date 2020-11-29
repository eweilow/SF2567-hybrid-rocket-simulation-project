class Lengths:
  m = 1
  dm = m/10
  cm = dm/10
  mm = cm/10

class Volume:
  m3 = 1
  dm3 = m3/(10**3)
  cm3 = dm3/(10**3)
  liter = dm3

class Pressure:
  kPa = 1e3
  MPa = 1e6
  bar = 100 * kPa

class Area:
  cm2 = Lengths.cm * Lengths.cm