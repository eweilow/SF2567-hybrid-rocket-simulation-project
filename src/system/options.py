printTime = False

printFailedInjector = False


enableTankWallHeatTransfer = True

enableCeaLookup = True

enableTankTemperatureInterpolation = True
tankInterpolant = 5 # 5th order
tankInterpolantPointCount = 1000

currentlySolvingWithDAE = False

enableDAESolver = True

solveTankVolumeContraintWithDAE = True

rootFindingType = "brentq"
# rootFindingType = "brentq"
# rootFindingType = "newton"
# rootFindingType = "bisect"

combustion_rtol = 1e-4
combustion_atol = 1e-4

flight_rtol = 1e-4
flight_atol = 1e-4