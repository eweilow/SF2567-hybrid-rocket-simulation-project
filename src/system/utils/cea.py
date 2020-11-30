
import numpy as np
import hashlib
import os, sys

import assumptions

import options

# ugh https://stackoverflow.com/a/45669280
class HiddenPrints:
    def __enter__(self):
        self._original_stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout.close()
        sys.stdout = self._original_stdout

hiddenPrints = HiddenPrints()


class NasaCEA:
  cacheHits = 0
  totalHits = 0
  inMemoryCache = {}

  @staticmethod
  def clear():
    NasaCEA.cacheHits = 0
    NasaCEA.totalHits = 0
    del NasaCEA.inMemoryCache
    NasaCEA.inMemoryCache = {}

  def getPerformanceParameters(self, combustionPressure, ambientPressure, oxidizerTemperature, expansionRatio, mixtureRatio):
    import CoolProp.CoolProp as CP
    from rocketcea.cea_obj import CEA_Obj, add_new_fuel, add_new_oxidizer, add_new_propellant
    from rocketcea.cea_obj_w_units import CEA_Obj as CEA_Obj_W_Units
    from rocketcea.blends import makeCardForNewTemperature

    combustionPressure = max(10, combustionPressure) # Limitation of the library...
    fuelTemperature = 293

    HMOLAR = CP.PropsSI('HMOLAR','T',oxidizerTemperature,'P',combustionPressure,"N2O") / 1e3
    oxidizerDensity = CP.PropsSI('D','T',oxidizerTemperature,'P',combustionPressure,"N2O")
    RHO = oxidizerDensity * 1e3 / (1e2 * 1e2 * 1e2)

    oxid_card_str = """
    oxid=NITROUS wt=1.0 t,K={:.2f} h,kj/mol={:.2f} rho,g/cc={:.2f} N 2 O 1
    """.format(oxidizerTemperature, HMOLAR, RHO)
    fuel_card_str = """
    fuel=C(gr) wt={:.4f} t,K={:.2f}
    fuel=SASOL907 wt={:.4f} t,K={:.2f} h,kj/mol={:.2f} rho,g/cc={:.3f} C 50 H 102
    """.format(
      assumptions.carbonBlackFraction.get(), 
      fuelTemperature, 
      1 - assumptions.carbonBlackFraction.get(),
      fuelTemperature, 
      assumptions.fuelEnthalpyOfFormation.get(), 
      assumptions.fuelDensityLiquid.get() / 1e3)

    problemString = "Pc={:.2f},Pa={:.2f},To={:.2f},Tf={:.2f},ER={:.2f},MR={:.2f},oxid={:},fuel={:}".format(
      combustionPressure,
      ambientPressure,
      oxidizerTemperature,
      fuelTemperature,
      expansionRatio,
      mixtureRatio,
      oxid_card_str,
      fuel_card_str
    )
    strHash = hashlib.md5(problemString.encode()).hexdigest()
    NasaCEA.totalHits = NasaCEA.totalHits + 1
    if options.enableCeaLookup and strHash in NasaCEA.inMemoryCache:
      # print("Hit cache")
      NasaCEA.cacheHits = NasaCEA.cacheHits + 1
      return NasaCEA.inMemoryCache[strHash]
    

    # print(HMOLAR, RHO)

    add_new_oxidizer('NITROUS_COOLPROP', oxid_card_str)
    add_new_fuel('SASOLWAX907_CARBONBLACK', fuel_card_str)

    with hiddenPrints:
      combustionPressure = round(combustionPressure, 2)
      mixtureRatio = round(mixtureRatio, 2)
      expansionRatio = round(expansionRatio, 2)
      ambientPressure = round(ambientPressure, 2)

      cea = CEA_Obj_W_Units(
        oxName="NITROUS_COOLPROP", 
        fuelName="SASOLWAX907_CARBONBLACK", 
        pressure_units='Pa', 
        cstar_units='m/s', 
        temperature_units='K', 
        sonic_velocity_units='m/s', 
        enthalpy_units='J/kg', 
        density_units='kg/m^3', 
        specific_heat_units='J/kg-K', 
        viscosity_units='poise', 
        thermal_cond_units="W/cm-degC"
      )
      Isp,mode = cea.estimate_Ambient_Isp(Pc=combustionPressure, MR=mixtureRatio, eps=expansionRatio, Pamb=ambientPressure)
      IspVac, Cstar, Tc, MW, gamma = cea.get_IvacCstrTc_ThtMwGam(Pc=combustionPressure, MR=mixtureRatio, eps=expansionRatio)
      Cp = cea.get_Chamber_Cp(Pc=combustionPressure, MR=mixtureRatio, eps=expansionRatio)
      rho = cea.get_Chamber_Density(Pc=combustionPressure, MR=mixtureRatio, eps=expansionRatio)
      CfVac, Cf, mode = cea.get_PambCf(Pc=combustionPressure, MR=mixtureRatio, eps=expansionRatio)
      PcOvPe = cea.get_PcOvPe(Pc=combustionPressure, MR=mixtureRatio, eps=expansionRatio)
      exitPressure = combustionPressure / PcOvPe

      result = (Isp, Cp, MW, Cstar, Tc, gamma, rho, Cf, exitPressure, oxidizerDensity)
      if options.enableCeaLookup:
        NasaCEA.inMemoryCache[strHash] = result
      return result
