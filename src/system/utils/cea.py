import CoolProp.CoolProp as CP
import numpy as np
from rocketcea.cea_obj import CEA_Obj, add_new_fuel, add_new_oxidizer, add_new_propellant
from rocketcea.cea_obj_w_units import CEA_Obj as CEA_Obj_W_Units
from rocketcea.blends import makeCardForNewTemperature
import os, sys

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
  def getPerformanceParameters(self, combustionPressure, ambientPressure, oxidizerTemperature, expansionRatio, mixtureRatio):
    combustionPressure = max(10, combustionPressure) # Limitation of the library...
    
    HMOLAR = CP.PropsSI('HMOLAR','T',oxidizerTemperature,'P',combustionPressure,"N2O") / 1e3
    RHO = CP.PropsSI('D','T',oxidizerTemperature,'P',combustionPressure,"N2O") * 1e3 / (1e2 * 1e2 * 1e2)

    # print(HMOLAR, RHO)

    card_str = """
    oxid=NITROUS wt=1.0 t,K={:.10f} h,kj/mol={:.10f} rho,g/cc={:.10f} N 2 O 1
    """.format(oxidizerTemperature, HMOLAR, RHO)
    add_new_oxidizer('NITROUS_COOLPROP', card_str.format(oxidizerTemperature))

    fuelTemperature = 293
    card_str = """
    fuel=C(gr) wt=0.02 t,K={:.4f}
    fuel=SASOL907 wt=0.98 t,K={:.4f} h,kj/mol=-1438.200 rho,g/cc=0.720 C 50 H 102
    """.format(fuelTemperature, fuelTemperature)
    add_new_fuel('SASOLWAX907_CARBONBLACK', card_str)

#    cea = CEA_Obj(
#      oxName="NITROUS_COOLPROP", 
#      fuelName="SASOLWAX907_CARBONBLACK"
#    )
#    s = cea.get_full_cea_output(Pc=combustionPressure/1e5, MR=mixtureRatio, eps=expansionRatio, pc_units="bar", output='siunits')
#    print(s)
#
#    pass

    with hiddenPrints:
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
      # cea = CEA_Obj(
      #   oxName=oxidizerCard, 
      #   fuelName="SASOLWAX907_CARBONBLACK"
      # )
      # s = cea.get_full_cea_output(Pc=combustionPressure/1e5, MR=mixtureRatio, eps=expansionRatio, pc_units="bar", output='siunits')
      # print(s)

      exitPressure = combustionPressure / PcOvPe

      return Isp, Cp, MW, Cstar, Tc, gamma, rho, Cf, exitPressure
