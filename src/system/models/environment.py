from models.base import Model

import assumptions

class EnvironmentModel(Model):
  def derivativesDependsOn(self, models):
    return []

  def derivedVariablesDependsOn(self, models):
    return []
  
  derived_ambientPressure = 0

  def initializeState(self):
    return [0]

  def computeDerivatives(self, t, state, derived, models):
    return [0]

  def computeDerivedVariables(self, t, state, models):
    return [assumptions.initialAtmosphericPressure.get()]