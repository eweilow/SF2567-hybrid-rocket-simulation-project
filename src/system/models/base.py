class Model:
  def derivativesDependsOn(self, models):
    return []

  def derivedVariablesDependsOn(self, models):
    return []
    
  def initializeSimplifiedModel(self, timeHistory, stateHistory, derivedVariablesHistory):
    return False

  def computeSimplifiedState(self, args, time):
    return []

  def computeSimplifiedDerivedVariables(self, args, time):
    return []

  def initializeState(self):
    return []

  def computeDerivatives(self, t, state, derived, models):
    return []

  def computeDerivedVariables(self, t, state):
    return []

