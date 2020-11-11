import numpy as np

"""
Create a matrix (M) + inverse matrix (invM) pair that represents the mapping between full system and model:
y* = M y
y = invM y*

Here we have the variables:
- y is the state of the full system
- y* is the state of the model in question
"""
def makeSystemMatrices(nextIndex, modelLength, fullSystemLength):
  matrix = np.zeros((fullSystemLength, modelLength))
  for i in range(modelLength):
    matrix[nextIndex][i] = 1
    nextIndex = nextIndex + 1
  invMatrix = np.transpose(matrix)

  return matrix, invMatrix, nextIndex
