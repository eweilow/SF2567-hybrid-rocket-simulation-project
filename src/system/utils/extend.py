import numpy as np

def extend(t, y):
  dt = t[-1] - t[-2]
  deriv = (y[-1] - y[-2]) / dt
  A = y[-1]
  k = deriv / (-A)

  # print("Extending...\nderiv = {:.4f}, A = {:.4f}, k = {:.4f}".format(deriv, A, k))
  # print("  t0 = {:.4f}, t1 = {:.4f}, y0 = {:.4f}, y1 = {:.4f}".format(t[-2], t[-1], y[-2], y[-1]))

  def get(T):
    return A * np.exp(-k *  (T - t[-1]))

  return get

def extendDerivative(t, y):
  dt = t[-1] - t[-2]
  deriv = (y[-1] - y[-2]) / dt
  A = y[-1]
  k = deriv / (-A)

  # print("Extending...\nderiv = {:.4f}, A = {:.4f}, k = {:.4f}".format(deriv, A, k))
  # print("  t0 = {:.4f}, t1 = {:.4f}, y0 = {:.4f}, y1 = {:.4f}".format(t[-2], t[-1], y[-2], y[-1]))

  def get(T):
    return -k * A * np.exp(-k *  (T - t[-1]))

  return get