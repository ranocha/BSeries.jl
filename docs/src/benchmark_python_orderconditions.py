import os
directory = os.path.dirname(os.path.abspath(__file__))
io_file = os.path.join(directory, "benchmark_python_orderconditions.txt")

import sys
with open(io_file, 'w') as io:
  print("Python version", sys.version, file=io)

import time
from orderConditions import BSeries
import nodepy.runge_kutta_method as rk

midpoint_method = rk.loadRKM("Mid22")
up_to_order = 9


with open(io_file, 'a') as io:
  print("\nModified equation", file=io)

start_time = time.time()
BSeries.set_order(up_to_order)
Y1 = BSeries.y() + midpoint_method.A[1,0] * BSeries.hf()
rk2 = BSeries.y() + midpoint_method.b[0] * BSeries.hf() + midpoint_method.b[1] * BSeries.compo_hf(Y1)
series = BSeries.modified_equation(rk2)
result = series.sum()
end_time = time.time()
with open(io_file, 'a') as io:
  print(result, file=io)
  print("", end_time - start_time, "seconds", file=io)

start_time = time.time()
BSeries.set_order(up_to_order)
Y1 = BSeries.y() + midpoint_method.A[1,0] * BSeries.hf()
rk2 = BSeries.y() + midpoint_method.b[0] * BSeries.hf() + midpoint_method.b[1] * BSeries.compo_hf(Y1)
series = BSeries.modified_equation(rk2)
result = series.sum()
end_time = time.time()
with open(io_file, 'a') as io:
  print(result, file=io)
  print("", end_time - start_time, "seconds", file=io)


print("\nModifying integrator")

start_time = time.time()
BSeries.set_order(up_to_order)
Y1 = BSeries.y() + midpoint_method.A[1,0] * BSeries.hf()
rk2 = BSeries.y() + midpoint_method.b[0] * BSeries.hf() + midpoint_method.b[1] * BSeries.compo_hf(Y1)
series = BSeries.modifying_integrator(rk2)
result = series.sum()
end_time = time.time()
with open(io_file, 'a') as io:
  print(result, file=io)
  print("", end_time - start_time, "seconds", file=io)

start_time = time.time()
BSeries.set_order(up_to_order)
Y1 = BSeries.y() + midpoint_method.A[1,0] * BSeries.hf()
rk2 = BSeries.y() + midpoint_method.b[0] * BSeries.hf() + midpoint_method.b[1] * BSeries.compo_hf(Y1)
series = BSeries.modifying_integrator(rk2)
result = series.sum()
end_time = time.time()
with open(io_file, 'a') as io:
  print(result, file=io)
  print("", end_time - start_time, "seconds", file=io)
