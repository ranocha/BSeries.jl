import os
directory = os.path.dirname(os.path.abspath(__file__))
io_file = os.path.join(directory, "benchmark_python_bseries.txt")

import sys
from importlib.metadata import version
with open(io_file, 'w') as io:
  print("Python version", sys.version, file=io)
  print("Package version", version('pybs'), file=io)

import time
import BSeries.bs as bs
import nodepy.runge_kutta_method as rk

midpoint_method = rk.loadRKM("Mid22")
up_to_order = 9

start_time = time.time()
series = bs.modified_equation(None, None,
                              midpoint_method.A, midpoint_method.b,
                              up_to_order, True)
result = sum(series.values())
end_time = time.time()
with open(io_file, 'a') as io:
  print(result, file=io)
  print("", end_time - start_time, "seconds", file=io)

start_time = time.time()
series = bs.modified_equation(None, None,
                              midpoint_method.A, midpoint_method.b,
                              up_to_order, True)
result = sum(series.values())
end_time = time.time()
with open(io_file, 'a') as io:
  print(result, file=io)
  print("", end_time - start_time, "seconds", file=io)
