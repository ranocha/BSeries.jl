import os
directory = os.path.dirname(os.path.abspath(__file__))
io_file = os.path.join(directory, "benchmark_python_pybs.txt")

import sys
from importlib.metadata import version
with open(io_file, 'w') as io:
  print("Python version", sys.version, file=io)
  print("Package version", version('pybs'), file=io)

import time
import pybs
from pybs.rungekutta import methods as rk_methods

midpoint_method = rk_methods.RKmidpoint
up_to_order = 9
number_of_terms = pybs.unordered_tree.number_of_trees_up_to_order(up_to_order+1)

from itertools import islice
def first_values(f, n):
  return (f(tree) for tree in islice(pybs.unordered_tree.tree_generator(), 0, n))


with open(io_file, 'a') as io:
  print("\nModified equation", file=io)

start_time = time.time()
midpoint_series = midpoint_method.phi()
series = pybs.series.modified_equation(midpoint_series)
result = sum(first_values(series, number_of_terms))
end_time = time.time()
with open(io_file, 'a') as io:
  print(result, file=io)
  print("", end_time - start_time, "seconds", file=io)

start_time = time.time()
midpoint_series = midpoint_method.phi()
series = pybs.series.modified_equation(midpoint_series)
result = sum(first_values(series, number_of_terms))
end_time = time.time()
with open(io_file, 'a') as io:
  print(result, file=io)
  print("", end_time - start_time, "seconds", file=io)


with open(io_file, 'a') as io:
  print("\nEnergy preservation", file=io)

start_time = time.time()
a = pybs.series.AVF
b = pybs.series.modified_equation(a)
result = pybs.series.energy_preserving_upto_order(b, up_to_order)
end_time = time.time()
with open(io_file, 'a') as io:
  print(result, file=io)
  print("", end_time - start_time, "seconds", file=io)

start_time = time.time()
a = pybs.series.AVF
b = pybs.series.modified_equation(a)
result = pybs.series.energy_preserving_upto_order(b, up_to_order)
end_time = time.time()
with open(io_file, 'a') as io:
  print(result, file=io)
  print("", end_time - start_time, "seconds", file=io)


with open(io_file, 'a') as io:
  print("\nSymplecticity (conservation of quadratic invariants)", file=io)

start_time = time.time()
a = rk_methods.RKimplicitMidpoint.phi()
result = pybs.series.symplectic_up_to_order(a, up_to_order)
end_time = time.time()
with open(io_file, 'a') as io:
  print(result, file=io)
  print("", end_time - start_time, "seconds", file=io)

start_time = time.time()
a = rk_methods.RKimplicitMidpoint.phi()
result = pybs.series.symplectic_up_to_order(a, up_to_order)
end_time = time.time()
with open(io_file, 'a') as io:
  print(result, file=io)
  print("", end_time - start_time, "seconds", file=io)
