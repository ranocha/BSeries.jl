import time
import pybs
from pybs.rungekutta import methods as rk_methods

midpoint_method = rk_methods.RKmidpoint
up_to_order = 9
number_of_terms = pybs.unordered_tree.number_of_trees_up_to_order(up_to_order+1)

from itertools import islice
def first_values(f, n):
  return (f(tree) for tree in islice(pybs.unordered_tree.tree_generator(), 0, n))

start_time = time.time()
midpoint_series = midpoint_method.phi()
series = pybs.series.modified_equation(midpoint_series)
result = sum(first_values(series, number_of_terms))
end_time = time.time()
print(result)
print("", end_time - start_time, "seconds")

start_time = time.time()
midpoint_series = midpoint_method.phi()
series = pybs.series.modified_equation(midpoint_series)
result = sum(first_values(series, number_of_terms))
end_time = time.time()
print(result)
print("", end_time - start_time, "seconds")
