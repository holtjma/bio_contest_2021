# Qualification - Problem 2.3
This problem is basically to find the best sum of two values from different sets (metabolites and adducts) that gets closest to the given value.
The naive way (see my initial solution of `solveProblem_old`) is to calculate _every_ possible combination, sort that list, and then search for the closest match in the list.
Unfortunately, as the list of metabolites and adducts gets larger, this list becomes too big to efficiently calculate and store in memory.

My revised solution does not pre-calculate anything at all.
Instead, for each observation (o), it iterates through the list of adducts (which was fixed to K<=10000), and calculates the best possible value for the metabolites to match the observation (e.g. `o-k`) and then searches the sorted list of metabolites to find the best possible metabolite that is available.
While the number of metabolites could get quite large (max 10**6), it is very efficient to search that list for the best possible value
