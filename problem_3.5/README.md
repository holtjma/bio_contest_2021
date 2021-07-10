# Finals - Problem 3.5
This problem is really a two part problem:
1. Figure out if a datasets has one or two viral populations and...
2. If there are two, figure out the minor sequence

For this problem in particular, I believe my approach is suboptimal (although still close to the best submitted answers).
So with that in mind...

## My approach
1. Calculate a Hamming distance matrix for each pair of viral sequences
2. Use those distances plus a clustering algorithm to separate the sequences into two groups (`AgglomerativeClustering` seemed to be the only one I got to work)
3. Generate a consensus sequence for each of the two groups
4. If the number of differences between the two sequences was less than the provided threshold (`k` in the problem), then only one population was considered present
5. Compare the fraction of reads in the smaller population and the total "error rate" to some tuned numbers.  If the fraction of reads is too low and the number of presumed "errors" is low, then it was marked as only having one sequence. 
6. Otherwise, we assumed two sequences were present and returned the other consensus.

The exact values for the population fraction and error rates are basically "magic" numbers that I tuned by hand until I got a relatively optimal answer.  The _correct_ thing to do likely factors in the actual statistical probabilities provided by the problem, but I was not able to figure this out during the 24hrs.
