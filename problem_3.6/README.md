# Finals - Problem 3.6
Compared to all the other problems, this would seemed to be the largest in scale.
The naive solution for both problem types is to compare each read to all the isoforms.
However, for the large problems, this became an extremely expensive task and made the naive approach borderline impossible without multi-processing.
I tried a couple methods for pruning the list of isoforms that needed to be compared, but none of them were particularly effective due to the unconstrained lengths of isoforms/reads.
I instead ended up taking a parallel approach when calculating the solutions, solving each problem individually and joining the results when all the computation was completed.
The results in this repo represent the single-threaded versions of my solutions, but parallelization is basically just restricting the execution to a subset of the problems.

I suspect there is some algorithmic trick to solving this more efficiently, or perhaps using something other than python for the language may make single-threaded feasible.
If someone has an efficient approach that doesn't _require_ parallelization to finish in the contest time, I'd be interested in seeing the approach.

## Easy problem (`main.py`)
This problem is to find the smallest index isoform that matches the read, and then to also report the total number of isoforms that could match the read.  
Naively, this seems to require comparing each read to each isoform if for no other purpose than the total count.

My approach was to loop over each isoform in an outer loop, and then each read in an inner loop.
The comparison involved doing a binary search in the isoform for the closest start exon, and then iterating through the rest of the read verifying that the exon endpoints overlapped sufficiently as described by the problem.
I then set the lowest isoform for that read if it was unset, and incremented the total isoform count for that read.

This parallelizes over each isoform and can be combined by just keeping the minimum isoform for all solved subproblems and adding up the total number of matching problems across all isoforms also.

## Hard problem (`main_v3.py`)
In contrast to the "easy" problem, my outer loop iterated over the reads and compared to each isoform in an inner loop.
I then calculated an overlap score for the exons and introns (same function with different inputs, `overlapRanges(...)`).
Combining the result and then just keep the best one for the read to solve the problem.
This trivially parallelizes by treating each read (or group of reads) independently.
