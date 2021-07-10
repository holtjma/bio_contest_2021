# Finals - Problem 3.2
This problem primarily focuses on genotype imputation, which means we are trying to predict missing genotypes from an individual based on a limited set of observations. 
The problem is complicated by providing a set of known haplotype pairs in the population, but provided observed diplotypes (e.g. 0, 1, 2) for the predictions.

## Initial approach (`main.py`)
My initial approach was to find the "best" combination of haplotypes that was closest to the observations and then inferring the missing by simply adding the best haplotype combination together.  
However, it became fairly obvious early in the contest that this approach was flawed.
Specifically, this approach didn't account for the concept of recombination (this is not _explicitly_ mentioned in the problem description either).
However, since the problem was simulated from the sim1000G package, it seemed likely to be relevant, so...

## Final approach (`main_v2.py`)
I tried a second approach based on machine learning.
The main idea was to turn all of the known haplotype combinations into diplotypes for training.
Then, for each missing locus, I trained a model to predict the diplotype value of that locus based on the nearby loci.
I generally found that adding more loci to the training (parameter `NUM_FEAT`), increased performance.
I also found that `GradientBoostingClassifier` tended to perform better than the `RandomForestClassifier` (this is generally expected in my experience).
Looking back, the solution (while clearly imperfect) is relatively simple.  
Strangely, it also seems to perform quite well when compared to the other comparators (with a cursory glance, it seems like the second best algorithm in the competition by scoring).
