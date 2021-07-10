# Qualification - Problem 2.4
This was, by far, my favorite problem of the qualification.
As you can see from my very messy code, I used several methods to solve the various problems (three different approaches currently exist in the code).

I was not able to find a way to avoid iterating over all the terms in the patient or over all the candidates diseases in the list.
Because of that limitation, it means the problem is figuring out how to _efficiently_ find the LCA for each phenotype in the patient compared to some set of terms in a disease.

For sub-problems 1-4, I basically used a sparse table built on a traversal algorithm.  It's quite messy, but is loosely based on concepts found here: "Method 3 (Sparse Table Algorithm)" from https://www.geeksforgeeks.org/range-minimum-query-for-static-array/.  I don't really recommend it.

For sub-problem 5, this was a special case of the general algorithm where the number of phenotype terms for each patient/disease was exactly 1.
Given this condition, you can actually assigned a disease ID to it's single phenotype and every parent of that phenotype.
Then, solving the best disease for a patient is just iterating up from the patient phenotype until you encounter a node with any assigned disease.

For sub-problems 6 and 7, the scale of the problems drastically increase and we don't have the special case from problem 5.
However, they _are_ still limited with a max set of six different phenotypes per patient and disease.
My solution pre-computed a sparse matrix where each row was a phenotype and each column a disease.
Then, for each set of patient phenotypes, I created a heap queue that initial contained exactly the patient phenotypes and was sorted such that the highest IC would pop off the queue first.
For each popped set of phenotypes, we then compare to the pre-computed matrix to see if any diseases have _all_ of those terms.
If any do, you have found the highest IC disease for the patient phenotypes.
If not, then for each phenotype in the popped set, it creates a new phenotype set where only that phenotype is modified to its parent, and the IC is calculated for that new set.
The process then repeats.

The end result is that the phenotype terms for the patient progressively modified one at a time (by moving up the tree) until the closest matching disease is found, and that disease is guaranteed to be the highest IC due to the heapq ordering.
