# Finals - Problem 3.2
This problem is about finding who is most likely to infect the most people (i.e. who is the biggest superspreader) in a community.
I have two versions of a solution, neither of which were perfect solutions to the problem.

## Initial approach (fast version)
My initial approach (`solveProblem(...)`) involved working _backwards_ through the timeline and assigning increasing scores to each individual based on their future contacts.  
This method was extremely inaccurate, primarily because there was no concept of whether an individual had already been infected or not.
However, it provided a decent estimate and got points on the board.

## Second approach (slow version)
I later revisited the problem to try and refine my answer using a slower, forward approach (problems 5-7, mainly).
This method actually used the initial approach to prioritize candidates that were _likely_ to be superspreaders.
Then, I did a calculation of the probability of infections with the best _X_ candidates (to limit runtime).
The most infectious of that candidate list was then used as the result.
This method was certainly imperfect, but did decently well.
I suspect there may be a subtle bug (or nuance I missed in the problem) that led to imperfect answers for problems 5 and 6.
For problem 7, I instead suspect that my approach as estimating the infection simply did not work well under those particular conditions.
