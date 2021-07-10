# Finals - Problem 3.3
This problem primarily focuses on the classic problem of finding the "causitive" mutation with a set of affected and unaffected individuals.  
We were also limited in the number of submissions we could attempt before getting locked out.
The classic solution would be to perform a GWAS (or similar) study, but I am unfortunately not too familiar with the details of how exactly that's done.
So the following is my attempt at a "GWAS-like" solution, that surprisingly worked.

## Easy problem (`main.py`)
My solution was to build a series of machine learning models based on adjacent haploid observations.
The number of observations used was tunable (see parameter `k`).
I attempted to train a model to predict whether a patient was affected or not based on those `k` observations.
For many problems, this very quickly led to an obvious peak in the accuracy of the models when plotted (see comment block in the code that is currently disabled).
However, this peak was not always correct, but it was usually in the "ballpark" of correct.
Subsequent attempts usually involved refining my answers until all of the answers were correct.

## Hard problem (`main_hard.py`)
The solution for this problem was practically identical, just with twice the input feature due to the diploid nature of the inputs.
The only other difference was refinement of the final answers took a little longer because the peak was not always as blatantly obvious as with the haploid problem.
