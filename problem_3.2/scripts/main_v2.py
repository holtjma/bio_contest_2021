
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

import os
import numpy as np
import json

from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import GradientBoostingClassifier

PROBLEM_DIR = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
DATA_DIR = f'{PROBLEM_DIR}/data'
RESULTS_DIR = f'{PROBLEM_DIR}/results'

def loadProblems(fn):
    ret = {}
    fp = open(fn, 'r')
    
    num_known, num_partial = [int(x) for x in fp.readline().rstrip().split(' ')]
    
    observed_dips = []
    for x in range(0, num_known):
        fp.readline()
        new_hap = fp.readline().rstrip()
        new_hap2 = fp.readline().rstrip()
        #I don't think this matters
        #assert(new_hap not in known_haps)
        #assert(new_hap2 not in known_haps)
        c1 = np.array([int(c) for c in new_hap], dtype='<i1')
        c2 = np.array([int(c) for c in new_hap2], dtype='<i1')
        observed_dips.append(c1+c2)

    unknown_haps = []
    for x in range(0, num_partial):
        fp.readline()
        unknown_haps.append(fp.readline().rstrip())

    fp.close()

    ret['observed_dips'] = observed_dips
    ret['unknown_haps'] = unknown_haps

    return ret

def solveProblem(problem):
    #Load the problem dict (or list of problem dicts)
    #metabolites = problem['metabolites']
    observed_dips = problem['observed_dips']
    unknown_haps = problem['unknown_haps']

    hap_len = len(observed_dips[0])
    print(f'hap_len: {hap_len}')
    print(f'num_known: {len(observed_dips)}')
    print(f'num_unknown: {len(unknown_haps)}')
    num_unknown = len(unknown_haps)

    print('Preprocessing...')
    #figure out where our unknowns are
    known_indices = [(c != '?') for c in unknown_haps[0]]
    #print(known_indices)

    unknown_ints = []
    for dip in unknown_haps:
        unknown_ints.append(
            np.array([int(c) if c != '?' else 3 for c in dip], dtype='<i1')
        )
    
    results = []
    for r in unknown_haps:
        results.append('')

    #NUM_FEAT = 31
    NUM_FEAT = 127
    for ind, ki in enumerate(known_indices):
        print('\t', ind, ki, len(known_indices))
        if ki:
            #this one is known, so just add the value
            for x in range(0, num_unknown):
                results[x] += unknown_haps[x][ind]
        else:
            #find the closest indices
            delta = 1
            found = []
            while len(found) < NUM_FEAT and (ind - delta > 0 or ind+delta < hap_len):
                if ind - delta > 0 and known_indices[ind-delta]:
                    found.append(ind-delta)
                if ind+delta < hap_len and known_indices[ind+delta]:
                    found.append(ind+delta)
                delta += 1
            
            lookups = np.array(found)
            
            #print(f'found for {ind}:', lookups)

            #fit a model
            xs = [od[lookups] for od in observed_dips]
            ys = [od[ind] for od in observed_dips]
            #forest = RandomForestClassifier(random_state=0)
            y_sum = np.sum(ys)
            if y_sum == 0 or y_sum == 2*len(ys):
                predictions = [ys[0]]*len(unknown_ints)
            else:
                forest = GradientBoostingClassifier(
                    random_state=0, max_depth=6, n_estimators=40, subsample=0.5
                )
                forest.fit(xs, ys)

                #now do some predictions
                xp = [ui[lookups] for ui in unknown_ints]
                predictions = forest.predict(xp)

            for x in range(0, num_unknown):
                results[x] += str(predictions[x])

    '''
    #now do the final imputations
    imputed = []
    for j, c in enumerate(uh):
        if c == '?':
            h1 = int(hap_arr[best_x][j])
            h2 = int(hap_arr[best_y][j])
            imputed.append(str(h1+h2))
        else:
            imputed.append(c)
    results.append(''.join(imputed))
    '''

    return results
    
def writeResults(fn, all_results):
    fp = open(fn, 'w+')
    for result in all_results:
        fp.write(result+'\n')
    fp.close()

if __name__ == '__main__':
    #there are usually multiple per problem
    starting_problem = 3
    ending_problem = 3

    if not os.path.exists(RESULTS_DIR):
        os.makedirs(RESULTS_DIR)

    #go through each ones
    for problem in range(starting_problem, ending_problem+1):
        #filenames below might need to change per problem
        print(f'Analyzing problem set #{problem}...')
        fn = f'{DATA_DIR}/test{problem}.txt'
        fno = f'{RESULTS_DIR}/{problem}.txt'

        #load the problems for this set
        problems = loadProblems(fn)

        #generate results for each one
        print(f'\tSolving problem {problem}...')
        all_results = solveProblem(problems)
        
        #finally save the results
        writeResults(fno, all_results)
