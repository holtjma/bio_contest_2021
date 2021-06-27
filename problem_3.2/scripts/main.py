
import os
import numpy as np
import json

PROBLEM_DIR = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
DATA_DIR = f'{PROBLEM_DIR}/data'
RESULTS_DIR = f'{PROBLEM_DIR}/results'

def loadProblems(fn):
    ret = {}
    fp = open(fn, 'r')
    
    num_known, num_partial = [int(x) for x in fp.readline().rstrip().split(' ')]
    
    known_haps = set([])
    for x in range(0, num_known):
        fp.readline()
        new_hap = fp.readline().rstrip()
        new_hap2 = fp.readline().rstrip()
        #I don't think this matters
        #assert(new_hap not in known_haps)
        #assert(new_hap2 not in known_haps)
        known_haps.add(new_hap)
        known_haps.add(new_hap2)

    unknown_haps = []
    for x in range(0, num_partial):
        fp.readline()
        unknown_haps.append(fp.readline().rstrip())

    fp.close()

    ret['known_haps'] = list(known_haps)
    ret['unknown_haps'] = unknown_haps

    return ret

def solveProblem(problem):
    #Load the problem dict (or list of problem dicts)
    #metabolites = problem['metabolites']
    known_haps = problem['known_haps']
    unknown_haps = problem['unknown_haps']

    hap_len = len(known_haps[0])
    print(f'hap_len: {hap_len}')
    print(f'num_known: {len(known_haps)}')
    print(f'num_unknown: {len(unknown_haps)}')

    print('Preprocessing...')
    #figure out where our unknowns are
    known_indices = []
    for i, c in enumerate(unknown_haps[0]):
        if c != '?':
            known_indices.append(i)
    known_indices = np.array(known_indices)
    print(f'len(known_indices): {known_indices.shape[0]}')

    #truncate the search space to the small haplotypes
    hap_arr = []
    small_arr = []
    for kh in known_haps:
        new_array = np.array([int(c) for c in kh], dtype='<i1')
        hap_arr.append(new_array)
        small_arr.append(new_array[known_indices])

    print('Getting results...')
    results = []
    
    for i, uh in enumerate(unknown_haps):
        print(f'\tSolving {i}')
        #print(uh)
        uh_arr = np.array([int(c) for c in uh.replace('?', '')], dtype='<i1')
        #print(uh_arr)
        best_x = []
        best_y = []
        #best_delta = 2*uh_arr.shape[0]
        best_score = 0

        for x in range(0, len(small_arr)):
            for y in range(x, len(small_arr)):
                total = small_arr[x]+small_arr[y]
                #delta = np.sum(np.absolute(total - uh_arr))
                delta_arr = np.absolute(total - uh_arr)
                consec = np.zeros(dtype='<u1', shape=delta_arr.shape)
                consec[0:-1] += (delta_arr[0:-1] == 0) & (delta_arr[1:] == 0)
                consec[1:] += consec[0:-1]

                #print(small_arr[x], small_arr[y], delta_arr, consec)
                score = np.sum((2-delta_arr) + 0.5*consec)
                    

                #if delta < best_delta:
                if score > best_score:
                    best_x = x
                    best_y = y
                    best_score = score
                '''
                if delta < best_delta:
                    best_x = [x]
                    best_y = [y]
                    best_delta = delta
                elif delta == best_delta:
                    best_x.append(x)
                    best_y.append(y)
                '''

                '''
                    print('\t', x, y, small_arr[x], small_arr[y])
                    print('\t', total)
                    print('\t', delta)
                elif delta == best_delta:
                    print('TIE')
                    print('\t', x, y, small_arr[x], small_arr[y])
                    print('\t', total)
                    print('\t', delta)
                '''
                    
        #print('best', best_x, best_y, known_haps[x], known_haps[y])
        
        #print('best_delta', best_delta, len(best_x))
        #assert(len(best_x) == len(best_y))
        '''
        if len(best_x) > 1:
            #tie-breaking
            best_found = (best_x[0], best_y[0])
            best_consec = 0

            for ind in range(0, len(best_x)):
                x = best_x[ind]
                y = best_y[ind]

                total = small_arr[x]+small_arr[y]
                delta_arr = np.absolute(total - uh_arr)

                curr_index = known_indices[0]-1
                curr_inc = 0
                consec_score = 0
                for j in range(0, known_indices.shape[0]):
                    if known_indices[j] == curr_index+1 and delta_arr[j] == 0:
                        #it is consecutive AND equal
                        curr_inc += 1
                        curr_index += 1
                        consec_score += curr_inc
                    elif known_indices[j] == curr_index+1:
                        #it is consecutive but NOT equal
                        curr_index += 1
                        curr_inc = 0
                    elif delta_arr[j] == 0:
                        #it is not consecutive but IS equal
                        curr_inc = 1
                        curr_index = known_indices[j]
                        consec_score += curr_inc
                    else:
                        #it is neither consecutive nor equal
                        curr_inc = 0
                        curr_index = known_indices[j]
                
                if consec_score > best_consec:
                    best_found = (x, y)
                    best_consec = consec_score
            
            best_x, best_y = best_found

        else:
            best_x = best_x[0]
            best_y = best_y[0]
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
