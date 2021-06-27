
import os
import json
import matplotlib.pyplot as plt
import numpy as np

PROBLEM_DIR = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
DATA_DIR = f'{PROBLEM_DIR}/data'
RESULTS_DIR = f'{PROBLEM_DIR}/results'

def loadProblems(fn):
    fp = open(fn, 'r')
    
    num_tests = int(fp.readline())
    tests = []
    for x in range(0, num_tests):
        num_org, num_blocks = [int(s) for s in fp.readline().rstrip().split(' ')]
        affected = []
        unaffected = []
        for y in range(0, num_org):
            dtype = fp.readline().rstrip()
            if dtype == '+':
                affected.append((
                    fp.readline().rstrip(),
                    fp.readline().rstrip()
                ))
            elif dtype == '-':
                unaffected.append((
                    fp.readline().rstrip(),
                    fp.readline().rstrip()
                ))
            else:
                raise Exception('')

        tests.append({
            'affected' : affected,
            'unaffected' : unaffected
        })

    fp.close()
    return tests

convert = {
    'A' : 0,
    'C' : 1,
    'G' : 2,
    'T' : 3
}

def convertKmer(kmer):
    i = 0
    for c in kmer:
        i *= 4
        i += convert[c]
    return i

def solveProblem_rf2(problem):
    #Load the problem dict (or list of problem dicts)
    #metabolites = problem['metabolites']
    affected = problem['affected']
    unaffected = problem['unaffected']
    block_len = len(affected[0][0])

    from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
    from sklearn.metrics import accuracy_score

    k = 5
    accs = []
    multi_accs = []
    multi_imps = []
    for r in range(0, k-1):
        multi_accs.append([])
        multi_imps.append([])

    for i in range(0, block_len-k+1):
        xs = []
        ys = []
        multi_accs.append([])
        multi_imps.append([])
        for aff in affected:
            xarr = []
            for j in range(i, i+k):
                xarr.append(convertKmer(aff[0][j]))
                xarr.append(convertKmer(aff[1][j]))
            xs.append(xarr)
            ys.append(1)
        for unaff in unaffected:
            xarr = []
            for j in range(i, i+k):
                xarr.append(convertKmer(unaff[0][j]))
                xarr.append(convertKmer(unaff[1][j]))
            xs.append(xarr)
            ys.append(0)
    
        forest = GradientBoostingClassifier(random_state=0)
        forest.fit(xs, ys)
        
        pred = forest.predict(xs)
        importances = forest.feature_importances_

        acc = accuracy_score(ys, pred)
        print(i, acc, importances)
        accs.append(acc)
        for offset in range(0, k):
            multi_accs[i+offset].append(acc)
            multi_imps[i+offset].append(importances[2*offset])
            multi_imps[i+offset].append(importances[2*offset+1])

    multi_accs = np.array([np.mean(s) for s in multi_accs])
    multi_imps = np.array([np.mean(s) for s in multi_imps])

    rs_max = np.argmax(multi_accs)
    '''
    plt.figure()
    plt.plot(range(0 + k // 2, block_len + k // 2 - k+1), accs, label='accs')
    plt.plot(multi_accs, label='multi')
    plt.plot(multi_imps, label='imps')
    plt.xlim(rs_max-10, rs_max+10)
    plt.axvline(rs_max)
    plt.grid()
    plt.legend()
    plt.show()
    #'''

    #for some stupid reason, you need to add a +1
    buffer = 0
    return (rs_max-buffer, rs_max+buffer+1)

def writeResults(fn, all_results):
    fp = open(fn, 'w+')
    print(fn)
    for result in all_results:
        print(result)
        fp.write(f'{result[0]} {result[1]}\n')
    fp.close()

if __name__ == '__main__':
    #there are usually multiple per problem
    starting_problem = 7
    ending_problem = 7

    if not os.path.exists(RESULTS_DIR):
        os.makedirs(RESULTS_DIR)

    #go through each ones
    for problem in range(starting_problem, ending_problem+1):
        #filenames below might need to change per problem
        print(f'Analyzing problem set #{problem}...')
        fn = f'{DATA_DIR}/0{problem}.txt'
        fno = f'{RESULTS_DIR}/{problem}.txt'

        #load the problems for this set
        problems = loadProblems(fn)

        #generate results for each one
        all_results = []
        for i, problem in enumerate(problems):
            print(f'\tSolving problem {i}...')
            #if i != 8:
            #    continue
            #result = solveProblem(problem)
            result = solveProblem_rf2(problem)
            all_results.append(result)

        #finally save the results
        writeResults(fno, all_results)
