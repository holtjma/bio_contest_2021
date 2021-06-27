
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
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
                affected.append(fp.readline().rstrip())
            elif dtype == '-':
                unaffected.append(fp.readline().rstrip())
            else:
                raise Exception('')

        tests.append({
            'affected' : affected,
            'unaffected' : unaffected
        })

    fp.close()
    return tests

CHARS = ['A', 'C', 'G', 'T']
def solveProblem(problem):
    #Load the problem dict (or list of problem dicts)
    #metabolites = problem['metabolites']
    affected = problem['affected']
    unaffected = problem['unaffected']
    block_len = len(affected[0])

    arrays = {
        'A' : [0]*block_len,
        'C' : [0]*block_len,
        'G' : [0]*block_len,
        'T' : [0]*block_len
    }

    for aff in affected:
        for i, c in enumerate(aff):
            arrays[c][i] += 1
    
    unarrays = {
        'A' : [0]*block_len,
        'C' : [0]*block_len,
        'G' : [0]*block_len,
        'T' : [0]*block_len
    }
    for unaff in unaffected:
        for i, c in enumerate(unaff):
            unarrays[c][i] += 1
    
    raw_scoring = []
    for i in range(0, block_len):
        best_score = 0
        bestc = 'A'
        bestc2 = 'C'
        for c in CHARS:
            for c2 in CHARS:
                if c == c2:
                    pass
                else:
                    score = arrays[c][i] + unarrays[c2][i]
                    if score > best_score:
                        best_score = score
                        bestc = c
                        bestc2 = c2
        
        raw_scoring.append((best_score, i, bestc, bestc2))
    
    raw_scoring.sort()
    for rs in raw_scoring:
        print(rs)

    rs_max = raw_scoring[-1][1]

    print('affected')
    for aff in affected:
        print('\t', aff[rs_max])
    print('unaffected')
    for unaff in unaffected:
        print('\t', unaff[rs_max])

    
    #'''
    plt.figure()
    for c in arrays:
        plt.plot(arrays[c], label=c)
    for c in unarrays:
        plt.plot(0.0-np.array(unarrays[c]), label=f'un-{c}')
    plt.grid()
    plt.legend()
    plt.xlim(rs_max-3, rs_max+3)
    plt.show()
    #'''
    
    return (rs_max, rs_max)

def solveProblem2(problem):
    #Load the problem dict (or list of problem dicts)
    #metabolites = problem['metabolites']
    affected = problem['affected']
    unaffected = problem['unaffected']
    block_len = len(affected[0])

    arrays = {}

    k = 5
    for aff in affected:
        for i in range(0, block_len-k):
            kmer = aff[i:i+k]
            if kmer not in arrays:
                arrays[kmer] = [0]*block_len
            arrays[kmer][i] += 1
    
    unarrays = {}
    for unaff in unaffected:
        for i in range(0, block_len-k):
            kmer = unaff[i:i+k]
            if kmer not in unarrays:
                unarrays[kmer] = [0]*block_len
            unarrays[kmer][i] += 1
    
    raw_scoring = []
    for i in range(0, block_len-k):
        best_score = 0
        s1 = 0
        bestc = ''
        s2 = 0
        bestc2 = ''
        for kmer in arrays:
            if arrays[kmer][i] > s1:
                bestc = kmer
                s1 = arrays[kmer][i]
        
        for unkmer in unarrays:
            if unarrays[unkmer][i] > s2:
                bestc2 = unkmer
                s2 = unarrays[kmer][i]

        assert(kmer != unkmer)
        
        score = s1+s2
        if score > best_score:
            best_score = score
            bestc = kmer
            bestc2 = unkmer
        
        raw_scoring.append((best_score, i, bestc, bestc2))
    raw_scoring.sort()
    for rs in raw_scoring:
        print(rs)
    rs_max = raw_scoring[-1][1]

    #'''
    plt.figure()
    for c in arrays:
        plt.plot(arrays[c], label=c)
    for c in unarrays:
        plt.plot(0.0-np.array(unarrays[c]), label=f'un-{c}')
    plt.grid()
    #plt.legend()
    plt.xlim(rs_max-3, rs_max+3)
    plt.show()
    #'''
    
    return (rs_max, rs_max)

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

def solveProblem_rf(problem):
    #Load the problem dict (or list of problem dicts)
    #metabolites = problem['metabolites']
    affected = problem['affected']
    unaffected = problem['unaffected']
    block_len = len(affected[0])

    x = []
    y = []
    k = 5
    for aff in affected:
        arr = []
        for i in range(0, len(aff)-k+1):
            arr.append(convertKmer(aff[i:i+k]))
        x.append(arr)
        y.append(1)
    
    for unaff in unaffected:
        arr = []
        for i in range(0, len(unaff)-k+1):
            arr.append(convertKmer(unaff[i:i+k]))
        x.append(arr)
        y.append(0)
    
    #works when k = 1
    #assert(len(x[0]) == len(affected[0]))

    from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier

    forest = GradientBoostingClassifier(random_state=0)
    forest.fit(x, y)

    importances = forest.feature_importances_
    base_importances = [0]*len(aff)
    for ind, i_val in enumerate(importances):
        for z in range(0, k):
            base_importances[ind+z] += i_val

    #rs_max = np.argmax(importances)
    rs_max = np.argmax(base_importances)
    buffer = 0
    '''
    plt.figure()
    plt.plot(importances)
    plt.plot(base_importances)
    plt.xlim(rs_max-10, rs_max+10)
    #plt.xlim(0, 10)
    plt.grid()
    plt.show()
    #'''
    return (rs_max-buffer, rs_max+buffer)

def solveProblem_rf2(problem):
    #Load the problem dict (or list of problem dicts)
    #metabolites = problem['metabolites']
    affected = problem['affected']
    unaffected = problem['unaffected']
    block_len = len(affected[0])

    from sklearn.ensemble import RandomForestClassifier,GradientBoostingClassifier
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
            xs.append([convertKmer(c) for c in aff[i:i+k]])
            ys.append(1)
        for unaff in unaffected:
            xs.append([convertKmer(c) for c in unaff[i:i+k]])
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
            multi_imps[i+offset].append(importances[offset])

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
    starting_problem = 4
    ending_problem = 4

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
