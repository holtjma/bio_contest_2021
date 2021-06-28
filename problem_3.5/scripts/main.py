
import os
import json

PROBLEM_DIR = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
DATA_DIR = f'{PROBLEM_DIR}/data'
RESULTS_DIR = f'{PROBLEM_DIR}/results'

from sklearn.cluster import AgglomerativeClustering
import numpy as np

def loadProblems(fn):
    ret = {}
    fp = open(fn, 'r')
    
    num_test, max_score = [int(s) for s in fp.readline().rstrip().split(' ')]
    beta = float(fp.readline())
    bonus_points = (fp.readline().rstrip() == '+')

    tests = []
    for x in range(0, num_test):
        assert(fp.readline().rstrip() == '---')
        ref_genome = fp.readline().rstrip()
        min_morph, num_reads = [int(s) for s in fp.readline().rstrip().split(' ')]
        hap_rate, error_rate = [float(s) for s in fp.readline().rstrip().split(' ')]
        reads = [fp.readline().rstrip() for r in range(0, num_reads)]
        tests.append({
            'ref' : ref_genome,
            'min_polymorphs' : min_morph,
            'hap_rate' : hap_rate,
            'error_rate' : error_rate,
            'reads' : reads
        })
    
    ret = {
        'max_score' : max_score,
        'beta' : beta,
        'bonus_points' : bonus_points,
        'tests' : tests
    }

    fp.close()
    return ret

def countDelta(s1, s2):
    c = 0
    for i in range(0, len(s1)):
        if s1[i] != s2[i]:
            c += 1
    return c

def generateConsensus(reads, ref):
    s = ''
    for x in range(0, len(reads[0])):
        d = {}
        for i in range(0, len(reads)):
            c = reads[i][x]
            d[c] = d.get(c, 0)+1
        max_value = 0
        max_c = ''
        for c in d:
            if d[c] > max_value:
                max_value = d[c]
                max_c = c
            elif d[c] == max_value and c == ref[x]:
                max_c = c
        s += max_c
    return s

def solveProblem(problem):
    #Load the problem dict (or list of problem dicts)
    #metabolites = problem['metabolites']
    #print(json.dumps(problem))
    #exit()

    fracs = []
    bonus_points = problem['bonus_points']
    print(f'Bonus points? {bonus_points}')
    for i, test in enumerate(problem['tests']):
        ref = test['ref']
        min_polymorphs = test['min_polymorphs']
        hap_rate = test['hap_rate']
        error_rate = test['error_rate']
        reads = test['reads']
        num_reads = len(reads)

        print('', i, min_polymorphs, hap_rate, error_rate, len(reads), sep='\t')
        
        distance_matrix = np.zeros(dtype='<u4', shape=(num_reads, num_reads))
        for i in range(0, num_reads):
            for j in range(i+1, num_reads):
                d = countDelta(reads[i], reads[j])
                distance_matrix[i, j] = d
                distance_matrix[j, i] = d
        
        #single was awful
        #clusterer = AgglomerativeClustering(n_clusters=2, affinity='precomputed', linkage='average')
        clusterer = AgglomerativeClustering(n_clusters=2, affinity='precomputed', linkage='complete')
        labels = clusterer.fit_predict(distance_matrix)
        #print(labels)
        
        lc2 = np.sum(labels)
        lc1 = num_reads - lc2
        
        #assert(lc1 >= lc2)
        
        '''
        if lc1 < lc2:
            l = lc1
            lc1 = lc2
            lc2 = l
        '''
    
        #TODO: check the haplotype differences
        rs1 = []
        rs2 = []
        for x, l in enumerate(labels):
            if l == 0:
                rs1.append(reads[x])
            else:
                rs2.append(reads[x])
            
        con1 = generateConsensus(rs1, ref)
        delta1 = countDelta(con1, ref)
        rs1_count = len(rs1)
        
        con2 = generateConsensus(rs2, ref)
        delta2 = countDelta(con2, ref)
        rs2_count = len(rs2)

        '''
        print('\td1', delta1, len(rs1))
        print('\td2', delta2, len(rs2))
        #'''
        #assert(delta2 > delta1)
        if delta2 < delta1:
            delta2 = delta1
            con2 = con1
            rs2_count = rs1_count

        if delta2 >= min_polymorphs:
            #probability time?
            #first lets do the chance of seq error against ref
            errors = [countDelta(r, ref) for r in reads]
            total_errors = np.sum(errors)
            total_len = len(ref)*num_reads
            #print((min(len(rs1), len(rs2)) / num_reads))
            frac_val = (min(len(rs1), len(rs2)) / num_reads)
            #error_morph_ratio = total_errors / min_polymorphs
            error_morph_ratio = abs(total_errors - num_reads*delta2) / num_reads
            fracs.append((frac_val, error_morph_ratio))
            print('\ttotal_error', total_errors, total_len, total_errors/total_len, frac_val)
            
            #change this for each test
            #if frac_val < 0.01:
            if frac_val < 0.07 and error_morph_ratio > 2.75:# and total_errors < 1000:
                yield ("1", None)
            else:
                yield("2", con2)
        else:
            yield ("1", None)
        print()

    sf = sorted(fracs)
    print('fracs')#, sorted(fracs))
    for f in sf:
        print(f)
        
def writeResults(fn, all_results):
    fp = open(fn, 'w+')
    for result in all_results:
        if result[0] == '1':
            fp.write('1\n')
        else:
            fp.write(f'2 {result[1]}\n')
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
        fn = f'{DATA_DIR}/in{problem}.txt'
        fno = f'{RESULTS_DIR}/{problem}.txt'

        #load the problems for this set
        problems = loadProblems(fn)

        #generate results for each one
        all_results = solveProblem(problems)
        
        #finally save the results
        writeResults(fno, all_results)
