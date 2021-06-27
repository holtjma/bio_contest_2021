
import bisect
import datetime
import os
import json
import glob
import numpy as np

PROBLEM_DIR = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
DATA_DIR = f'{PROBLEM_DIR}/data'
RESULTS_DIR = f'{PROBLEM_DIR}/results'

def loadProblems(fn):
    ret = {}
    fp = open(fn, 'r')
    
    num_isoforms, delta = [int(s) for s in fp.readline().rstrip().split(' ')]
    #assert(delta == 0)
    isoforms = []
    for x in range(0, num_isoforms):
        isoforms.append(fp.readline().rstrip())
    
    num_reads = int(fp.readline())
    reads = []
    for x in range(0, num_reads):
        reads.append(fp.readline().rstrip())

    fp.close()

    ret['isoforms'] = isoforms
    ret['reads'] = reads
    ret['delta'] = delta

    return ret

def getMaxEnd(data):
    m = max(int(d.split('-')[-1]) for d in data)
    return m

def getMaxLength(data):
    m = np.max(np.array(
        [int(d.split('-')[-1]) - int(d.split('-')[0]) for d in data]
    ))
    return m

def parseRanges(s):
    regions = s.split(',')
    ret = []
    for r in regions:
        s, e = r.split('-')
        ret.append((int(s), int(e)))
    return ret

def solveProblem(problem):
    #Load the problem dict (or list of problem dicts)
    isoforms = problem['isoforms']
    reads = problem['reads']
    delta = problem['delta']

    #figure out just how big our sparse matrix needs to be
    max_end = max(getMaxEnd(isoforms), getMaxEnd(reads))+1
    max_length = getMaxLength(isoforms)
    max_length2 = getMaxLength(reads)
    print(len(isoforms), len(reads), max_end, max_length, max_length2)

    #different strategy, we are now looping over reads and need to pre-process isoforms
    print('Preprocessing isoforms...')
    iso_ranges = [
        parseRanges(iso) for iso in isoforms
    ]

    print('Getting mins/maxs...')
    min_starts = np.array([ir[0][0] for ir in iso_ranges])
    max_ends = np.array([ir[-1][1] for ir in iso_ranges])

    print('Sorting isoforms...')
    start_inds = np.argsort(min_starts)
    start_ordered = min_starts[start_inds]
    #end_inds = np.argsort(max_ends)
    #end_ordered = max_ends[end_inds]

    print('Processing reads...')
    num_reads = len(reads)
    for rc, r in enumerate(reads):
        parsed_read = parseRanges(r)

        best_match = 0
        best_score = 0

        start = parsed_read[0][0]
        end = parsed_read[-1][1]
        print(f'\t[{datetime.datetime.now()}] Read #{rc} / {num_reads}, l={end-start}...')
        pr_arr = np.zeros(dtype='bool', shape=(end - start, ))
        for s, e in parsed_read:
            pr_arr[s-start:e-start] = True
        intronic_arr = np.logical_not(pr_arr)

        exon_bp = np.sum(pr_arr)
        intron_bp = np.sum(intronic_arr)
        
        #TODO: make this more efficiently search, for now just do all of them
        #search_space = range(0, len(iso_ranges))
        #'''
        end_search = bisect.bisect_left(start_ordered, end)
        assert(end_search == 0 or start_ordered[end_search-1] < end)
        search_space = start_inds
        #'''

        #now search the search space
        #for iso_ind, iso_range in enumerate(iso_ranges):
        for iso_ind in search_space:
            iso_range = iso_ranges[iso_ind]
            iso_start = iso_range[0][0]
            iso_end = iso_range[-1][1]
            #if iso_start > end or iso_end < start:
            #    score = 0
            if iso_start > end:
                break
            elif iso_end < start:
                #score = 0
                pass
            else:
                iso_arr = np.zeros(dtype='bool', shape=(end - start, ))
                start_ind = max(0, bisect.bisect(iso_range, (start, ))-2)
                for s, e in iso_range[start_ind:]:
                    if e < start:
                        pass
                    elif s > end:
                        break
                    else:
                        smod = max(0, s - start)
                        emod = min(end-start, e - start)
                        iso_arr[smod:emod] = True
                
                exon_count = np.sum(pr_arr & iso_arr)
                intron_count = np.sum(intronic_arr & np.logical_not(iso_arr))
                score = 2*(exon_count / exon_bp) + (intron_count / intron_bp)
                if score > best_score:
                    best_score = score
                    best_match = iso_ind

        #print(best_score / 3.0)
        yield best_match
    
def writeResults(fn, all_results):
    fp = open(fn, 'w+')
    for bm in all_results:
        fp.write(f'{bm}\n')
    fp.close()

if __name__ == '__main__':
    #there are usually multiple per problem
    all_filenames = sorted(glob.glob(f'{DATA_DIR}/*.txt'))
    #starting_problem = 0
    #ending_problem = 0

    #only for these
    starting_problem = 8
    ending_problem = 8
    #ending_problem = 9

    if not os.path.exists(RESULTS_DIR):
        os.makedirs(RESULTS_DIR)

    #go through each ones
    for problem in range(starting_problem, ending_problem+1):
        #filenames below might need to change per problem
        print(f'Analyzing problem set #{problem}...')
        #fn = f'{DATA_DIR}/{problem}.txt'
        fn = all_filenames[problem]
        fno = f'{RESULTS_DIR}/{problem}.txt'

        #load the problems for this set
        problems = loadProblems(fn)

        #generate results for each one
        all_results = solveProblem(problems)

        #finally save the results
        writeResults(fno, all_results)
