
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

def makeIntrons(reg):
    introns = []
    for x in range(0, len(reg)-1):
        introns.append((reg[x][1], reg[x+1][0]))
    return introns

def overlapRanges(i, r):
    r_start = r[0]
    curr_ind = max(0, bisect.bisect_left(i, r_start)-1)
    i_s, i_e = i[curr_ind]
    overlap = 0

    r_ind = 0
    i_ind = curr_ind

    while r_ind < len(r) and i_ind < len(i):
        r_s, r_e = r[r_ind]
        i_s, i_e = i[i_ind]
        
        max_start = max(i_s, r_s)
        min_end = min(i_e, r_e)
        if max_start < min_end:
            overlap += min_end - max_start
            
        if r_e < i_e:
            r_ind += 1
        else:
            i_ind += 1

    return overlap

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
    iso_introns = [
        makeIntrons(iso) for iso in iso_ranges
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
        read_introns = makeIntrons(parsed_read)

        best_match = 0
        best_score = 0

        start = parsed_read[0][0]
        end = parsed_read[-1][1]
        print(f'\t[{datetime.datetime.now()}] Read #{rc} / {num_reads}, l={end-start}...')
        pr_arr = np.zeros(dtype='bool', shape=(end - start, ))
        
        exon_bp = np.sum([a[1]-a[0] for a in parsed_read])
        intron_bp = (end - start) - exon_bp
        
        #TODO: make this more efficiently search, for now just do all of them
        search_space = range(0, len(iso_ranges))
        
        #now search the search space
        #for iso_ind, iso_range in enumerate(iso_ranges):
        for iso_ind in search_space:
            iso_range = iso_ranges[iso_ind]
            iso_start = iso_range[0][0]
            iso_end = iso_range[-1][1]
            
            if iso_start > end or iso_end < start:
                #score = 0
                pass
            else:
                exon_count = overlapRanges(iso_range, parsed_read)
                intron_count = overlapRanges(iso_introns[iso_ind], read_introns)
                score = 2*(exon_count / exon_bp) + (intron_count / intron_bp)
                
                if score > best_score:
                    best_score = score
                    best_match = iso_ind
                    #print(exon_count, exon_bp, intron_count, intron_bp)

        #print(best_score / 3.0)
        print(best_match, best_score)
        #exit()
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
    starting_problem = 7
    ending_problem = 7
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
