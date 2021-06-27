
import bisect
import os
import json
import glob

PROBLEM_DIR = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
DATA_DIR = f'{PROBLEM_DIR}/data'
RESULTS_DIR = f'{PROBLEM_DIR}/results'

def loadProblems(fn):
    ret = {}
    fp = open(fn, 'r')
    
    num_isoforms, delta = [int(s) for s in fp.readline().rstrip().split(' ')]
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
    num_isoforms = len(isoforms)

    #initial answers are all -1 and 0
    first_match = [-1]*len(reads)
    num_match = [0]*len(reads)

    #strategy: iterate through each isoform, compare to the read to mark as a match or not; first find goes into find bucket
    print('Preprocessing reads...')
    parsed_reads = [parseRanges(r) for r in reads]

    print('Isoform iterating...')
    for iso_num, isoform in enumerate(isoforms):
        print(f'\tIsoform #{iso_num} / {num_isoforms}...')
        iso_ranges = parseRanges(isoform)
        #print(iso_ranges)
        for read_num, pr in enumerate(parsed_reads):
            #print(pr)
            #compare them here
            pr_start = pr[0]
            ind = bisect.bisect(iso_ranges, pr_start)
            
            initial_check = max(0, ind-1)
            maximal_check = min(len(iso_ranges), ind+1)
            
            initial_range = -1
            for si in range(initial_check, maximal_check):
                curr_range = iso_ranges[si]
                if (pr_start[0] >= curr_range[0]-delta and
                    abs(pr_start[1] - curr_range[1]) <= delta):
                    #this is the initial range
                    initial_range = si
                    break
            
            if initial_range != -1 and (initial_range+len(pr) <= len(iso_ranges)):
                #this is a candidate, check all the in-between regions
                mismatch = False
                for offset in range(1, len(pr)-1):
                    pr_range = pr[offset]
                    curr_range = iso_ranges[initial_range+offset]
                    if (abs(pr_range[0] - curr_range[0]) <= delta and
                        abs(pr_range[1] - curr_range[1]) <= delta):
                        #still good
                        pass
                    else:
                        mismatch = True
                        break

                #TODO: finally check the tail region
                pr_range = pr[-1]
                curr_range = iso_ranges[initial_range+len(pr)-1]
                if (abs(pr_range[0] - curr_range[0]) <= delta and
                    pr_range[1] <= curr_range[1]+delta):
                    pass
                else:
                    mismatch = True

                if not mismatch:
                    num_match[read_num] += 1
                    if first_match[read_num] == -1:
                        first_match[read_num] = iso_num

            #exit()

    print('Getting results...')
    results = {
        'first_match' : first_match,
        'num_match' : num_match
    }
    return results
    
def writeResults(fn, all_results):
    fp = open(fn, 'w+')
    first_match = all_results['first_match']
    num_match = all_results['num_match']
    for i, fm in enumerate(first_match):
        fp.write(f'{fm} {num_match[i]}\n')
    fp.close()

if __name__ == '__main__':
    #there are usually multiple per problem
    all_filenames = sorted(glob.glob(f'{DATA_DIR}/*.txt'))
    starting_problem = 6
    ending_problem = 6

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
