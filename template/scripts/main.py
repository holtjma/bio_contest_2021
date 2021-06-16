
import os
import json

PROBLEM_DIR = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
DATA_DIR = f'{PROBLEM_DIR}/data'
RESULTS_DIR = f'{PROBLEM_DIR}/results'

def loadProblems(fn):
    ret = {}
    fp = open(fn, 'r')
    
    #TODO: fill in loading script, return is either a dictionary or a list of dictionaries typically

    fp.close()
    return ret

def solveProblem(problem):
    #Load the problem dict (or list of problem dicts)
    #metabolites = problem['metabolites']
    
    print('Preprocessing...')
    #TODO: any data structure setup

    print('Getting results...')
    results = []
    
    #TODO: calculate results

    return results
    
def writeResults(fn, all_results):
    fp = open(fn, 'w+')
    for result in all_results:
        #TODO: something problem dependent
        pass
    fp.close()

if __name__ == '__main__':
    #there are usually multiple per problem
    starting_problem = 1
    ending_problem = 1

    if not os.path.exists(RESULTS_DIR):
        os.makedirs(RESULTS_DIR)

    #go through each ones
    for problem in range(starting_problem, ending_problem+1):
        #filenames below might need to change per problem
        print(f'Analyzing problem set #{problem}...')
        fn = f'{DATA_DIR}/{problem}.txt'
        fno = f'{RESULTS_DIR}/{problem}.txt'

        #load the problems for this set
        problems = loadProblems(fn)

        #generate results for each one
        all_results = []
        for i, problem in enumerate(problems):
            print(f'\tSolving problem {i}...')
            result = solveProblem(problem)
            all_results.append(result)

        #finally save the results
        writeResults(fno, all_results)
