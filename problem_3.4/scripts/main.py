
import os
import copy
import json
import numpy as np

PROBLEM_DIR = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
DATA_DIR = f'{PROBLEM_DIR}/data'
RESULTS_DIR = f'{PROBLEM_DIR}/results'

def loadProblems(fn):
    ret = []
    fp = open(fn, 'r')
    
    num_tests = int(fp.readline())
    for x in range(0, num_tests):
        num_people, days = [int(s) for s in fp.readline().rstrip().split(' ')]
        assert(days == 7)
        contact_data = []
        max_person = 1
        for d in range(0, days):
            day_contacts = []
            contact_count = int(fp.readline())
            for c_id in range(0, contact_count):
                c1, c2, p = fp.readline().rstrip().split(' ')
                c1 = int(c1)
                c2 = int(c2)
                p = float(p)
                day_contacts.append((c1, c2, p))

                if c1 > max_person:
                    max_person = c1
                if c2 > max_person:
                    max_person = c2

            contact_data.append(day_contacts)
    
        ret.append({
            'contact_data' : contact_data,
            'comm_size' : max_person
        })

    fp.close()
    return ret

def solveProblem(problem):
    print('Preprocessing...')
    contact_data = problem['contact_data']
    comm_size = problem['comm_size']
    print(f'comm_size: {comm_size}')
    
    #these are all 1-indexed for some reason
    current_infected = [1.0]*(comm_size+1)
    
    '''
    print('days:')
    for d_id, day in enumerate(contact_data):
        print('\t', day)
    '''
    for d_id in range(len(contact_data)-1, -1, -1):
        prev_infected = copy.deepcopy(current_infected)
        day = contact_data[d_id]
        #print('\t', d_id, day)
        for source, dest, prob in day:
            prev_infected[source] += current_infected[dest]*prob

        current_infected = prev_infected

    print('Getting results...')
    am = np.argmax(current_infected)

    return am

def solveProblem_forward(problem):
    print('Preprocessing...')
    contact_data = problem['contact_data']
    comm_size = problem['comm_size']
    print(f'comm_size: {comm_size}')
    
    #these are all 1-indexed for some reason
    current_infected = [1.0]*(comm_size+1)
    
    for d_id in range(len(contact_data)-1, -1, -1):
        prev_infected = copy.deepcopy(current_infected)
        day = contact_data[d_id]
        for source, dest, prob in day:
            prev_infected[source] += current_infected[dest]*prob
        current_infected = prev_infected

    print('Getting results...')
    rev_order = np.argsort(current_infected)

    infected_counts = [0.0]*(comm_size+1)

    #for x in range(1, comm_size+1):
    candidates = 256
    for i in range(0, min(candidates, comm_size)):
        x = rev_order[-i-1]
        
        #one person starts out infected
        infected = [0.0]*(comm_size+1)
        infected[x] = 1.0

        #lets see how much damage they do
        for d_id in range(0, len(contact_data)):
            next_infected = copy.deepcopy(infected)
            day = contact_data[d_id]
            for source, dest, prob in day:
                #probability already infected + (prob NOT infected * prob source is infected * transmission prob)
                next_infected[dest] += (1.0-next_infected[dest])*infected[source]*prob
            infected = next_infected
        
        expected_infections = np.sum(infected)
        #infected_counts.append(expected_infections)
        infected_counts[x] = expected_infections
        print(f'candidate {i} -> {x} = {expected_infections}')

    print('Getting results...')
    am = np.argmax(infected_counts)

    return am
    
def writeResults(fn, all_results):
    fp = open(fn, 'w+')
    for result in all_results:
        fp.write(f'{result}\n')
    fp.close()

if __name__ == '__main__':
    #there are usually multiple per problem
    starting_problem = 5
    ending_problem = 7

    if not os.path.exists(RESULTS_DIR):
        os.makedirs(RESULTS_DIR)

    #go through each ones
    for problem in range(starting_problem, ending_problem+1):
        #filenames below might need to change per problem
        print(f'Analyzing problem set #{problem}...')
        fn = f'{DATA_DIR}/test{problem}'
        fno = f'{RESULTS_DIR}/{problem}.txt'

        #load the problems for this set
        problems = loadProblems(fn)

        #generate results for each one
        all_results = []
        for i, problem in enumerate(problems):
            print(f'\tSolving problem {i}...')
            result = solveProblem_forward(problem)
            all_results.append(result)

        #finally save the results
        writeResults(fno, all_results)
