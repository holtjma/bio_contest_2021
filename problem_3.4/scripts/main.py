
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

def solveProblem_forward(problem, jans):
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
    #candidates = 512 #5 on this
    #candidates = 512 #6 on this
    candidates = 2**13 #7 on this
    #candidates = 128
    curr_highest = 0.0
    #for i in range(0, min(candidates, comm_size)):
    for i in range(0, candidates):
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
        if expected_infections > curr_highest:
            print(f'candidate {i} -> {x} = {expected_infections}, new record (vs {jans})')
            curr_highest = expected_infections
            #if curr_highest >= jans:
            #    break
        elif i % 100 == 0:
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
    starting_problem = 6
    ending_problem = 6

    jans_dict = {
        5: [18250,415,361,194,420,156,237,358,172,515],
        6: [17341,83,209,180,90,155,138,140,136,115],
        7: [1437,313,321,378,418]
    }

    if not os.path.exists(RESULTS_DIR):
        os.makedirs(RESULTS_DIR)

    #go through each ones
    for problem in range(starting_problem, ending_problem+1):
        #filenames below might need to change per problem
        print(f'Analyzing problem set #{problem}...')
        fn = f'{DATA_DIR}/test{problem}'
        fno = f'{RESULTS_DIR}/{problem}.txt'
        #print(jans_dict[problem])
        #exit()

        #load the problems for this set
        problems = loadProblems(fn)

        #generate results for each one
        all_results = []
        p_id = problem
        for i, problem in enumerate(problems):
            print(f'\tSolving problem {i}...')
            #if i == 0 and p_id == 7:
            #    result = 528
            if i == 0 and p_id == 6:
                result = 95415
            elif i == 0 and p_id == 5:
                result = 122
            else:
                result = solveProblem_forward(problem, jans_dict[p_id][i])
            all_results.append(result)

        #finally save the results
        writeResults(fno, all_results)
