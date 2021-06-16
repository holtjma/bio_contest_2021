
import json

def loadProblems(fn):
    ret = []
    fp = open(fn, 'r')
    num_problems = int(fp.readline())
    for x in range(0, num_problems):
        fields = fp.readline().strip().split(' ')
        num_entries = int(fields[0])
        entry_len = int(fields[1])
        entries = []
        for y in range(0, num_entries):
            entries.append(fp.readline().rstrip())
        
        ret.append({
            'num_entries' : num_entries,
            'entry_len' : entry_len,
            'entries' : entries
        })

    fp.close()
    return ret

def solveProblem(problem):
    #print(json.dumps(problem, indent=4))
    
    num_entries = problem['num_entries']
    entry_len = problem['entry_len']
    entries = problem['entries']

    modes = {}
    next_id = 1

    result = []
    for x in range(0, entry_len):
        join_mode = ''.join([entries[y][x] for y in range(0, num_entries)])
        if join_mode not in modes:
            modes[join_mode] = next_id
            next_id += 1
        result.append(modes[join_mode])

    return {
        'num_modes' : next_id - 1,
        'ordered' : result
    }
    
def writeResults(fn, all_results):
    fp = open(fn, 'w+')
    for result in all_results:
        num_modes = result['num_modes']
        modes = result['ordered']
        fp.write(f'{num_modes}\n')
        fp.write(' '.join([str(s) for s in modes])+'\n')
    fp.close()

if __name__ == '__main__':
    #problem = '1'
    problem = '2'
    fn = f'/Users/matt/githubProjects/bio_contest_quals_2021/problem_2.2/data/{problem}.txt'
    fno = f'/Users/matt/githubProjects/bio_contest_quals_2021/problem_2.2/results/{problem}.txt'

    problems = loadProblems(fn)
    all_results = []
    for i, problem in enumerate(problems):
        print(f'Solving problem {i}...')
        result = solveProblem(problem)
        #print(json.dumps(result, indent=4))
        all_results.append(result)

    writeResults(fno, all_results)
