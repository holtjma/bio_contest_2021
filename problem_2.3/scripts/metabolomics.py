
import bisect
import json

def loadProblems(fn):
    ret = []
    fp = open(fn, 'r')
    num_problems = int(fp.readline())
    for x in range(0, num_problems):
        fields = fp.readline().strip().split(' ')
        num_meta = int(fields[0])
        num_adducts = int(fields[1])
        num_obs = int(fields[2])
        
        metabolites = [
            int(m.replace('.', '')) for m in fp.readline().rstrip().split(' ')
        ]
        assert(len(metabolites) == num_meta)

        adducts = [
            int(a.replace('.', '')) for a in fp.readline().rstrip().split(' ')
        ]
        assert(len(adducts) == num_adducts)

        obs = [
            int(o.replace('.', '')) for o in fp.readline().rstrip().split(' ')
        ]
        assert(len(obs) == num_obs)

        ret.append({
            'metabolites' : metabolites,
            'adducts' : adducts,
            'observations' : obs
        })

    fp.close()
    return ret

def solveProblem_old(problem):
    #print(json.dumps(problem, indent=4))
    
    metabolites = problem['metabolites']
    adducts = problem['adducts']
    observations = problem['observations']

    print('Calculating sums...')
    combos = []
    for m in range(0, len(metabolites)):
        if m % 1000 == 0:
            print(f'\t{m}...')
        for k in range(0, len(adducts)):
            d = metabolites[m]+adducts[k]
            combos.append((d, m, k))
    print('Sorting...')
    combos.sort()

    print('Getting results...')
    results = []
    for i, obs in enumerate(observations):
        if i % 1000 == 0:
            print(f'\t{i}...')
        
        best_delta = 0xFFFFFFFFFFFFFFFF
        best_m = 1
        best_k = 1
        insert_point = bisect.bisect(combos, (obs, 0, 0))
        
        for x in range(
            max(0, insert_point-1), min(insert_point+2, len(combos))
        ):
            delta = abs(obs - combos[x][0])
            if delta < best_delta:
                best_delta = delta
                best_m = combos[x][1]+1
                best_k = combos[x][2]+1
        
        results.append((best_m, best_k))

    return results

def solveProblem(problem):
    #print(json.dumps(problem, indent=4))
    
    metabolites = problem['metabolites']
    adducts = problem['adducts']
    observations = problem['observations']

    met_sort = []
    for i, m in enumerate(metabolites):
        met_sort.append((m, i))
    
    print('Sorting metabolites...')
    met_sort.sort()

    print('Running observations...')
    results = []
    for i, obs in enumerate(observations):
        if i % 100 == 0:
            print(f'\t{i}...')
        
        best_delta = 0xFFFFFFFFFFFFFFFF
        best_m = 1
        best_k = 1

        for k, adduct in enumerate(adducts):
            target = obs - adduct
            insert_point = bisect.bisect(met_sort, (target, 0))
            
            for x in range(
                max(0, insert_point-1), min(insert_point+2, len(met_sort))
            ):
                delta = abs(met_sort[x][0] - target)
                if delta < best_delta:
                    best_delta = delta
                    best_m = met_sort[x][1]+1
                    best_k = k+1
        
        results.append((best_m, best_k))

    return results
    
def writeResults(fn, all_results):
    fp = open(fn, 'w+')
    for result in all_results:
        for tup in result:
            fp.write(f'{tup[0]} {tup[1]}\n')
    fp.close()

if __name__ == '__main__':
    for problem in range(5, 6):
        print(f'Analyzing problem set #{problem}...')
        fn = f'/Users/matt/githubProjects/bio_contest_quals_2021/problem_2.3/data/{problem}.txt'
        fno = f'/Users/matt/githubProjects/bio_contest_quals_2021/problem_2.3/results/{problem}.txt'

        problems = loadProblems(fn)
        all_results = []
        for i, problem in enumerate(problems):
            print(f'\tSolving problem {i}...')
            result = solveProblem(problem)
            #print(json.dumps(result, indent=4))
            all_results.append(result)

        writeResults(fno, all_results)
