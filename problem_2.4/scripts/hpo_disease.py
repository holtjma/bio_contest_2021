
import bisect
import copy
import heapq
import json
import datetime
import numpy as np

from scipy.sparse import csr_matrix, lil_matrix

#from custom_segment_tree import CustomSegmentTree
from sparse_table import STARunner

def loadProblems(fn):
    fp = open(fn, 'r')
    
    #tree desc
    num_vertices = int(fp.readline())
    parent_ids = [int(x) for x in fp.readline().rstrip().split(' ')]
    assert(len(parent_ids) == num_vertices-1)

    #get IC
    ic = [int(x) for x in fp.readline().rstrip().split(' ')]
    assert(len(ic) == num_vertices)

    #pull out all the diseases
    diseases = []
    num_diseases = int(fp.readline())
    for x in range(0, num_diseases):
        disease_desc = [int(x) for x in fp.readline().rstrip().split(' ')]
        assert(disease_desc[0] == len(disease_desc)-1)
        diseases.append(disease_desc[1:])
    
    #now the patients
    num_patients = int(fp.readline())
    patients = []
    for x in range(0, num_patients):
        patient_desc = [int(x) for x in fp.readline().rstrip().split(' ')]
        assert(patient_desc[0] == len(patient_desc)-1)
        patients.append(patient_desc[1:])
    fp.close()
    
    ret = {
        'num_vertices' : num_vertices,
        'parent_ids' : parent_ids,
        'ic' : ic,
        'diseases' : diseases,
        'patients' : patients
    }
    return ret

def solveProblem(problem):
    #print(json.dumps(problem, indent=4))
    results = []

    #patients & diseases
    patients = problem['patients']
    diseases = problem['diseases']
    
    #graph stuff
    parent_ids = problem['parent_ids']
    ic = problem['ic']

    print('Test stats:')
    print(f'\tic: {len(ic)}')
    print(f'\tPatients: {len(patients)}')
    print(f'\tDiseases: {len(diseases)}')
    
    return

    print('Pre-processing diseases...')
    disease_sets = {}
    for d_id, disease in enumerate(diseases):
        if d_id % 1000 == 0:
            print(f'\t{d_id}...')
        disease_set = set([1])
        for d_pheno in disease:
            current = d_pheno
            while current != 1:
                disease_set.add(current)
                current = parent_ids[current-2]
        disease_sets[d_id] = disease_set

    print('Processing...')
    for p_id, patient in enumerate(patients):
        if True or p_id % 100 == 0:
            print(f'\t{p_id}...')
        best_disease = 0
        best_score = 0

        for d_id, disease in enumerate(diseases):
            print(f'\t\t{d_id}...')
            disease_score = 0
            for p_pheno in patient:
                '''
                disease_set = set([1])
                for d_pheno in disease:
                    current = d_pheno
                    while current != 1:
                        disease_set.add(current)
                        current = parent_ids[current-2]
                '''
                disease_set = disease_sets[d_id]
                current = p_pheno
                while current not in disease_set:
                    current = parent_ids[current-2]
                
                disease_score += ic[current-1]
            
            if disease_score > best_score:
                best_score = disease_score
                best_disease = d_id+1
        
        results.append(best_disease)
    
    return results

def solveProblem_st(problem):
    #print(json.dumps(problem, indent=4))
    results = []

    #patients & diseases
    patients = problem['patients']
    diseases = problem['diseases']
    
    #graph stuff
    parent_ids = problem['parent_ids']
    ic = problem['ic']

    avg_p_pheno = np.mean([len(p) for p in patients])
    min_p_pheno = np.min([len(p) for p in patients])
    max_p_pheno = np.max([len(p) for p in patients])
    avg_d_pheno = np.mean([len(d) for d in diseases])
    min_d_pheno = np.min([len(d) for d in diseases])
    max_d_pheno = np.max([len(d) for d in diseases])

    print('Test stats:')
    print(f'\tic: {len(ic)}')
    print(f'\tPatients: {len(patients)}, {avg_p_pheno}, {min_p_pheno}, {max_p_pheno}')
    print(f'\tDiseases: {len(diseases)}, {avg_d_pheno}, {min_d_pheno}, {max_d_pheno}')
    #return
    
    #first build the segment tree
    #print('Building SegmentTree...')
    #segment_tree = CustomSegmentTree(parent_ids)

    print('Building SparseTable...')
    sparse_table = STARunner(parent_ids)

    print('Building disease orders...')
    disease_ordered = []
    for i, disease in enumerate(diseases):
        d_id = i+1
        first_occs = [(sparse_table.first_occ[pheno], pheno) for pheno in disease]
        first_occs.sort()
        disease_ordered.append(first_occs)

    pheno_cache = {}
    cache_hits = 0
    cache_miss = 0
    print('Processing...')
    for p_id, patient in enumerate(patients):
        if True or p_id % 100 == 0:
            print(f'[{datetime.datetime.now()}]\t{p_id} {len(patient)} Cache (hit/miss): {cache_hits} {cache_miss}...')
        
        '''
        best_disease = 0
        best_score = 0
        for d_id, disease in enumerate(diseases):
            #print(f'\t\t{d_id} {len(disease)}...')
            disease_score = 0
            for p_pheno in patient:
                max_ic = max([
                    ic[sparse_table.lca_node(p_pheno, d_pheno)-1] for d_pheno in disease
                ])
                disease_score += max_ic
            
            if disease_score > best_score:
                best_score = disease_score
                best_disease = d_id+1
        '''
        '''
        total_scores = np.zeros(dtype='float', shape=(len(diseases), ))
        for p_pheno in patient:
            if p_pheno in pheno_cache:
                cache_hits += 1
                disease_array = pheno_cache[p_pheno]
            else:
                cache_miss += 1
                disease_array = np.zeros(dtype='float', shape=(len(diseases), ))
                for d_id, disease in enumerate(diseases):
                    max_ic = max([
                        ic[sparse_table.lca_node(p_pheno, d_pheno)-1] for d_pheno in disease
                    ])
                    disease_array[d_id] = max_ic
                pheno_cache[p_pheno] = disease_array
            
            total_scores += disease_array
            
        #add one for offset
        best_disease = np.argmax(total_scores)+1
        '''
        best_disease = 0
        best_score = 0
        for d_id, disease in enumerate(diseases):
            #print(f'\t\t{d_id} {len(disease)}...')
            disease_score = 0
            for p_pheno in patient:
                #max_ic = max([
                #    ic[sparse_table.lca_node(p_pheno, d_pheno)-1] for d_pheno in disease
                #])
                lookup = (sparse_table.first_occ[p_pheno], p_pheno)
                ind = bisect.bisect(disease_ordered[d_id], lookup)
                if ind == len(disease_ordered[d_id]):
                    max_ic = ic[sparse_table.lca_node(p_pheno, disease_ordered[d_id][-1][1])-1]
                elif disease_ordered[d_id][ind] == lookup:
                    max_ic = ic[p_pheno-1]
                elif ind == 0:
                    max_ic = ic[sparse_table.lca_node(p_pheno, disease_ordered[d_id][0][1])-1]
                else:
                    max_ic = max(
                        ic[sparse_table.lca_node(p_pheno, disease_ordered[d_id][ind-1][1])-1],
                        ic[sparse_table.lca_node(p_pheno, disease_ordered[d_id][ind][1])-1]
                    )
                
                disease_score += max_ic
            
            #print(d_id, disease_score)
            
            if disease_score > best_score:
                best_score = disease_score
                best_disease = d_id+1

        results.append(best_disease)
    
    return results

def solveProblem5(problem):
    #print(json.dumps(problem, indent=4))
    results = []

    #patients & diseases
    patients = problem['patients']
    diseases = problem['diseases']
    
    #graph stuff
    parent_ids = problem['parent_ids']
    ic = problem['ic']

    avg_p_pheno = np.mean([len(p) for p in patients])
    min_p_pheno = np.min([len(p) for p in patients])
    max_p_pheno = np.max([len(p) for p in patients])
    avg_d_pheno = np.mean([len(d) for d in diseases])
    min_d_pheno = np.min([len(d) for d in diseases])
    max_d_pheno = np.max([len(d) for d in diseases])

    print('Test stats:')
    print(f'\tic: {len(ic)}')
    print(f'\tPatients: {len(patients)}, {avg_p_pheno}, {min_p_pheno}, {max_p_pheno}')
    print(f'\tDiseases: {len(diseases)}, {avg_d_pheno}, {min_d_pheno}, {max_d_pheno}')
    #return

    assert(max_d_pheno == 1)
    assert(max_p_pheno == 1)

    print('Pre-computing disease nodes...')
    disease_nodes = np.zeros(dtype='<u4', shape=(len(ic)+1, ))
    for d_id, dis in enumerate(diseases):
        assert(len(dis) == 1)
        curr_id = dis[0]
        while curr_id != 1:
            if disease_nodes[curr_id] == 0:
                disease_nodes[curr_id] = d_id+1
                curr_id = parent_ids[curr_id-2]
            else:
                break

    print('Processing...')
    for p_id, patient in enumerate(patients):
        if p_id % 10000 == 0:
            print(f'[{datetime.datetime.now()}]\t{p_id}...')
        
        curr_id = patient[0]
        while disease_nodes[curr_id] == 0:
            curr_id = parent_ids[curr_id-2]
        
        results.append(disease_nodes[curr_id])
    
    return results

def solveProblem67(problem):
    #print(json.dumps(problem, indent=4))
    results = []

    #patients & diseases
    patients = problem['patients']
    diseases = problem['diseases']
    
    #graph stuff
    parent_ids = problem['parent_ids']
    ic = problem['ic']

    avg_p_pheno = np.mean([len(p) for p in patients])
    min_p_pheno = np.min([len(p) for p in patients])
    max_p_pheno = np.max([len(p) for p in patients])
    avg_d_pheno = np.mean([len(d) for d in diseases])
    min_d_pheno = np.min([len(d) for d in diseases])
    max_d_pheno = np.max([len(d) for d in diseases])

    print('Test stats:')
    print(f'\tic: {len(ic)}')
    print(f'\tPatients: {len(patients)}, {avg_p_pheno}, {min_p_pheno}, {max_p_pheno}')
    print(f'\tDiseases: {len(diseases)}, {avg_d_pheno}, {min_d_pheno}, {max_d_pheno}')
    #return

    #assert(max_d_pheno <= 6)
    #assert(max_p_pheno <= 6)

    print('Pre-computing sparse matrix...')
    #matrix = csr_matrix((len(ic), len(diseases)), dtype='?')
    matrix = lil_matrix((len(ic), len(diseases)), dtype='<u1')
    for d_id, disease in enumerate(diseases):
        if d_id % 1000 == 0:
            print(d_id)
        for p_id in disease:
            #print(f'\t{p_id}')
            while p_id != 1:
                matrix[p_id-1, d_id] = 1
                p_id = parent_ids[p_id - 2]
            matrix[0, d_id] = 1

    matrix = matrix.tocsr()
    #print((matrix[0, :] + matrix[1, :]))

    print('Processing...')
    for p_id, patient in enumerate(patients):
        best_disease = None
        if True or p_id % 10000 == 0:
            print(f'[{datetime.datetime.now()}]\t{p_id}...')
        
        starting_match = np.array(patient)
        starting_ic = sum([ic[pheno-1] for pheno in patient])

        heap = []
        heapq.heappush(heap, (0-starting_ic, starting_match))

        while len(heap) > 0:
            curr_ic, phenos = heapq.heappop(heap)
            curr_ic = 0-curr_ic
            phenos = np.array(phenos)

            #for p in phenos:
            #    print('\t', p, np.max(matrix[p-1]))
            sumtrix = np.sum(matrix[phenos-1], axis=0)

            #print(sumtrix.shape)
            am = np.argmax(sumtrix)
            max_value = sumtrix[0, am]
            if max_value == len(patient):
                #it's been found
                best_disease = am+1
                break
            else:
                for i, p in enumerate(phenos):
                    new_phenos = copy.deepcopy(phenos)
                    
                    if p == 1:
                        #we can't go "up" anymore on this one
                        pass
                    else:
                        parent_id = parent_ids[p-2]
                        new_phenos[i] = parent_id
                        new_ic = sum([ic[pheno-1] for pheno in new_phenos])
                        heapq.heappush(heap, (0-new_ic, list(new_phenos)))

        results.append(best_disease)
    
    return results
    
def writeResults(fn, all_results):
    fp = open(fn, 'w+')
    for result in all_results:
        fp.write(f'{result}\n')
    fp.close()

if __name__ == '__main__':
    for prob_id in range(0, 8):
        if prob_id == 0:
            problem = 'example'
        else:
            problem = f'test{prob_id}'
        print(f'Analyzing problem set #{problem}...')
        fn = f'/Users/matt/githubProjects/bio_contest_quals_2021/problem_2.4/data/{problem}'
        fno = f'/Users/matt/githubProjects/bio_contest_quals_2021/problem_2.4/results/{problem}.txt'
        problems = loadProblems(fn)
        
        print(f'Solving problems...')
        if prob_id == 5:
            all_results = solveProblem5(problems)
        elif prob_id == 6 or prob_id == 7:
            all_results = solveProblem67(problems)
        else:
            all_results = solveProblem_st(problems)
            
        writeResults(fno, all_results)
        