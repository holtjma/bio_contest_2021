
import math
import numpy as np
import time

'''
Base on "Method 3 (Sparse Table Algorithm)" from https://www.geeksforgeeks.org/range-minimum-query-for-static-array/
'''

class STARunner:
    def __init__(self, parent_ids):
        self.parent_ids = parent_ids
        
        #parse into children
        children = {
            1 : []
        }
        self.depth = [-1]*(len(self.parent_ids)+2)
        self.depth[0:2] = [0, 0]
        for x, p_id in enumerate(self.parent_ids):
            c_id = x+2
            children[p_id].append(c_id)
            children[c_id] = []
            self.depth[c_id] = self.depth[p_id]+1
        assert(-1 not in self.depth)
        
        #now create an in-order traversal
        start_node = 1
        self.traversal = []
        self.first_occ = [-1]*(len(self.parent_ids)+2)

        visited = set([])
        stack = [start_node]
        while len(stack) > 0:
            next_node = stack.pop()
            self.traversal.append(next_node)
            if next_node not in visited:
                self.first_occ[next_node] = len(self.traversal)-1
                visited.add(next_node)
                for c_id in sorted(children[next_node])[::-1]:
                    stack.append(next_node)
                    stack.append(c_id)
        assert(-1 not in self.first_occ[1:])
        
        #self.min_traversal = [self.depth[t] for t in self.traversal]
        #self.min_traversal = self.traversal
        
        #print(self.min_traversal)

        self.traversal = np.array(self.traversal)
        
        trav_len = self.traversal.shape[0]
        log_size = int(np.floor(np.log2(trav_len)))+1
        self.lookup_table = np.zeros(dtype='<u4', shape=(self.traversal.shape[0], log_size))

        self.lookup_table[:, 0] = self.traversal
        for x in range(1, log_size):
            log_factor = 2**x
            offset = 2**(x-1)
            for y in range(0, trav_len-log_factor+1):
                min_check_left = self.lookup_table[y][x-1]
                min_check_right = self.lookup_table[y+offset][x-1]
                self.lookup_table[y, x] = min(min_check_left, min_check_right)
    
    def lca_node(self, x, y):
        fc1 = self.first_occ[x]
        fc2 = self.first_occ[y]
        if fc2 < fc1:
            f = fc1
            fc1 = fc2
            fc2 = f
        
        j_offset = int(np.floor(np.log2(fc2-fc1+1)))
        lookup1 = self.lookup_table[fc1, j_offset]
        lookup2 = self.lookup_table[fc2-(2**j_offset)+1, j_offset]
        return min(lookup1, lookup2)
