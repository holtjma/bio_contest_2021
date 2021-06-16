
import math
import numpy as np
import time

from segment_tree import *

class CustomSegmentTree:
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
        #print(self.traversal)
        #print(self.first_occ)

        #self.min_traversal = [self.depth[t] for t in self.traversal]
        self.min_traversal = self.traversal
        self.search_tree = SegmentTree(self.min_traversal)

        #print(self.min_traversal)
        '''

        #now build the search tree as a static array
        log_size = math.ceil(math.log(len(self.min_traversal), 2))+1
        print('log_size:', log_size)
        arr_size = 2**log_size - 1
        self.search_tree = np.zeros(dtype='<u8', shape=(arr_size, ))
        self.search_tree[:] = 0xFFFFFFFFFFFFFFFF

        copy_point = 2**(log_size-1)-1
        self.search_tree[copy_point:copy_point+len(self.min_traversal)] = self.min_traversal

        for x in range(copy_point-1, -1, -1):
            c1 = 2*x+1
            c2 = 2*x+2
            self.search_tree[x] = min(self.search_tree[c1], self.search_tree[c2])

        #print(self.search_tree)
        '''
    
    def lca_node(self, x, y):
        fc1 = self.first_occ[x]
        fc2 = self.first_occ[y]
        if fc2 < fc1:
            f = fc1
            fc1 = fc2
            fc2 = f
        
        #super naive slow mode
        #min_index = int(np.min(self.traversal[fc1:fc2+1]))
        
        curr_min = self.search_tree.query(fc1, fc2, 'min')

        '''
        #node, start_range, end_range
        stack = [(0, 0, self.search_tree.shape[0] // 2)]
        curr_min = 0xFFFFFFFFFFFFFFFF
        
        while len(stack) > 0:
            node, start, end = stack.pop()
            if start >= fc1 and end <= fc2:
                #fully contained
                curr_min = min(curr_min, self.search_tree[node])
            elif start > fc2 or end < fc1:
                #fully outside
                pass
            else:
                #not fully contained
                midpoint = (start+end) // 2
                stack.append(
                    (2*node+1, start, midpoint)
                )
                stack.append(
                    (2*node+2, midpoint+1, end)
                )
                #print('\t', (2*node+1, start, midpoint), (2*node+2, midpoint, end))
        '''

        return int(curr_min)
