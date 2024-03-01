#!/usr/bin/python3


from CS312Graph import *
import time

# helper methods
def print_distances(distances):
    print('Distances\n')
    for node in distances:
        print(f'Node: {node.node_id}, Distance: {distances[node]}')
def print_prev(prev):
    print('Previous\n')
    for node in prev:
        if(prev[node] == None):
            print(f'Node: {node.node_id}, Previous: None')
        else:
            print(f'Node: {node.node_id}, Previous: {prev[node].node_id}')

# Priority Queue using Array
class PriorityQueueArray:
    def __init__(self):
        self.queue = {}

    # Doesn't care about the order of the queue just throws it in.
    def insert(self, node): # O(1) operation
        self.queue[node] = float('inf')

    def decrease_key(self, node, new_val):
        self.queue[node] = new_val # O(1) operation

    def delete_min(self):
        min_node = min(self.queue, key=self.queue.get) # O(n) operation
        del self.queue[min_node]
        return min_node

    def is_empty(self):
        return len(self.queue) == 0 # O(1) operation

    def makequeue(self, nodes, distances): # O(n) operation
        for node in nodes:
            self.insert(node)

    def printqueue(self): # O(n) operation
        for node in self.queue:
            print(f'Node: {node.node_id}, Distance: {self.queue[node]}')

class PriorityQueueHeap:
    def __init__(self):
        self.heap = []  # List of CS312GraphNode
        self.heap_node_positions = {}  # Node_id: index in heap
        self.priorities = {}  # Node_id: priority

    def insert(self, node, priority): # O(logn)
        self.heap.append(node)
        self.heap_node_positions[node.node_id] = len(self.heap) - 1
        self.priorities[node.node_id] = priority
        self.bubble_up(len(self.heap) - 1)

    def decrease_key(self, node, new_priority): # O(logn)
        index = self.heap_node_positions[node.node_id]
        self.priorities[node.node_id] = new_priority
        self.bubble_up(index)
        
    def delete_min(self): # O(logn)
        min_node = self.heap[0]
        last_node = self.heap.pop()
        if self.heap:
            self.heap[0] = last_node
            self.heap_node_positions[last_node.node_id] = 0
            self.sift_down(0)
        self.heap_node_positions.pop(min_node.node_id)
        self.priorities.pop(min_node.node_id)
        return min_node

    def is_empty(self): # O(1)
        return len(self.heap) == 0

    def makequeue(self, nodes, distances): # O(nlogn)
        for node in nodes:
            self.insert(node, distances[node])

    def bubble_up(self, index): # O(logn)
        parent_index = (index - 1) // 2
        while index > 0 and self.compare_nodes(index, parent_index):
            self.swap(index, parent_index)
            index = parent_index
            parent_index = (index - 1) // 2

    def sift_down(self, index): # O(logn)
        smallest = index
        left_index = 2 * index + 1
        right_index = 2 * index + 2

        if left_index < len(self.heap) and self.compare_nodes(left_index, smallest):
            smallest = left_index
        if right_index < len(self.heap) and self.compare_nodes(right_index, smallest):
            smallest = right_index

        if smallest != index:
            self.swap(index, smallest)
            self.sift_down(smallest)

    def compare_nodes(self, index1, index2): # O(1)
        node_id1 = self.heap[index1].node_id
        node_id2 = self.heap[index2].node_id
        return self.priorities[node_id1] < self.priorities[node_id2]

    def swap(self, index1, index2): # O(1)
        self.heap[index1], self.heap[index2] = self.heap[index2], self.heap[index1]
        self.heap_node_positions[self.heap[index1].node_id] = index1
        self.heap_node_positions[self.heap[index2].node_id] = index2
    
    def printqueue(self): # O(n)
        for node in self.heap:
            print(f'Node: {node.node_id}, Distance: {self.priorities[node.node_id]}')


class NetworkRoutingSolver:
    def __init__( self):
        pass

    def initializeNetwork( self, network ):
        assert( type(network) == CS312Graph )
        self.network = network

    def getShortestPath( self, destIndex ):
        self.dest = destIndex
        path_edges = []
        total_length = 0
        node = self.network.nodes[self.dest]

        while node != self.network.nodes[self.source]:
            # No path to node
            if self.prev[node] == None:
                return {'cost':float('inf'), 'path':[]}
            
            prev_node = self.prev[node]
            for edge in prev_node.neighbors:
                if edge.dest == node:
                    path_edges.append( (prev_node.loc, edge.dest.loc, '{:.0f}'.format(edge.length)) )
                    total_length += edge.length
                    break
            node = prev_node
        path_edges.reverse() # Reverse the list to get the correct order

        return {'cost':total_length, 'path':path_edges}

    def computeShortestPaths( self, srcIndex, use_heap=False ):
        self.source = srcIndex
        if use_heap:
            pq = PriorityQueueHeap()
        else:
            pq = PriorityQueueArray()
        
        t1 = time.time()
        self.distances, self.prev = self.dijkstra(self.network.nodes[srcIndex], pq)
        #print_distances(self.distances)
        #print_prev(self.prev)
        t2 = time.time()
        return (t2-t1)
    
    def dijkstra(self, src, pq): # O((V+E)logV) or O(V^2) depending on the priority queue
        # Initialize all distances to infinity
        distances = {}
        prev = {}
        for node in self.network.nodes:
            distances[node] = float('inf')
            prev[node] = None
        distances[src] = 0
        # Insert all nodes into the priority queue
        pq.makequeue(self.network.nodes, distances)
        #pq.printqueue()
        while not pq.is_empty():
            u = pq.delete_min()
            #print(f'Node from delete_min: {u.node_id}, Distance: {distances[u]}')
            for edge in u.neighbors:
                v = edge.dest
                if distances[v] > distances[u] + edge.length:
                    distances[v] = distances[u] + edge.length
                    prev[v] = u
                    pq.decrease_key(v, distances[v])

        return distances, prev
