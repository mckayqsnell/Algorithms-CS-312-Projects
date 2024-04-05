#!/usr/bin/python3

from which_pyqt import PYQT_VER
if PYQT_VER == 'PYQT5':
	from PyQt5.QtCore import QLineF, QPointF
elif PYQT_VER == 'PYQT4':
	from PyQt4.QtCore import QLineF, QPointF
elif PYQT_VER == 'PYQT6':
	from PyQt6.QtCore import QLineF, QPointF
else:
	raise Exception('Unsupported Version of PyQt: {}'.format(PYQT_VER))


import time
import numpy as np
from TSPClasses import *
import heapq
import itertools


class TSPSolver:
	def __init__( self, gui_view ):
		self._scenario = None

	def setupWithScenario( self, scenario ):
		self._scenario = scenario


	''' <summary>
		This is the entry point for the default solver
		which just finds a valid random tour.  Note this could be used to find your
		initial BSSF.
		</summary>
		<returns>results dictionary for GUI that contains three ints: cost of solution,
		time spent to find solution, number of permutations tried during search, the
		solution found, and three null values for fields not used for this
		algorithm</returns>
	'''

	def defaultRandomTour( self, time_allowance=60.0 ):
		results = {}
		cities = self._scenario.getCities()
		ncities = len(cities)
		foundTour = False
		count = 0
		bssf = None
		start_time = time.time()
		while not foundTour and time.time()-start_time < time_allowance:
			# create a random permutation
			perm = np.random.permutation( ncities )
			route = []
			# Now build the route using the random permutation
			for i in range( ncities ):
				route.append( cities[ perm[i] ] )
			bssf = TSPSolution(route)
			count += 1
			if bssf.cost < np.inf:
				# Found a valid route
				foundTour = True
		end_time = time.time()
		results['cost'] = bssf.cost if foundTour else math.inf
		results['time'] = end_time - start_time
		results['count'] = count
		results['soln'] = bssf
		results['max'] = None
		results['total'] = None
		results['pruned'] = None
		return results


	''' <summary>
		This is the entry point for the greedy solver, which you must implement for
		the group project (but it is probably a good idea to just do it for the branch-and
		bound project as a way to get your feet wet).  Note this could be used to find your
		initial BSSF.
		</summary>
		<returns>results dictionary for GUI that contains three ints: cost of best solution,
		time spent to find best solution, total number of solutions found, the best
		solution found, and three null values for fields not used for this
		algorithm</returns>
	'''

	# O(n^2) where n is the number of cities
	def greedy( self, time_allowance=60.0 ):
		start_time = time.time()
		cities = self._scenario.getCities()
		ncities = len(cities)
		visted = [False] * ncities # list of visited cities
		route = []

		# start at the first city (could be any city)
		currrent_city_index = 0
		visted[currrent_city_index] = True
		route.append(cities[currrent_city_index])

		# loop for visting all cities
		# O(n^2) where n is the number of cities
		while len(route) < ncities:
			nearest_city_index = None
			min_distance = np.inf

			# Loop for finding the nearest city
			#O(n) where n is the number of cities
			for i in range(ncities):
				if not visted[i]:
					# get the distance from the current city to the city i
					distance = cities[currrent_city_index].costTo(cities[i])
					# update the nearest city if the distance is smaller
					if distance < min_distance:
						min_distance = distance
						nearest_city_index = i

			# If no nearest city found, break the loop
			if nearest_city_index is None:
				break # No Solution
			else: 
				# update the current city and add it to the route
				currrent_city_index = nearest_city_index
				visted[currrent_city_index] = True
				route.append(cities[currrent_city_index])

		# Check if we made it back to the start
		if route[-1].costTo(route[0]) == np.inf:
			bssf = None
		else: 
			bssf = TSPSolution(route)

		end_time = time.time()
		results = {
			'cost': bssf.cost if bssf is not None else math.inf,
			'time': end_time - start_time,
			'count': 1 if bssf else 0,
			'soln': bssf,
			'max': None,
			'total': None,
			'pruned': None
		}

		return results



	''' <summary>
		This is the entry point for the branch-and-bound algorithm that you will implement
		</summary>
		<returns>results dictionary for GUI that contains three ints: cost of best solution,
		time spent to find best solution, total number solutions found during search (does
		not include the initial BSSF), the best solution found, and three more ints:
		max queue size, total number of states created, and number of pruned states.</returns>
	'''

	def branchAndBound( self, time_allowance=60.0):
		start_time = time.time()
		self.bssf = self.initialize_bssf()

		cities = self._scenario.getCities()
		n = len(cities)
		initial_matrix = self.create_initial_matrix()
		initial_reduced_matrix, lower_bound = self.reduce_matrix(initial_matrix)

		# make the root node
		root_node = Node(level=0, path=[0], reduced_matrix=initial_reduced_matrix, cost=lower_bound)

		# create a priority queue
		pq = []
		heapq.heappush(pq, (lower_bound, -1, root_node))

		# initialize counters
		num_solutions = 0
		max_queue_size = 0
		num_states_created = 1 # including root node
		num_pruned_states = 0

		# loop until the priority queue is empty or time is up
		while pq and (time.time() - start_time) < time_allowance:
			# get the node with the lowest bound
			_ , _ , current_node = heapq.heappop(pq)

			# if current node's cost is less than the bssf, explore it
			if current_node.cost < self.bssf.cost:
				# if the path is complete, update the bssf
				if len(current_node.path) == n:
					# add the cost of the edge back to the start city
					cost_of_edge_back_to_start = current_node.reduced_matrix[current_node.path[-1], current_node.path[0]]
					complete_cost =  current_node.cost + cost_of_edge_back_to_start
					# update the bssf
					if complete_cost < self.bssf.cost or self.bssf is None:
						new_path = current_node.path + [current_node.path[0]]
						self.bssf = TSPSolution([cities[i] for i in new_path])
						num_solutions += 1
				else:
					# create children for each city not yet visited
					for new_city_index in range(1, n):
						# if the city is not in the path
						if new_city_index not in current_node.path:
							# create a child node
							child_node = self.create_child_node(current_node, new_city_index)
							num_states_created += 1
							# if the child's cost is less than the bssf, add it to the queue
							if child_node.cost < self.bssf.cost:
								heapq.heappush(pq, (child_node.cost, -child_node.level, child_node))
							else:
								num_pruned_states += 1
				
				# update the max queue size
				max_queue_size = max(max_queue_size, len(pq))
			else:
				num_pruned_states += 1

		end_time = time.time()
		results = {
			'cost': self.bssf.cost,
			'time': end_time - start_time,
			'count': num_solutions,
			'soln': self.bssf,
			'max': max_queue_size,
			'total': num_states_created,
			'pruned': num_pruned_states,
		}

		return results

	def create_initial_matrix(self):
		cities = self._scenario.getCities()
		n = len(cities)
		matrix = np.full((n, n), np.inf)

		for i in range(n):
			for j in range(i+1, n):
				if i != j:
					cost = cities[i].costTo(cities[j])
					matrix[i, j] = matrix[j, i] = cost

		return matrix

	def initialize_bssf(self):
		# Run the greedy algorithm to get the initial BSSF
		results = self.greedy()
		bssf = results['soln']

		# If the greedy algorithm did not find a solution, run the default random tour
		if bssf is None or bssf.cost == np.inf:
			results = self.defaultRandomTour()
			bssf = results['soln']
		
		return bssf
	
	def reduce_matrix(self, matrix):
		# Reduce rows
		# n = matrix.shape[0]

		# # Reduce rows
		# # O(n^2) where n is the number of cities (rows)
		# for i in range(n):
		# 	row_min = np.min(matrix[i, :]) # O(n)
		# 	# prevent extra operations
		# 	if(row_min != np.inf and row_min != 0):
		# 		matrix[i, :] -= row_min
		# 		bound += row_min

		# # Reduce columns
		# # O(n) where n is the number of cities (columns)
		# for j in range(n):
		# 	col_min = np.min(matrix[:, j]) # O(n)
		# 	# prevent extra operations
		# 	if col_min != np.inf and col_min != 0:
		# 		matrix[:, j] -= col_min
		# 		bound += col_min
  
		# Reduce rows
		bound = 0
		row_min = np.min(matrix, axis=1, keepdims=True)  # Find min in rows
		for i in range(len(row_min)):
			if np.isfinite(row_min[i]):
				matrix[i, :] -= row_min[i]
		# Subtract row min only from finite values
		bound = np.sum(row_min[np.isfinite(row_min)])  # Sum of row reductions, excluding inf

		# Reduce columns
		col_min = np.min(matrix, axis=0)  # Find min in columns
		for j in range(len(col_min)):
			if np.isfinite(col_min[j]):
				matrix[:, j] -= col_min[j]
		bound += np.sum(col_min[np.isfinite(col_min)])  # Sum of column reductions, excluding inf

		return matrix, bound
	
	def create_child_node(self, parent_node, new_city_index):
		# Copy over the parent's reduced matrix
		new_matrix = np.copy(parent_node.reduced_matrix) # O(n^2)
		# print("New Matrix before reduction: ")
		# print(new_matrix)

		# where I came from (row: from)
		new_matrix[parent_node.path[-1], :] = np.inf # O(n)

		# can't go back to where I came from (column: to)
		new_matrix[:, new_city_index] = np.inf # O(n)

		# for selected edge (i, j), set (j, i) to infinity
		# prevents going back to the same city prematurely
		if(len(parent_node.path) > 1):
			new_matrix[new_city_index, parent_node.path[-1]] = np.inf # O(n)

		# Reduce the new matrix and update the bound
		new_matrix, reduction_cost = self.reduce_matrix(new_matrix) # O(n^2)
		# cost of og bound + new cost + cost of taking the edge
		cost_of_edge = parent_node.reduced_matrix[parent_node.path[-1], new_city_index] if parent_node.path[-1] != new_city_index else 0
		new_bound = parent_node.cost + reduction_cost + cost_of_edge

		# new node with new info
		new_path = parent_node.path[:] + [new_city_index] # O(n) since it's a list append
		# new_path = parent_node.path.copy().append(new_city_index)
		# print("Creating child node with path:", new_path)
		new_node = Node(level=parent_node.level + 1, path=new_path, reduced_matrix=new_matrix, cost=new_bound, parent=parent_node)

		return new_node
  


# For easibility with keeping track of the branch and bound tree
class Node:
	def __init__(self, level, path, reduced_matrix, cost, parent=None):
		self.level = level  # Level of node in the state space tree
		self.path = path  # Path taken to reach this node
		self.reduced_matrix = reduced_matrix  # Reduced cost matrix for this node
		self.cost = cost  # Lower bound of the path cost
		self.parent = parent  # Reference to parent node
	def __lt__(self, other):
		return self.cost < other.cost

# For testing purposes

if __name__ == "__main__":
	# Testing reduce matrix
	matrix = np.array([[np.inf, 7, 3, 12], [3, np.inf, 6, 14], [5, 8, np.inf, 6], [9, 3, 5, np.inf]])
	print(matrix)
	
	new_matrix, bound = TSPSolver(None).reduce_matrix(matrix)
	print(new_matrix)
	print(bound)

	# Testing create child node from a parent node with that matrix
	parent_node = Node(0, [0], new_matrix, bound)
	new_city_index = 2
	n = 4
	new_node = TSPSolver(None).create_child_node(parent_node, new_city_index)
	print("New Node Level: ", new_node.level)
	print("New Node Path: ", new_node.path)
	print("New Node Reduced Matrix: ")
	print(new_node.reduced_matrix)
	print("New Node Cost: ", new_node.cost)
	print("Parent Node Path: ", parent_node.path)

	# Test making another child from that child
	new_city_index = 3
	new_node2 = TSPSolver(None).create_child_node(new_node, new_city_index)
	print("New Node Level: ", new_node2.level)
	print("New Node Path: ", new_node2.path)
	print("New Node Reduced Matrix: ")
	print(new_node2.reduced_matrix)
	print("New Node Cost: ", new_node2.cost)
	print("Parent Node Path: ", new_node.path)

	# Test making another child from that child
	new_city_index = 1
	new_node3 = TSPSolver(None).create_child_node(new_node2, new_city_index)
	print("New Node Level: ", new_node3.level)
	print("New Node Path: ", new_node3.path)
	print("New Node Reduced Matrix: ")
	print(new_node3.reduced_matrix)
	print("New Node Cost: ", new_node3.cost)
	print("Parent Node Path: ", new_node2.path)
