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
		overall_best_cost = np.inf
		overall_best_solution = None
		
		for starting_city_index in range(ncities):
			visted = [False] * ncities # list of visited cities
			route = []
			currrent_city_index = starting_city_index
			visted[currrent_city_index] = True
			route.append(cities[currrent_city_index])

			# loop for visting all cities
			# O(n^2)
			while len(route) < ncities and time.time() - start_time < time_allowance:
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
					break # no solution from city, try another city
				else: 
					# update the current city and add it to the route
					currrent_city_index = nearest_city_index
					visted[currrent_city_index] = True
					route.append(cities[currrent_city_index])
			if len(route) == ncities and route[-1].costTo(route[0]) != np.inf:
				bssf = TSPSolution(route)
				print(f"Found a solution in greedy with cost: {bssf.cost}")
				if bssf.cost < overall_best_cost:
					overall_best_cost = bssf.cost
					overall_best_solution = bssf
		
		if(overall_best_solution is None):
			print("No solution found in greedy with any starting city")
			bssf = None
		else:
			bssf = overall_best_solution

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

	def branchAndBound( self, time_allowance=60.0):
		start_time = time.time()
		bssf = self.initialize_bssf()

		cities = self._scenario.getCities()
		n = len(cities)
		initial_matrix = self.create_initial_matrix()
		initial_reduced_matrix, lower_bound = self.reduce_matrix(initial_matrix)
		
		# make the root node
		root_node = Node(level=0, path=[0], reduced_matrix=initial_reduced_matrix, cost=lower_bound)

		# create a priority queue
		pq = []
		heapq.heappush(pq, (lower_bound, 0, root_node))

		# initialize counters
		num_solutions = 0
		max_queue_size = 0
		num_states_created = 1 # including root node
		num_pruned_states = 0

		# loop until the priority queue is empty or time is up
		while pq and (time.time() - start_time) < time_allowance:
			# get the node with the lowest bound
			_ , _ , current_node = heapq.heappop(pq)

			if current_node.cost >= bssf.cost or len(current_node.path) > n:
				num_pruned_states += 1
				continue

			if len(current_node.path) == n and current_node.reduced_matrix[current_node.path[-1], 0] != np.inf:
				complete_cost = current_node.cost + current_node.reduced_matrix[current_node.path[-1], 0]
				if complete_cost < bssf.cost:
					new_path = current_node.path
					bssf = TSPSolution([cities[i] for i in new_path])
					num_solutions += 1
				continue

			# O(n)
			for new_city_index in range(1, n):
				if new_city_index not in current_node.path:
					child_node = self.create_child_node(current_node, new_city_index)
					num_states_created += 1

					if child_node.cost < bssf.cost and child_node.cost != np.inf:
						# inverted the level to prioritize nodes with lower cost and higher depth
						heapq.heappush(pq, (child_node.cost, -child_node.level, child_node))
					else:
						num_pruned_states += 1
			# O(n)?
			max_queue_size = max(max_queue_size, len(pq))

		end_time = time.time()
		results = {
			'cost': bssf.cost if bssf is not None else math.inf,
			'time': end_time - start_time,
			'count': num_solutions,
			'soln': bssf,
			'max': max_queue_size,
			'total': num_states_created,
			'pruned': num_pruned_states,
		}

		return results

	# O(n^2)
	def create_initial_matrix(self):
		cities = self._scenario.getCities()
		n = len(cities)
		# if an edge doesn't exist, it's infinity
		matrix = np.full((n, n), np.inf) 

		# dont assume symmetric matrix
		for i in range(n):
			for j in range(n):
				if i != j:
					cost = cities[i].costTo(cities[j])
					matrix[i, j] = cost

		return matrix
	
	# O(n^2) --> because of greedy
	def initialize_bssf(self):
		# Run the greedy algorithm to get the initial BSSF
		results = self.greedy()
		bssf = results['soln']

		return bssf
	
	# O(n^2)
	def reduce_matrix(self, matrix):
		# Reduce rows
		bound = 0
		# row_min matrix (1D array)
		row_min = np.min(matrix, axis=1, keepdims=True)
		for i in range(len(row_min)):
			# Subtract row min only from finite values
			if np.isfinite(row_min[i]):
				matrix[i, :] -= row_min[i]
		# Sum of row reductions, excluding inf
		bound = np.sum(row_min[np.isfinite(row_min)])

		# Reduce columns
		# col_min matrix(1D array)
		col_min = np.min(matrix, axis=0)
		for j in range(len(col_min)):
			if np.isfinite(col_min[j]):
				matrix[:, j] -= col_min[j]
		# Sum of column reductions, excluding inf
		bound += np.sum(col_min[np.isfinite(col_min)])

		return matrix, bound
	
	# O(n^2)?
	def create_child_node(self, parent_node, new_city_index):
		# Copy over the parent's reduced matrix
		new_matrix = np.copy(parent_node.reduced_matrix) # O(n^2)
		# where I came from (row: from)
		new_matrix[parent_node.path[-1], :] = np.inf # O(n)
		# can't go back to where I came from (column: to)
		new_matrix[:, new_city_index] = np.inf # O(n)

		# for selected edge (i, j), set (j, i) to infinity
		# prevents going back to the same city prematurely
		if(len(parent_node.path) > 1):
			new_matrix[new_city_index, parent_node.path[-1]] = np.inf # O(n)

		new_matrix, reduction_cost = self.reduce_matrix(new_matrix) # O(n^2)

		cost_of_edge = parent_node.reduced_matrix[parent_node.path[-1], new_city_index]
		new_bound = parent_node.cost + reduction_cost + cost_of_edge
  		# O(n) since it's a list append
		new_path = parent_node.path[:] + [new_city_index]
		# new node with new info
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
		# priortize nodes with lower cost and higher dpeth when they have the same cost
		if self.cost == other.cost:
			return self.level > other.level
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
