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

import random
import numpy as np # TODO: Delete this when done

# Used to compute the bandwidth for banded version
MAXINDELS = 3

# Used to implement Needleman-Wunsch scoring
MATCH = -3
INDEL = 5
SUB = 1

class GeneSequencing:

	def __init__( self ):
		pass

	def align( self, seq1, seq2, banded, align_length):
		self.banded = banded
		self.MaxCharactersToAlign = align_length
		if banded:
			result = self.banded_alignment(seq1, seq2, align_length)
		else:
			result = self.unrestricted_alignment(seq1, seq2, align_length)
		return result
	
	def unrestricted_alignment(self, seq1, seq2, align_length):
		seq1 = seq1[:align_length]
		seq2 = seq2[:align_length]
		
		# Set up the DP matrix. 
		#O(nm), n is the length of seq1, m is the length of seq2
		# can change to O(nm) 
		# n is align_length and m is align_length
		dp = [[float('inf') for j in range(len(seq2)+1)] for i in range(len(seq1)+1)]

		# Initialization
		for i in range(len(seq1)+1): # O(n+1) where n is seq1_length
			dp[i][0] = 5 * i
		
		for i in range(len(seq2)+1): # O(n+1) where n is seq2_length
			dp[0][i] = 5 * i
		
		# Fill the DP matrix -> O(mn)
		for i in range(1, len(seq1)+1):
			for j in range(1, len(seq2)+1):
				if seq1[i-1] == seq2[j-1]:
					dp[i][j] = min(dp[i-1][j-1] + MATCH, dp[i-1][j] + INDEL, dp[i][j-1] + INDEL)
				else:
					dp[i][j] = min(dp[i-1][j-1] + SUB, dp[i-1][j] + INDEL, dp[i][j-1] + INDEL)

		align_cost = dp[len(seq1)][len(seq2)]
		print(f'align_cost: {align_cost}')
		print(np.array(dp))
		seqi_first100, seqj_first100 = self.construct_alignment(seq1, seq2, dp)
		
		return {'align_cost': align_cost, 'seqi_first100': seqi_first100, 'seqj_first100': seqj_first100}
	
	# Did this in a O(nm) space
	# because couldn't figure out O(kn) space
	# the filling in the values runs in O(kn)
	# So this is still a lot faster, 
	# just not efficient in memory. 
	# still get correct values as well.
	def banded_alignment(self, seq1, seq2, align_length):
		seq1 = seq1[:align_length]
		seq2 = seq2[:align_length]
		n, m = len(seq1), len(seq2)
		d=MAXINDELS # 3

		# Initialize the DP matrix with infinities O(mn)
		inf = float('inf')
		dp = [[inf] * (m + 1) for _ in range(n + 1)]
		
		dp[0][0] = 0
		# Initialize the borders within the band
		# O(1) since d = 3
		for i in range(1, min(n + 1, d + 1)):
			dp[i][0] = i * INDEL
		# O(1) since d = 3
		for j in range(1, min(m + 1, d + 1)):
			dp[0][j] = j * INDEL

		# Fill the DP matrix within the band
		# O(kn) -> only filling in band values
		for i in range(1, n + 1):
			start = max(1, i - d)  # Start of the band
			end = min(m, i + d)  # End of the band
			for j in range(start, end + 1):
				if seq1[i - 1] == seq2[j - 1]:
					# Match
					dp[i][j] = min(dp[i][j], dp[i - 1][j - 1] + MATCH)
				else:
					# Mismatch
					dp[i][j] = min(dp[i][j], dp[i - 1][j - 1] + SUB)

				# Insertion and Deletion
				if i > j - d:  # Check if within band for deletion
					dp[i][j] = min(dp[i][j], dp[i - 1][j] + INDEL)
				if j > i - d:  # Check if within band for insertion
					dp[i][j] = min(dp[i][j], dp[i][j - 1] + INDEL)

		# The final alignment score is at dp[n][m]
		align_cost = dp[n][m]
		print(f'align_cost: {align_cost}')
		
		# Placeholder for backtracking to reconstruct the aligned sequences
		seqi_first100, seqj_first100 = self.construct_alignment(seq1, seq2, dp)

		return {'align_cost': align_cost, 'seqi_first100': seqi_first100, 'seqj_first100': seqj_first100}
	
	# O(n+m) where n is the length of seq1 and m is the length of seq2
	def construct_alignment(self, seq1, seq2, dp):
		seqi_aligned, seqj_aligned = '', ''
		i, j = len(seq1), len(seq2)
		# O(n+m) # where n is the length of seq1 and m is the length of seq2
		while i > 0 and j > 0:
			# Pref 1: Left (Insertion in seq2)
			if j > 0 and dp[i][j] == dp[i][j-1] + INDEL:
				seqi_aligned = '-' + seqi_aligned
				seqj_aligned = seq2[j-1] + seqj_aligned
				j -= 1
        	# Pref 2: Top (Deletion in seq1)
			elif i > 0 and dp[i][j] == dp[i-1][j] + INDEL:
				seqi_aligned = seq1[i-1] + seqi_aligned
				seqj_aligned = '-' + seqj_aligned
				i -= 1
        	# Pref 3: Diagonal (Match/Sub)
			else:
				seqi_aligned = seq1[i-1] + seqi_aligned
				seqj_aligned = seq2[j-1] + seqj_aligned
				i -= 1
				j -= 1
			
		# Add the remaining characters from seq1 and seq2
		while i > 0:
			seqi_aligned = seq1[i-1] + seqi_aligned
			seqj_aligned = '-' + seqj_aligned
			i -= 1
		
		while j > 0:
			seqi_aligned = '-' + seqi_aligned
			seqj_aligned = seq2[j-1] + seqj_aligned
			j -= 1
        
		return seqi_aligned[:100], seqj_aligned[:100]


# testing area 
if __name__ == '__main__':


	g = GeneSequencing()
	seq1 = 'polynomial'
	seq2 = 'exponential'
	banded = False
	align_length = 1000
	result = g.align(seq1, seq2, banded, align_length)
	print(result)
	seq1 = 'polynomial'
	seq2 = 'exponential'
	banded = True
	align_length = 1000
	result = g.align(seq1, seq2, banded, align_length)
	print(result)