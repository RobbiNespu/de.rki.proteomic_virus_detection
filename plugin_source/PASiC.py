"""
Adaption of similarity correction method from GASiC (plus normalized abundances)

@author: Anke Penzlin, June 2013
"""
import numpy as np
import scipy.optimize as opt

def similarity_correction(sim, reads, N):
	""" Calculate corrected abundances given a similarity matrix and observations using optimization.
	Copied from GASiC, for metaproteomic experiments replace the term 'reads' with 'spectra'.
	
	Input:
	sim:   [numpy.array (M,M)] with pairwise similarities between species
	reads: [numpy.array (M,)] with number of observed reads for each species
	N:     [int] total number of reads
	
	Output:
	abundances: [numpy.array (M,)] estimated abundance of each species in the sample
	"""
	
	# transform reads to abundances and rename similarity matrix
	A = np.transpose(sim)
	r = reads.astype(np.float) / N
	
	rng = range(len(reads))
	
	# Now solve the optimization problem: min_c |Ac-r|^2 s.t. c_i >= 0 and 1-sum(c_i) >= 0
	# construct objective function
	def objective (c):
		# numpy implementation of objective function |A*c-r|^2
		return np.sum(np.square(np.dot(A,c)-r))
	
	# construct constraints
	def build_con(k):
		# constraint: k'th component of x should be >0
		return lambda c: c[k]
	cons = [build_con(k) for k in rng]
	# constraint: the sum of all components of x should not exceed 1
	cons.append( lambda c: 1-np.sum(c) )
	
	# initial guess
	c_0 = np.array([0.5 for i in range(len(reads))])
	
	# finally: optimization procedure
	abundances = opt.fmin_cobyla(objective, c_0, cons, disp=0, rhoend=1e-10, maxfun=10000)
	
	return abundances


def PASiC(simiMat, counts, N):
	"""
	Call similarity_correction and additionally calculate normalized abundances.
	
	Input:
	simiMat: [numpy.array (M,M)] with pairwise similarities between species
	counts:  [numpy.array (M,)] with number of observed PSMs for each species
	N:       [int] total number of spectra
	
	Output:
	abundances: [numpy.array (M,)] estimated abundance of each species in the sample
	normAbundances: [numpy.array (M,)] abundances normalized to one
	"""
	abundances = similarity_correction(sim=simiMat, reads=counts, N=N)
	normAbundances = np.divide(abundances, float(sum(abundances)))
	return abundances, normAbundances

'''if __name__ == "__main__":
	# Idee: np-Dateien einlesen von Eingabepfaden -d common directory
	#c = np.load("/data/NG4/anke/spectra/acid_mine/counts_Biofilm.dat")
	#M = np.load("/data/NG4/anke/spectra/acid_mine/unw_similarityMatrix_seq.dat")
	#M_wtd = np.load("/data/NG4/anke/spectra/acid_mine/wtd_similarityMatrix_seq.dat")
	corrAbunds, normCorr = PASiC(simiMat=M, counts=c, N=10001)
	wtdCorr, n_wtdCorr = PASiC(simiMat=M_wtd, counts=c, N=10001)
	print "observed  abundances: ", np.round(np.divide(c, 10001),5)
	print "observed, normalized: ", np.round(np.divide(c, float(sum(c))),5)
	print "corrected  abundances:", np.round(corrAbunds,5)
	print "corrected, normalized:", np.round(normCorr,5)
	print "with weighting"
	print "corrected  abundances:", np.round(wtdCorr,5)
	print "corrected, normalized:", np.round(n_wtdCorr,5)
'''