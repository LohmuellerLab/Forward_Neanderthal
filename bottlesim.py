#! /usr/bin/env python

"""
usage: demoselsim.py outfilename EUR/ASN

We need to define
    relative fitnesses: w11, w12, w22
    population sizes: Nanc, Nb1, N1, Nb2, N2
    times: tB1, t1, tB2, t2, tcurr
    starting allele count: ac
    number of iterations: niter
    output file name: outfilename
"""

import numpy
import sys
import multiprocessing

def pnext(p, N, Nout, h, s):
	p = float(p); N = float(N); h = float(h); s = float(s)
	w11 = 1+s; w12 = 1+h*s; w22 = 1
	wbar = ((p**2) * w11) + (2*p*(1-p)*w12) + (((1-p)**2) * w22)
	pdet = p*(p*w11 + (1-p)*w12)/wbar
	#print '{0}: {1} {2} {3}'.format(step, N, AC,pdet)
	x=numpy.random.binomial(2*Nout, pdet)
	return x/(2*N)

def trajecvec(h, s, popn):
	if popn == "EUR": #initialize EUR
		tB = 100
		NB = 549
	elif popn == "ASN": #initialize ASN
		tB = 100
		NB = 407
	
	AC = [0.02]
	
	for t in range(0, tB): #current population
		AC.append(pnext(AC[t], NB, NB, h, s))
		
		if (AC[t+1] == 0) or (AC[t+1] == 1):
			return AC
	
	return AC
        
def main():
	popn = sys.argv[2]
	
	niter = 1000000
	hvec = [0., 0.5, 2.]
	svec = numpy.linspace(-1e-5,-1e-2)
	svec = numpy.sort(numpy.concatenate([svec,[0]]))
	outfilename = sys.argv[1]
	outfile = open(outfilename, 'w')
	outfile.write('popn,h,s,mean_freq,seg,lost,fixed\n')
	
	for h in hvec:
		for s in svec:
			results_vec = []
			
			for i in xrange(niter):
				results_vec.append(trajecvec(h, s, popn)[-1])

			lost = results_vec.count(0)
			fixed = results_vec.count(1)
			seg = niter - lost - fixed
			outstr = '{0},{1},{2},{3},{4},{5},{6}\n'.format(popn, h, s, numpy.mean(results_vec),seg,lost,fixed)
			outfile.write(outstr)
	
	outfile.close()
	
if __name__ == "__main__":
	main()