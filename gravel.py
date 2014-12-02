#! /usr/bin/env python

"""
usage: gravel2.py outfilename h s f

h = dominance coefficient
s = selection coefficient (not population-scaled)
f = starting frequency of admixture
"""

import numpy
import sys

def pnext_twopop(AC, N, Nout, h, s, m):
	"""
	AC is a tuple (pAf, pAsEu)
	NAf is a constant, population size taken from Gravel et al. 2011

	function returns a tuple of allele frequencies (pAf, pAsEu)
	"""
	
	#Define 
	NAf = 14474.
	pAf, pAsEu = AC
	NAsEu = float(N); Nout=float(Nout); h = float(h); s = float(s)
	m = float(m)
	
	#get relative fitnesses
	w11 = 1+s; w12 = 1+h*s; w22 = 1
	wbarAf = ((pAf**2) * w11) + (2*pAf*(1-pAf)*w12) + (((1-pAf)**2) * w22)
	wbarAsEu = ((pAsEu**2) * w11) + (2*pAsEu*(1-pAsEu)*w12) + (((1-pAsEu)**2) * w22)
	
	#selection
	pdetAf = pAf*(pAf*w11 + (1-pAf)*w12)/wbarAf
	pdetAsEu = pAsEu*(pAsEu*w11 + (1-pAsEu)*w12)/wbarAsEu
	
	#migration
	pmigAf = pdetAf*(1-m) + pdetAsEu*m
	pmigAsEu = pdetAsEu*(1-m) + pdetAf*m
	
	#drift
	ps = (numpy.random.binomial(2*NAf, pmigAf)/(2*NAf),numpy.random.binomial(2*Nout, pmigAsEu)/(2*Nout))
	
	return ps
	
def pnext_threepop(AC, Ns, Nouts, ms, h, s):
	"""
	AC is a tuple (pAf, pAs, pEu)
	Ns is a tuple (NAs, NEu)
	Nouts is a tuple(Nas_out, NEu_out)
	NAf is a constant
	ms is a tuple(mAfAs, mAfEu, mAsEu)
	
	function returns a tuple of allele frequencies (pAf, pAsEu)
	"""
	
	NAf = 14474.
	
	#make sure everything is float
	pAf, pAs, pEu = map(float, AC)
	NAs, NEu = map(float, Ns)
	NAs_out, NEu_out = map(float, Nouts)
	mAfAs, mAfEu, mAsEu = map(float, ms)
	h = float(h); s = float(s)
	
	#get relative fitnesses
	w11 = 1+s; w12 = 1+h*s; w22=1
	wbarAf = ((pAf**2) * w11) + (2*pAf*(1-pAf)*w12) + (((1-pAf)**2) * w22)
	wbarAs = ((pAs**2) * w11) + (2*pAs*(1-pAs)*w12) + (((1-pAs)**2) * w22)
	wbarEu = ((pEu**2) * w11) + (2*pEu*(1-pEu)*w12) + (((1-pEu)**2) * w22)
	
	#selection
	pdetAf = pAf*(pAf*w11 + (1-pAf)*w12)/wbarAf
	pdetAs = pAs*(pAs*w11 + (1-pAs)*w12)/wbarAs
	pdetEu = pEu*(pEu*w11 + (1-pEu)*w12)/wbarEu
	
	#migration
	pmigAf = pdetAf*(1-(mAfAs+mAfEu)) + pdetAs*mAfAs + pdetEu*mAfEu
	pmigAs = pdetAs*(1-(mAfAs+mAsEu)) + pdetAf*mAfAs + pdetEu*mAsEu
	pmigEu = pdetEu*(1-(mAfEu+mAsEu)) + pdetAf*mAfEu + pdetAs*mAsEu
	
	#drift
	ps = (numpy.random.binomial(2*NAf, pmigAf)/(2*NAf),numpy.random.binomial(2*NAs_out, pmigAs)/(2*NAs_out),numpy.random.binomial(2*NEu_out, pmigEu)/(2*NEu_out))
	
	return ps

def twopop_vec(h, s, pct):
	"""
	takes in args, returns either a string or a tuple (ASN, EUR)
	"""
	#parameters of the Gravel model
	#two population
	Nb = 1861
	mAfB = 15e-5
	tAdmix = 1900
	
	#three population
	tEuAs = 920
	NEu0 = 1032
	rEu = 0.0038
	NAs0 = 554
	rAs = 0.0048
	
	mAfAs = 0.78e-5; mAfEu = 2.5e-5; mAsEu = 3.11e-5
	
	#get time interval points
	t_onepop = tAdmix - tEuAs
	t_twopop = tEuAs
	
	pct = float(pct)/100
	
	AC = (0,pct)
	
	#bottleneck period
	for t in xrange(t_onepop):
		AC = pnext_twopop(AC, Nb, Nb, h, s, mAfB)
	
	#population split, one step to split into two vectors.
	#vector of population sizes = (AC, AC)
	
	#use lambda functions to describe N
	NEu = lambda t:int(round(NEu0 * numpy.exp(rEu*t)))
	NAs = lambda t:int(round(NAs0 * numpy.exp(rAs*t)))
	
	AC = pnext_threepop((AC[0], AC[1], AC[1]), (Nb, Nb), (NAs(0), NEu(0)), (mAfAs, mAfEu, mAsEu), h, s)
	
	#second pulse of admixture, important to know that allele frequencies are capped at 1. Because we assume that all alleles are deleterious or neutral this cap shouldn't really matter.
	if ((AC[1] + pct * 0.2) > 1):
		injectAC = 1
	else:
		injectAC = AC[1] + (pct*0.2)
	
	AC = (AC[0],injectAC,AC[2])
	
	for i in xrange(1, t_twopop):
		AC = pnext_threepop(AC, (NAs(i-1), NEu(i-1)), (NAs(i), NEu(i)), (mAfAs, mAfEu, mAsEu), h, s)
		if (AC == (0,0,0)) or (AC == (1,1,1)):
			return (AC[1], AC[2])
	
	#return pAs, pEu
	return (AC[1], AC[2])

def summary(ps, Ns):
	"""
	ps = tuple of allele frequncies
	Ns = tuple of population sizes
	
	This is the step where we sample down to the 1000 Genomes sizes.
	"""
	finalp = []
	for p, N in zip(ps, Ns):
		finalp.append(numpy.random.binomial(2*N, p)/(2*float(N)))
	return tuple(finalp)	
	
def main():
	outfilename = sys.argv[1]
	h = float(sys.argv[2])
	s = float(sys.argv[3])
	pct = int(sys.argv[4])
	
	niter = 1000000
	#final # of Asn, Eur
	Nf = (97, 85)
	
	outfile = open(outfilename, 'w')
	
	results_vec = []
	
	for i in xrange(niter):
		results_vec.append(summary(twopop_vec(h, s, pct), Nf))
		sys.stdout.write('{0}\r'.format(i+1))
	
	As, Eu = zip(*results_vec)
	means = (numpy.mean(As), numpy.mean(Eu))
	stdevs = (numpy.std(As), numpy.std(Eu))
	numfixed = (sum([x==1 for x in As]), sum([x==1 for x in Eu]))
	numlost = (sum([x==0 for x in As]), sum([x==0 for x in Eu]))
	
	#popn h s f mean_freqx2 R stdevx2 #fixedx2 #lostx2
	outstr = 'h s f ASN EUR R sd_ASN sd_EUR fix_ASN fix_EUR lost_ASN lost_EUR\n'
	outfile.write(outstr)
	outstr = '{0} {1} {2} {3} {4} {5} {6} {7} {8} {9} {10} {11}'.format(h, s, pct, means[0], means[1], means[0]/means[1], stdevs[0], stdevs[1], numfixed[0], numfixed[1], numlost[0], numlost[1])
	outfile.write(outstr)
	
	outfile.close()

if __name__ == "__main__":
	main()