#! /usr/bin/env python

"""
usage: demoselsim.py outfilename popn h pct
"""

import numpy
import sys

def pnext(AC, N, Nout, h, s):
	AC = float(AC); N = float(N); h = float(h); s = float(s)
	p = AC/(2*N)
	w11 = 1+s; w12 = 1+h*s; w22 = 1
	wbar = ((p**2) * w11) + (2*p*(1-p)*w12) + (((1-p)**2) * w22)
	pdet = p*(p*w11 + (1-p)*w12)/wbar
	#print '{0}: {1} {2} {3}'.format(step, N, AC,pdet)
	x=numpy.random.binomial(2*Nout, pdet)
	return x

def trajecvec(h, s, popn, pct):
	if popn == "EUR": #initialize EUR
		tNcurr = 620
		tB1 = 720
		tAdmix = 1900
		N = 10085
		NB1 = 549
		Nfinal = 85
	elif popn == "ASN": #initialize ASN
		N = 10063
		tNcurr = 540 
		NB1 = 407
		tAdmix = 1900
		tB1 = 640
		Nfinal = 97
	
	#get interval points
	first_point = tAdmix-tB1 
	second_point = tAdmix-tNcurr
	AC = int(round(N*2*(pct/100.)))
	AC = [AC]
	
	for t in range(0, first_point): #current population
		AC.append(pnext(AC[t], N, N, h, s))
		
		if (AC[t+1] == 0) or (AC[t+1] == 2*N):
			AC[-1] = 'FIXED:{0}'.format(AC[-1])
			return AC
	
	AC.append(pnext(AC[first_point], N, NB1, h, s))
	
	if (AC[first_point+1] == 0) or (AC[first_point+1] == 2*NB1):
		AC[-1] = 'FIXED:{0}'.format(AC[-1])
		return AC
		
	for t in range(first_point+1, second_point):
		AC.append(pnext(AC[t], NB1, NB1, h, s))
		
		if (AC[t+1] == 0) or (AC[t+1] == 2*NB1):
			AC[-1] = 'FIXED:{0}'.format(AC[-1])
			return AC
			
	AC.append(pnext(AC[second_point], NB1, N, h, s))
	
	if (AC[second_point+1] == 0) or (AC[second_point+1] == 2*N):
		AC[-1] = 'FIXED:{0}'.format(AC[-1])
		return AC
	
	for t in range(second_point+1, tAdmix):
		AC.append(pnext(AC[t], N, N, h, s))
		
		if (AC[t+1] == 0) or (AC[t+1] == 2*N):
			AC[-1] = 'FIXED:{0}'.format(AC[-1])
			return AC
	
	AC.append(pnext(AC[t], N, Nfinal, h, s))
	
	if (AC[-1] == 0) or (AC[t+1] == 2*Nfinal):
		AC[-1] = 'FIXED:{0}'.format(AC[-1])
	
	return AC

def checkint(str):
	try:
		int(str)
		return True
	except ValueError:
		return False

def summary(AC):
	if checkint(AC[-1]):
		return 'SEG:{0}'.format(AC[-1])
	elif AC[-1] == 'FIXED:0':
		time = len(AC)
		return 'LOST:{0}'.format(time)
	else:
		time = len(AC)
		return 'FIXED:{0}'.format(time)

def main():
	outfilename = sys.argv[1]
	popn = sys.argv[2]
	h = float(sys.argv[3])
	pct = int(sys.argv[4])
	
	niter = 1000000
	
	outfile = open(outfilename, 'w')
	
	results_vec = []
	
	for i in xrange(niter):
		s = numpy.random.gamma(0.184, 8200)
		s = s/(2.*25636)
		results_vec.append(summary(trajecvec(h, s, popn, pct)))
		sys.stdout.write('{0}\r'.format(i))
	lost_count = 0; fixed_count = 0; seg_count = 0
	seg_vec = []; lost_vec = []; fixed_vec = []
	
	for result in results_vec:
		state, count = result.split(':')
		if state == 'LOST':
			#print state
			lost_count += 1
			lost_vec.append(int(count))
		elif state == 'SEG':
			#print state
			seg_count += 1
			seg_vec.append(int(count))
		elif state == 'FIXED':
			#print state
			fixed_count += 1
			fixed_vec.append(int(count))
	
	if popn == "EUR": #initialize EUR
		N = 85
	elif popn == "ASN": #initialize ASN
		N = 97
	
	wfreq = numpy.mean([0]*lost_count + seg_vec + [1]*fixed_count)/(2.*N)
	if seg_vec != []:
		segfreq = numpy.mean(seg_vec)/(2.*N)
	else:
		segfreq = 0
	
	s = "dist"
	outstr = '{0} {1} {2} {3} {4}'.format(popn, h, s, wfreq, segfreq)
	outfile.write(outstr)
	
	outfile.close()
	
if __name__ == "__main__":
	main()
