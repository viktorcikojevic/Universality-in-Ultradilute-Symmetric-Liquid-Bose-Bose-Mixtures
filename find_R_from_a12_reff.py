import numpy as np
import math




def Gamma(R, V0):
	return R*np.power(V0/2.0, 0.5)

def sclen_func_help(R, gamma):
    return R*(1.0 - np.tan(gamma)/gamma)

def sclen_func(R, V0):
	g = Gamma(R, V0)
	return sclen_func_help(R, g)

def reff_func_help(R, gamma):
	return R*(1.0 + (3.*np.tan(gamma) - gamma*(3.0 + gamma**2)) / (3.*gamma*(gamma-np.tan(gamma))**2))

def reff_func(R, V0):
	g = Gamma(R, V0)
	return reff_func_help(R, g)
	
	
def findU0(aWanted, r0, ms):
	pi = np.pi	
	a = 1.E-08 	* (np.pi/r0)**2 / 2.
	b = (1.-1.E-08) 	* (np.pi/r0)**2 / 2.
	c = (a + b)/2.0
	a_a = 0.0
	a_b = 0.0
	a_c = 0.0
	EPS = 1.E-06
	while(True):
		a_a = sclen_func(r0, a)- aWanted
		a_b = sclen_func(r0, b)- aWanted
		a_c = sclen_func(r0, c)- aWanted
		if( a_a * a_c < 0.0):
			b = c
		if(a_c * a_b < 0.0):
			a = c    	
		ct =  (a + b)/2.0
		if( np.abs( (ct-c)/c  ) < EPS):
			break
		c = ct
	return c
	
def findUR(aWanted, reff_wanted):
	EPS = 1.E-06
	r0_min = 1.E-08
	r0_max = 1.E+08
	r0_mid = (r0_min + r0_max) / 2
	while(True):
		ua = findU0(aWanted, r0_min, 1.)
		ub = findU0(aWanted, r0_max, 1.)
		uc = findU0(aWanted, r0_mid, 1.)
		ra = reff_func(r0_min, ua)	
		rb = reff_func(r0_max, ub)	
		rc = reff_func(r0_mid, uc)
		aa = ra - reff_wanted
		ab = rb - reff_wanted
		ac = rc - reff_wanted
		if( aa * ac < 0.0):
			r0_max = r0_mid
		if(ac * ab < 0.0):
			r0_min = r0_mid		
		r0_mid_t =  (r0_min + r0_max)/2.0
		if( np.abs( (r0_mid_t-r0_mid)/r0_mid  ) < EPS):
			break
		r0_mid = r0_mid_t
	
	return(uc, r0_mid)
	
	
	

a12 = -1.0555
reff = 500000
print("Input params: a12= ", a12, " reff = ", reff)
print("Calculating ... ")
V0, R = findUR(a12, reff)
print("V0 = ", V0, " R= ", R,  "reff = ", reff_func(R, V0), "a12 = ", sclen_func(R, V0))

   	
   	
   	

	
	
	
	
	

