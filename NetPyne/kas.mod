TITLE Potassium A-type current for nucleus accumbens (Kv1.2)

COMMENT
Jason Moyer 2004 - jtmoyer@seas.upenn.edu

Shen W, Hernandez-Lopez S, Tkatch T, Held JE, Surmeier DJ (2004).
Kv1.2-containing k+ channels regulate subthreshold excitability of
striatal medium spiny neurons. J Neurophys 91: 1337-1349.

Some of the parameters were published incorrectly in the original paper - these
are detailed below

ENDCOMMENT

UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
        (S)  = (siemens)
}
 
NEURON {
        SUFFIX kas
        USEION k READ ek WRITE ik
        RANGE  gkbar, ik, qfact, mshift, hshift
}
 
PARAMETER {
    gkbar   =   0.01 (mho/cm2)	: 0.01 soma&prox; 0.00091483 mid&dist

	qfact = 9				: qfact = 3 after equations were set but before 
							:	temp correction for fig 6E; QF = 9 is 35 degC

	vmh = -27.0	(mV)		: Shen 2004
	vmc = -16	(mV)		: Shen 2004	
	
	vhh = -33.5	(mV)		: Shen 2004 
	vhc = 21.5	(mV)		: Shen 2004
	
	taum0 = 3.4	(ms)		: Shen 2004
	Cm = 89.2	(ms)		: Shen 2004
	vthm = -34.3	(mV)	: Shen 2004
	vtcm = 30.1	(mV)		: Shen 2004
	
	alpha = 1				: correspondence with josh held
	vth1 = -0.96	(mV)	: Shen 2004
	vtc1 = 29.01	(mV)	: Shen 2004
	
	beta = 1				: josh held
	vth2 = -0.96	(mV)	: Shen 2004
	vtc2 = 100 	(mV)		: Shen 2004
	
	Ch = 9876.6	(ms)		: 548.7 * 18, set to match tauh = 2 s with div by qfact = 3
	a = 0.996				: josh held
	hshift = 0		(mV)
	mshift = 0		(mV)
	htaushift = -90	(mV)	: to correct kinetics
}
 
STATE { m h }
 
ASSIGNED {
		ek				(mV)
        v 				(mV)
        ik 				(mA/cm2)
        gk				(S/cm2)
        minf
	hinf
        mtau		(ms)
        htau		(ms)
   }
  
INITIAL {
	settables(v)
	m = minf
	h = hinf

}

BREAKPOINT {
        SOLVE state METHOD cnexp
        gk = gkbar * m * m * (0.996*h + (1-0.996)) 
        ik = gk * ( v - ek )
}

DERIVATIVE state { 
        settables(v)
	mtau = mtau / qfact
	htau = htau / qfact
        m' = (minf - m)/mtau
        h' = (hinf - h)/htau
}

PROCEDURE settables( v (mV) ) {
	LOCAL left, right

	TABLE minf, hinf, mtau, htau DEPEND hshift, Ch
		FROM -200 TO 200 WITH 201
		
	  	minf = 1 / (1+(exp( (v +27 - mshift) / (-16) )))
  		hinf = 1 / (1+(exp( (v +33.5 - hshift) / 21.5 )))
  		
 		mtau = 3.4  +  89.2 * exp( - ((v+34.3)/30.1)^2 )

		left = 1 * exp( -(v+0.96+90)/29.01 )	: originally exp((v-vth1)/vtc1)
		right = 1 * exp( (v+0.96+90)/100 )	: originally exp(-(v-vth2)/vtc2)
		htau = 9876.6  /  ( left + right )
}





 
