TITLE NaP - persistent sodium current for nucleus accumbens 

COMMENT
Implemented by HS
ENDCOMMENT

UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
}
 
NEURON {
        SUFFIX nap
        USEION na READ ena WRITE ina
        RANGE  gnabar, ina, qfact
}
 
PARAMETER {
	gnabar = 4e-05 (mho/cm2)	
	qfact = 3
	v    				(mV)
	ena					(mV)
	
	
}

ASSIGNED {
        ina		(mA/cm2)
        gna		(mho/cm2)
        minf	(1)
		hinf	(1)	
		mtau	(ms)			
}

 
STATE { m h }
 

BREAKPOINT {
        SOLVE states METHOD cnexp
        gna = gnabar * m * h  
        ina = gna * ( v - ena )

}
 

 
INITIAL {
	settables(v)
	m = minf
	h = hinf
}

FUNCTION_TABLE tauhnap(v(mV))  (ms)		

DERIVATIVE states { 
        settables(v)
        m' = (minf - m) / mtau
        h' = (hinf - h) / (tauhnap(v)/qfact)    
}

UNITSOFF
PROCEDURE settables(v) {  
	TABLE minf, hinf, mtau FROM -120 TO 40 WITH 641

		minf = 1 / (1 + exp( (v +52.6 +10) / (-4.6)))
		hinf = 1 / (1 + exp( (v +48.8 +10) / 10))
		
		if (v < -40) {			
			mtau = 0.025 + 0.14 * exp( (v +40 ) / 10)
		} else {
			mtau = 0.02 + 0.145 * exp( (-v -40) / 10)
		}

}
UNITSON
 
 
