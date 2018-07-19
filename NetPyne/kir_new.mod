TITLE Kir potassium current for nucleus accumbens (IRK1 = Kir 2.1 - see Mermelstein)

COMMENT 

implemented by HS

ENDCOMMENT


UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
}
 
NEURON {
        SUFFIX kir
        USEION k READ ek WRITE ik
        RANGE  gkbar, ik, mvhalf, mslope, mshift, qfact
}
 
PARAMETER {
	gkbar  = 0.00015 		(mho/cm2)	

	mvhalf = -52		(mV)	
	mslope = 13		(mV)	
	mshift = 30			(mV)	
						
	qfact = 0.5				
    v      (mV)
    ek     (mV)
}

ASSIGNED {
        ik              (mA/cm2)
        gk              (S/cm2)
        minf            (1)     
}
 
STATE { m }
 

BREAKPOINT {
        SOLVE state METHOD cnexp
        gk = gkbar * m
        ik = gk * ( v - ek )
}
  
INITIAL {
	settables(v)
	m = minf
}

FUNCTION_TABLE taumkir(v(mV))  (ms)		

DERIVATIVE state { 
        settables(v)
        m' = (minf - m) / ( taumkir(v)/qfact )
}
 
UNITSOFF
PROCEDURE settables(v) {  
	TABLE minf DEPEND mvhalf, mshift, mslope
		FROM -200 TO 200 WITH 201
			minf = 1  /  ( 1 + exp( (v - mvhalf + mshift) / mslope) )
}
UNITSON
 
 
