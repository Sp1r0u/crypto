######################################
# 

def verif(E, VerifK, Proof, IO):

    c0 = IO[0]
    c1 = IO[1]
    c2 = IO[2]
    c3 = IO[3]
    c4 = IO[4]
    c7 = IO[5]

    ###print c0, c1, c2, c3, c4, c7

    Fp = E.base_field()
    Fp = Fp.base()
    ###print 'Fp',Fp

    P = VerifK['P']
    n = P.order()
    p = Fp.order()
    k = GF(n)(p).multiplicative_order()

    print '*********************************'
    print '*** Verifier checks the proof ***'
    print '*********************************'
    
    ### EQ(1) ###
    Q = VerifK['Q']
    alphavQ = VerifK['alphavQ']
    rvVmidP = Proof['rvVmid(s)P']
    alphavrvVmidP = Proof['alphavrvVmid(s)P']

    RHS = alphavrvVmidP.tate_pairing(Q, n, k)
    LHS = rvVmidP.tate_pairing(alphavQ, n, k)

    if RHS==LHS: print '   EQ(1):OK'
    else:        print '   EQ(1):NOT OK'

    ###print 'RHS:',RHS
    ###print '    ',my_tate_pairing(E, alphavrvVmidP, Q)

    ###print 'LHS:',LHS
    ###print '    ',my_tate_pairing(E, rvVmidP, alphavQ)

    ### EQ(2) ###
    alphawP = VerifK['alphawP']
    rwWmidQ = Proof['rwWmid(s)Q']
    alphawrwWmidP = Proof['alphawrwWmid(s)P']
    
    RHS = alphawrwWmidP.tate_pairing(Q, n, k)
    LHS = alphawP.tate_pairing(rwWmidQ, n, k)

    if RHS==LHS: print '   EQ(2):OK'
    else:        print '   EQ(2):NOT OK'

    ### EQ(3) ###
    alphayQ = VerifK['alphayQ']
    ryYmidP = Proof['ryYmid(s)P']
    alphayryYmidP = Proof['alphayryYmid(s)P']

    RHS = alphayryYmidP.tate_pairing(Q, n, k)
    LHS = ryYmidP.tate_pairing(alphayQ, n, k)

    if RHS==LHS: print '   EQ(3):OK'
    else:        print '   EQ(3):NOT OK'

    ### EQ(4) ###
    ZmidP = Proof['Zmid(s)P']
    betaQ = VerifK['betaQ']
    betaP = VerifK['betaP']

    RHS  = (rvVmidP+ryYmidP).tate_pairing(betaQ, n, k)
    RHS *= betaP.tate_pairing(rwWmidQ, n, k)
    LHS  = ZmidP.tate_pairing(Q, n, k)
    
    if RHS==LHS: print '   EQ(4):OK'
    else:        print '   EQ(4):NOT OK'

    ### EQ(5) ###
    rytP = VerifK['ryt(s)P']
    hQ   = Proof['h(s)Q']

    c0rvV0P = c0*VerifK['rvV0(s)P']
    c2rvV2P = c2*VerifK['rvV2(s)P']

    c1rwW1Q = c1*VerifK['rwW1(s)Q']
    c3rwW3Q = c3*VerifK['rwW3(s)Q']
    c4rwW4Q = c4*VerifK['rwW4(s)Q']
    
    c7ryY7P = c7*VerifK['ryY7(s)P']

    rvVP = c0rvV0P + c2rvV2P + Proof['rvVmid(s)P']
    rwWQ = c1rwW1Q + c3rwW3Q + c4rwW4Q + Proof['rwWmid(s)Q']
    ryYP = c7ryY7P + Proof['ryYmid(s)P'] 

    RHS  = rvVP.tate_pairing(rwWQ, n, k)
    LHS  = ryYP.tate_pairing(Q, n, k)
    LHS *= rytP.tate_pairing(hQ, n, k)	

    if RHS==LHS: print '   EQ(5):OK'
    else:        print '   EQ(5):NOT OK'

    ###print '    ',my_tate_pairing(E, ryYP, Q)
    ###GT = [ ( P.tate_pairing(Q, n, k) )^i for i in range(P.order())]
    ###print '   GT:',GT

    return True  

####################################
# c: fInput
# h: polynomial s.t. p(x) = h(x)t(x)

def generatePOC(E, x5, x6, x7, c, c5, c6, c7, EvalK):
###    Fp = E.base_field()
###    R  = PolynomialRing(Fp, 'x')
###    p  = Fp.order()

    p  = int(P.order())
    Fp = GF(p)
    R  = PolynomialRing(Fp, 'x')

    xlist = []
    ylist = []
    ###print 'R:',R
    for i in range(4):
    	foo = False
	while foo == False:
		x_i = int(Fp.random_element())
		if (x_i!=0) and (x_i!=x5) and (x_i!=x6) and (x_i!=x7): foo = True
	dummy1 = getTarget(x_i, x5, x6, x7, p)
	dummy1 = inverse_mod(dummy1, p)
	dummy2 = getP(x_i, x5, x6, x7, p, c, c5, c6, c7)
	y_i  = dummy2 * dummy1
	y_i %= Fp.order()
	xlist.append(x_i)
	ylist.append(y_i) 
    h = R.lagrange_polynomial([(xlist[i], ylist[i]) for i in range(4)])
    ###h = R.lagrange_polynomial([(xlist[i], ylist[i]) for i in range(len(xlist))])
    h_coeff = h.coefficients()
    print '   h(x)=',h
    ###print '   h_i:',h_coeff[0], ' ', h_coeff[1]

    print '******************************'
    print '*** Worker generates proof ***'
    print '******************************'
    
    PI = {}
    PI['rvVmid(s)P']       = c5*EvalK['rvV5(s)P']
    PI['alphavrvVmid(s)P'] = c5*EvalK['alphavrvV5(s)P']

    PI['rwWmid(s)Q']	   = c6*EvalK['rwW6(s)Q']
    PI['alphawrwWmid(s)P'] = c6*EvalK['alphawrwW6(s)P']

    PI['ryYmid(s)P']       = c5*EvalK['ryY5(s)P'] + c6*EvalK['ryY6(s)P']
    PI['alphayryYmid(s)P'] = c5*EvalK['alphayryY5(s)P'] + c6*EvalK['alphayryY6(s)P']

    PI['Zmid(s)P'] = c5*EvalK['z5(s)P'] + c6*EvalK['z6(s)P']

    PI['h(s)Q'] = int(h_coeff[0])*EvalK['s^0Q'] + int(h_coeff[1])*EvalK['s^1Q']

    print '   PROOF PI:'
    print '      ',PI

    return PI

#############################################
# c5, c6, c7: values of the circuit wires 

def evaluateQAP(x5, x6, x7, E, c):
    ###Fp = E.base_field()
    ###p  = Fp.order()
    p  = P.order()
    Fp = GF(p)
    c5, c6, c7 = var('c5, c6, c7')

    print '****************************'
    print '*** Worker evaluates QAP ***'
    print '****************************'
    print '   input:',fInput
    eq1 = getP(x5, x5, x6, x7, p, c, c5, c6, c7) == 0
    eq2 = getP(x6, x5, x6, x7, p, c, c5, c6, c7) == 0
    eq3 = getP(x7, x5, x6, x7, p, c, c5, c6, c7) == 0
    print '   eq1:',eq1
    print '   eq2:',eq2
    print '   eq3:',eq3
    solns = solve([eq1, eq2, eq3], c5, c6, c7)
    print '   solns: c5 =',(int(c5.subs(solns)))%p
    print '          c6 =',(int(c6.subs(solns)))%p
    print '          c7 =',(int(c7.subs(solns)))%p
    return ( (int(c5.subs(solns)))%p, (int(c6.subs(solns)))%p, (int(c7.subs(solns)))%p )

#####################################################################
# P: generator of the n-torsion group G1
# Q: generator of the n-torsion group G2
# PubK: public key dictionary.
# 	PubK contains the evaluation and the public verification keys

def generatePublicParams(s, P, Q, x5, x6, x7, E):
    ###Fp = E.base_field()
    ###p  = Fp.order()
    p  = P.order()
    Fp = GF(p)

    print '**********************************'
    print '*** Client generates EK and VK ***'
    print '**********************************'
    
    rv=rw=alphav=alphaw=alphay=beta=gamma=0
    
    while( rv==0     ): rv     = int(Fp.random_element())
    while( rw==0     ): rw     = int(Fp.random_element())
    while( alphav==0 ): alphav = int(Fp.random_element())
    while( alphaw==0 ): alphaw = int(Fp.random_element())
    while( alphay==0 ): alphay = int(Fp.random_element())
    while( beta==0   ): beta   = int(Fp.random_element())
    while( gamma==0  ): gamma  = int(Fp.random_element())

    ry = (rv*rw)%p
        
    print '   SET of RANDOM ELEMENTS:'
    print '      (rv, rw, ry, alphav, alphaw, alphay, beta, gamma):',rv,rw,ry,alphav,alphaw,alphay,beta,gamma

    VK = {}

    VK['P'] = P
    VK['Q'] = Q

    VK['alphavQ']  = alphav*Q
    VK['alphawQ']  = alphaw*Q
    VK['alphawP']  = alphaw*P
    VK['alphayQ']  = alphay*Q

    VK['betaP']    = beta*P
    VK['betaQ']    = beta*Q
    
    VK['ryT(s)P']  = (ry*getTarget(s, x5, x6, x7, p)%p)*P
    
    VK['rvV0(s)P'] = (rv*getV0(s, x5, x6, x7, p)%p)*P
    VK['rvV2(s)P'] = (rv*getV2(s, x5, x6, x7, p)%p)*P

    VK['rwW1(s)Q'] = (rw*getW1(s, x5, x6, x7, p)%p)*Q
    VK['rwW3(s)Q'] = (rw*getW3(s, x5, x6, x7, p)%p)*Q
    VK['rwW4(s)Q'] = (rw*getW4(s, x5, x6, x7, p)%p)*Q 
  
    VK['ryY7(s)P'] = (ry*getY7(s, x5, x6, x7, p)%p)*P

    VK['ryt(s)P']  = (ry*getTarget(s, x5, x6, x7, p)%p)*P

    print '   PUBLIC VERIFICATION KEY VK:'
    print '      ',VK
    ###printDictionary( VK )

    EK = {}
    EK['rvV5(s)P']       = (rv*getV5(s, x5, x6, x7, p)%p)*P
    EK['alphavrvV5(s)P'] = (alphav*rv*getV5(s, x5, x6, x7, p)%p)*P
    
    EK['rwW6(s)Q']       = (rw*getW6(s, x5, x6, x7, p)%p)*Q
    EK['alphawrwW6(s)P'] = (alphaw*rw*getW6(s, x5, x6, x7, p)%p)*P
    
    EK['ryY5(s)P']       = (ry*getY5(s, x5, x6, x7, p)%p)*P
    EK['ryY6(s)P']       = (ry*getY6(s, x5, x6, x7, p)%p)*P
    EK['alphayryY5(s)P'] = (alphay*ry*getY5(s, x5, x6, x7, p)%p)*P
    EK['alphayryY6(s)P'] = (alphay*ry*getY6(s, x5, x6, x7, p)%p)*P

    EK['z5(s)P'] = (rv*beta*getV5(s, x5, x6, x7, p)+ry*beta*getY5(s, x5, x6, x7, p)%p)*P
    EK['z6(s)P'] = (rw*beta*getW6(s, x5, x6, x7, p)+ry*beta*getY6(s, x5, x6, x7, p)%p)*P

    EK['s^0Q'] = (s^0)*Q
    EK['s^1Q'] = (s^1)*Q

    print '   PUBLIC EVALUATION KEY EK:'
    print '      ',EK

    return (EK, VK)


###################################################
# Fp: base field
# (x5, x6, x7): set of multiplication gates labeled
#      	   	by the index of their output values

def generateQAP(E, s):
    ###Fp = E.base_field()
    p  = P.order()
    Fp = GF(p)
    x5 = getLabel(Fp)
    x6 = getLabel(Fp)
    x7 = getLabel(Fp)

    t  = getTarget(s, x5, x6, x7, Fp.order())

    print '****************************'
    print '*** Client generates QAP ***'
    print '****************************'
    print '   PRIVATE'
    print '      secret s:',s
    ###print '      target poly. t(s):',t
    print '   PUBLIC'
    print '      mult. gate labels x5:',x5
    print '                        x6:',x6
    print '                        x7:',x7
    
    return (x5, x6, x7)

####################################
# (V_i, W_i, Y_i) set of polynomials
# where i=0,..,7

def getV0(x, x5, x6, x7, p):
    foo1 = (x7-x5)*(x7-x6)
    foo1 = inverse_mod(foo1, p)
    foo2 = (x-x5)*(x-x6)
    return (foo1*foo2)%p
    
def getV2(x, x5, x6, x7, p):
    foo1 = (x5-x6)*(x5-x7)
    foo1 = inverse_mod(foo1, p)
    foo2 = (x-x6)*(x-x7)
    return (foo1*foo2)%p
    
def getV5(x, x5, x6, x7, p):
    foo1 = (x6-x5)*(x6-x7)
    foo1 = inverse_mod(foo1, p)
    foo2 = (x-x5)*(x-x7)
    return (foo1*foo2)%p
    
def getW1(x, x5, x6, x7, p):
    foo1 = (x7-x5)*(x7-x6)
    foo1 = inverse_mod(foo1, p)
    foo2 = (x-x5)*(x-x6)
    return (foo1*foo2)%p
    
def getW3(x, x5, x6, x7, p):
    foo1 = (x5-x6)*(x5-x7)
    foo1 = inverse_mod(foo1, p)
    foo2 = (x-x6)*(x-x7)
    foo3 = foo1*foo2
    dummy1 = (x6-x5)*(x6-x7)
    dummy1 = inverse_mod(dummy1, p)
    dummy2 = (x-x5)*(x-x7)
    dummy3 = dummy1*dummy2
    return (foo3+dummy3)%p
    
def getW4(x, x5, x6, x7, p):
    foo1 = (x6-x5)*(x6-x7)
    foo1 = inverse_mod(foo1, p)
    foo2 = (x-x5)*(x-x7)
    return (foo1*foo2)%p
    
def getW6(x, x5, x6, x7, p):
    foo1 = (x7-x5)*(x7-x6)
    foo1 = inverse_mod(foo1, p)
    foo2 = (x-x5)*(x-x6)
    return (foo1*foo2)%p
    
def getY5(x, x5, x6, x7, p):
    foo1 = (x5-x6)*(x5-x7)
    foo1 = inverse_mod(foo1, p)
    foo2 = (x-x6)*(x-x7)
    return (foo1*foo2)%p
    
def getY6(x, x5, x6, x7, p):
    foo1 = (x6-x5)*(x6-x7)
    foo1 = inverse_mod(foo1, p)
    foo2 = (x-x5)*(x-x7)
    return (foo1*foo2)%p
    
def getY7(x, x5, x6, x7, p):
    foo1 = (x7-x5)*(x7-x6)
    foo1 = inverse_mod(foo1, p)
    foo2 = (x-x5)*(x-x6)
    return (foo1*foo2)%p
    
def getP(x, x5, x6, x7, p, c, c5, c6, c7):
    V0 = getV0(x, x5, x6, x7, p)
    V2 = getV2(x, x5, x6, x7, p)
    V5 = getV5(x, x5, x6, x7, p)
    W1 = getW1(x, x5, x6, x7, p)
    W3 = getW3(x, x5, x6, x7, p)
    W4 = getW4(x, x5, x6, x7, p)
    W6 = getW6(x, x5, x6, x7, p)
    Y5 = getY5(x, x5, x6, x7, p)
    Y6 = getY6(x, x5, x6, x7, p)
    Y7 = getY7(x, x5, x6, x7, p)
    foo1 = c[0]*V0+c[2]*V2+c5*V5
    foo2 = c[1]*W1+c[3]*W3+c[4]*W4+c6*W6
    foo3 = c5*Y5+c6*Y6+c7*Y7
    return (foo1*foo2-foo3)

#################################################
# Fp: base field
# Notice labels are positive values GREATER THAN 0

def getLabel(Fp):
    foo = False
    while foo == False:
    	  dummy = int(Fp.random_element())
	  if dummy != 0: foo = True
    ###return dummy
    return dummy%(P.order())

######################
# t: target polynomial

def getTarget(x, x5, x6, x7, p):
    foo = (x-x5)*(x-x6)*(x-x7)
    return foo%p
    
######################
# Fp: base field
# number of inputs = 4

def setInputData(E):
    mylist = []
    ###Fp = E.base_field()
    p  = P.order()
    Fp = GF(p)
    mylist.append(1)
    for i in range(4):
    	foo = False
	while foo == False:
	      dummy = int(Fp.random_element())
	      if dummy != 0:
	      	 foo = True
		 mylist.append(dummy)
    return (mylist)

################
#

def printDictionary(D):
    for x in D:
    	print x,':',D[x]