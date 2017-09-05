######################################
# 

def verif(E, CRS, PoC, P):
    Fp = E.base_field()
    #print 'Fp',Fp

    VP   = PoC['V*P']
    WQ   = PoC['W*Q']

    n    = VP.order()
    p    = Fp.order()
    k    = GF(n)(p).multiplicative_order()

    ans1 = VP.tate_pairing(WQ, n, k)

    YP   = PoC['Y*P']
    Q    = CRS['Q']

    ans2 = YP.tate_pairing(Q, n, k)

    TP   = CRS['t(s)P']
    HQ   = PoC['h(s)*Q']

    ans3 = TP.tate_pairing(HQ, n, k)

    print 'ans1:',ans1
    print 'ans2*ans3:',ans2*ans3

    print 'ans1/ans2:',ans1/ans2
    print 'ans3:',ans3

    TQ   = CRS['t(s)Q']
    HP   = PoC['h(s)*P']

    ans3 = HP.tate_pairing(TQ, n, k)

    print 'ans1:',ans1
    print 'ans2*ans3:',ans2*ans3

    exit()

    #V = int(PoC['V*P']%Fp.order())
    #print 'V:',V, type(V)
    #W = int(PoC['W*Q']%Fp.order())
    #print 'W:',W, type(W)
    #VW = (V*W)%Fp.order()
    #print 'V*W:', VW

    #Y = PoC['Y*P']%Fp.order()
    #print 'Y:',Y, type(W)

    #T = CRS['t(s)P']%Fp.order()
    #print 'T:',T, type(T)
    #H = PoC['h(s)*Q']%Fp.order()
    #print 'H:',H, type(H)
    #YTH = (Y+T*H)%Fp.order()

    #print 'Y+T*H:',YTH
   
    #print '##############################################'

    VWP  = VW*CRS['P']
    n    = VWP.order()
    p    = Fp.order()
    print 'VWP.order():',n
    print 'Fp.order() :',p
    k    = GF(n)(p).multiplicative_order()
    e    = Integer((p^k-1)/n)
    ans1 = VWP.tate_pairing(CRS['Q'], n, k)
    print 'ans1:',ans1, ans1^n
    ans1a = VWP._miller_(CRS['Q'],n)
    VP    = (V%Fp.order())*CRS['P']
    WQ    = (W%Fp.order())*CRS['Q']
    ans1b = VP._miller_(WQ,n)
    print 'ans1(miller fn):', ans1a, ans1b
    print 'e:',e
    print 'ans1^e(tate):',ans1a^e, ans1b^(e) #ans12b^e

    S = E.random_point()
    S = EK(S)
    Q = CRS['Q']
    dummy1 = VP._miller_(WQ+S,n)
    dummy2 = VP._miller_(S,n)
    print 'S:',S
    print 'dummy1:',dummy1
    print 'dummy2:',dummy2    
    print 'result:',(dummy1/dummy2)^e
    dummy1 = VWP._miller_(Q+S,n)
    dummy2 = VWP._miller_(S,n)
    print 'dummy1:',dummy1
    print 'dummy2:',dummy2    
    print 'result:',(dummy1/dummy2)^e

    print 'VWP:',VWP
    print 'VWP:',Integer(PoC['V*P']*PoC['W*Q'])*CRS['P']
    print 'VWP:',(Integer(PoC['V*P']*PoC['W*Q'])%Fp.order())*CRS['P']
    print 'VWP:',Integer(PoC['V*P']*PoC['W*Q'])
    #for i in range (Integer(PoC['V*P']*PoC['W*Q'])): print i*CRS['P']

    YTHP = YTH*CRS['P']
    ans2 = YTHP.tate_pairing(CRS['Q'], n, k)
    print 'ans2:',ans2, ans2^n, YTHP._miller_(CRS['Q'],n)^e

    VP   = (V%Fp.order())*CRS['P']
    WQ   = (W%Fp.order())*CRS['Q']
    ans3 = VP.tate_pairing(WQ, n, k)
    ###ans31 = WQ.tate_pairing(VP, n, k)
    print 'ans3:',ans3, ans3^n, VP in G1_EK, WQ in G2, VP._miller_(WQ,n)^e
    print 'ans31:',ans31
    print 'ans3/ans31:',ans3/ans31

    P     = CRS['P']
    Q     = CRS['Q']
    ans4  = P.tate_pairing(CRS['Q'], n, k)
    ans41 = pow(ans4,(VW)%p)
    ans42 = pow(ans4,VW)
    ###print 'ans4:',ans4, ans4^n
    print 'ans41:',ans41, ans41^n
    print 'ans42:',ans42, ans42^n 

    ans51 = P.tate_pairing(W*Q, n, k)
    ans52 = (W*CRS['P']).tate_pairing(CRS['Q'], n, k)
    print 'ans51:',ans51,'ans52:',ans52

    ans51 = P.tate_pairing(V*Q, n, k)
    ans52 = (V*CRS['P']).tate_pairing(CRS['Q'], n, k)
    print 'ans51:',ans51,'ans52:',ans52

    ans51 = (V*W*P).tate_pairing(Q, n, k)
    ans52 = P.tate_pairing(V*W*Q, n, k)
    print 'ans51:',ans51,'ans52:',ans52

    ans51 = (((V*W)%Fp.order())*P).tate_pairing(Q, n, k)
    ans52 = P.tate_pairing(((V*W)%Fp.order())*Q, n, k)
    print 'ans51:',ans51,'ans52:',ans52


    pause

    print 'P:',P, type(P), CRS['P'], type(CRS['P'])
    print 'V:',V, PoC['V*P']
    print '130*P:',130*E(P)
    print '12368*P:',12368*E(P)
    V1 = V*CRS['P']
    print 'V1:',V1
    V2 = int((PoC['V*P']%Fp.order()))*CRS['P']
    print 'V2:',V2
    V3 = PoC['V*P']*CRS['P']
    print 'V3:',V3
    V4 = int(PoC['V*P'])*CRS['P']
    print 'V4:',V4
    pause
    ####W = int((PoC['W*Q']%Fp.order()))#*CRS['Q']
    ####print 'V:',V, V in G1_EK
    ####print 'W:',W, W in G2
    ###V = int(PoC['V*P'])*CRS['P']
    ####W = int(PoC['W*Q'])#*CRS['Q']
    ###print 'V:',V, V in G1_EK
    ###print 'W:',W, W in G2
    ###print 'type(CRS[P]):',type(CRS['P']), CRS['P'], CRS['P'].order()
    ###print '-----', 134*CRS['P']
    print '@@@@@', 53*CRS['P']
    print '%%%%%', (134+pow(53,4))*CRS['P']
    ###for i in range(50):
    ###	print 'i:',i,2+i*53,(2+i*53*53)*CRS['P']
    pause
    n = V.order()
    k = GF(n)(Fp.order()).multiplicative_order()
    foo1 = V.tate_pairing(W, n, k)
    print 'n:',n
    print 'k:',k
    print 'e(V*P,W*Q):',foo1
    
    T = CRS['t(s)P']*CRS['P']
    H = PoC['h(s)*Q']*CRS['Q']
    n = T.order()
    k = GF(n)(Fp.order()).multiplicative_order()
    foo2 = T.tate_pairing(H, n, k)
    print 'n:',n
    print 'k:',k
    print 'e(T*P,H*Q):',foo2

    H = PoC['h(s)*P']*CRS['P']
    T = CRS['t(s)P']*CRS['Q']
    n = H.order()
    k = GF(n)(Fp.order()).multiplicative_order()
    foo2 = H.tate_pairing(T, n, k)
    print 'n:',n
    print 'k:',k
    print 'e(H*P,T*Q):',foo2

    Y = PoC['Y*P']*CRS['P']
    Q = CRS['Q']
    n = Y.order()
    k = GF(n)(Fp.order()).multiplicative_order()
    foo3 = Y.tate_pairing(Q, n, k)
    print 'n:',n
    print 'k:',k
    print 'e(Y*P,Q):',foo3

    P = CRS['P']
    n = P.order()
    k = GF(n)(Fp.order()).multiplicative_order()
    foo4 = P.tate_pairing(Q, n, k)
    print 'n:',n
    print 'k:',k
    print 'e(P,Q):',foo4
    V = int(PoC['V*P'])
    W = int(PoC['W*Q'])
    Y = int(PoC['Y*P'])
    T = int(CRS['t(s)P'])
    H = int(PoC['h(s)*P'])
    print 'V:',V,'W:',W,'Y:',Y,'T:',T,'H:',T
    dummy1 = V*W
    dummy2 = Y+T*H
    print 'e(P,Q)^(dummy1):',pow(foo4,dummy1)
    print 'e(P,Q)^(dummy2):',pow(foo4,dummy2)
    print 'dummy1:',dummy1,'dummy2:',dummy2
    dummy1 %= Fp.order()
    dummy2 %= Fp.order()
    print 'e(P,Q)^(dummy1):',pow(foo4,dummy1)
    print 'e(P,Q)^(dummy2):',pow(foo4,dummy2)
    print 'dummy1:',dummy1,'dummy2:',dummy2
    print 'e(P,Q)^(V*W):',pow(foo4,PoC['V*P']*PoC['W*Q'])
    print 'e(P,Q)^(Y+T*H):',(pow(foo4,PoC['Y*P']+PoC['h(s)*P']*CRS['t(s)P']))
    print 'e(P,Q)^(Y)  :',pow(foo4,PoC['Y*P'])
    print 'e(P,Q)^(T*H):',pow(foo4,PoC['h(s)*P']*CRS['t(s)P'])

    print 'e(P,Q)^(V*W-Y)',pow(foo4,(PoC['V*P']*PoC['W*Q']-PoC['Y*P']))
    
    return True  

######################################
# 

def verif_old(E, CRS, PoC, P):
    Fp = E.base_field()
    print 'Fp',Fp

    V = int(PoC['V*P']%Fp.order())
    print 'V:',V, type(V)
    W = int(PoC['W*Q']%Fp.order())
    print 'W:',W, type(W)
    VW = (V*W)%Fp.order()
    print 'V*W:', VW

    Y = PoC['Y*P']%Fp.order()
    print 'Y:',Y, type(W)

    T = CRS['t(s)P']%Fp.order()
    print 'T:',T, type(T)
    H = PoC['h(s)*Q']%Fp.order()
    print 'H:',H, type(H)
    YTH = (Y+T*H)%Fp.order()

    print 'Y+T*H:',YTH
   
    print '##############################################'

    VWP  = VW*CRS['P']
    n    = VWP.order()
    p    = Fp.order()
    print 'VWP.order():',n
    print 'Fp.order() :',p
    k    = GF(n)(p).multiplicative_order()
    e    = Integer((p^k-1)/n)
    ans1 = VWP.tate_pairing(CRS['Q'], n, k)
    print 'ans1:',ans1, ans1^n
    ans1a = VWP._miller_(CRS['Q'],n)
    VP    = (V%Fp.order())*CRS['P']
    WQ    = (W%Fp.order())*CRS['Q']
    ans1b = VP._miller_(WQ,n)
    print 'ans1(miller fn):', ans1a, ans1b
    print 'e:',e
    print 'ans1^e(tate):',ans1a^e, ans1b^(e) #ans12b^e

    S = E.random_point()
    S = EK(S)
    Q = CRS['Q']
    dummy1 = VP._miller_(WQ+S,n)
    dummy2 = VP._miller_(S,n)
    print 'S:',S
    print 'dummy1:',dummy1
    print 'dummy2:',dummy2    
    print 'result:',(dummy1/dummy2)^e
    dummy1 = VWP._miller_(Q+S,n)
    dummy2 = VWP._miller_(S,n)
    print 'dummy1:',dummy1
    print 'dummy2:',dummy2    
    print 'result:',(dummy1/dummy2)^e

    print 'VWP:',VWP
    print 'VWP:',Integer(PoC['V*P']*PoC['W*Q'])*CRS['P']
    print 'VWP:',(Integer(PoC['V*P']*PoC['W*Q'])%Fp.order())*CRS['P']
    print 'VWP:',Integer(PoC['V*P']*PoC['W*Q'])
    #for i in range (Integer(PoC['V*P']*PoC['W*Q'])): print i*CRS['P']

    YTHP = YTH*CRS['P']
    ans2 = YTHP.tate_pairing(CRS['Q'], n, k)
    print 'ans2:',ans2, ans2^n, YTHP._miller_(CRS['Q'],n)^e

    VP   = (V%Fp.order())*CRS['P']
    WQ   = (W%Fp.order())*CRS['Q']
    ans3 = VP.tate_pairing(WQ, n, k)
    ###ans31 = WQ.tate_pairing(VP, n, k)
    print 'ans3:',ans3, ans3^n, VP in G1_EK, WQ in G2, VP._miller_(WQ,n)^e
    print 'ans31:',ans31
    print 'ans3/ans31:',ans3/ans31

    P     = CRS['P']
    Q     = CRS['Q']
    ans4  = P.tate_pairing(CRS['Q'], n, k)
    ans41 = pow(ans4,(VW)%p)
    ans42 = pow(ans4,VW)
    ###print 'ans4:',ans4, ans4^n
    print 'ans41:',ans41, ans41^n
    print 'ans42:',ans42, ans42^n 

    ans51 = P.tate_pairing(W*Q, n, k)
    ans52 = (W*CRS['P']).tate_pairing(CRS['Q'], n, k)
    print 'ans51:',ans51,'ans52:',ans52

    ans51 = P.tate_pairing(V*Q, n, k)
    ans52 = (V*CRS['P']).tate_pairing(CRS['Q'], n, k)
    print 'ans51:',ans51,'ans52:',ans52

    ans51 = (V*W*P).tate_pairing(Q, n, k)
    ans52 = P.tate_pairing(V*W*Q, n, k)
    print 'ans51:',ans51,'ans52:',ans52

    ans51 = (((V*W)%Fp.order())*P).tate_pairing(Q, n, k)
    ans52 = P.tate_pairing(((V*W)%Fp.order())*Q, n, k)
    print 'ans51:',ans51,'ans52:',ans52


    pause

    print 'P:',P, type(P), CRS['P'], type(CRS['P'])
    print 'V:',V, PoC['V*P']
    print '130*P:',130*E(P)
    print '12368*P:',12368*E(P)
    V1 = V*CRS['P']
    print 'V1:',V1
    V2 = int((PoC['V*P']%Fp.order()))*CRS['P']
    print 'V2:',V2
    V3 = PoC['V*P']*CRS['P']
    print 'V3:',V3
    V4 = int(PoC['V*P'])*CRS['P']
    print 'V4:',V4
    pause
    ####W = int((PoC['W*Q']%Fp.order()))#*CRS['Q']
    ####print 'V:',V, V in G1_EK
    ####print 'W:',W, W in G2
    ###V = int(PoC['V*P'])*CRS['P']
    ####W = int(PoC['W*Q'])#*CRS['Q']
    ###print 'V:',V, V in G1_EK
    ###print 'W:',W, W in G2
    ###print 'type(CRS[P]):',type(CRS['P']), CRS['P'], CRS['P'].order()
    ###print '-----', 134*CRS['P']
    print '@@@@@', 53*CRS['P']
    print '%%%%%', (134+pow(53,4))*CRS['P']
    ###for i in range(50):
    ###	print 'i:',i,2+i*53,(2+i*53*53)*CRS['P']
    pause
    n = V.order()
    k = GF(n)(Fp.order()).multiplicative_order()
    foo1 = V.tate_pairing(W, n, k)
    print 'n:',n
    print 'k:',k
    print 'e(V*P,W*Q):',foo1
    
    T = CRS['t(s)P']*CRS['P']
    H = PoC['h(s)*Q']*CRS['Q']
    n = T.order()
    k = GF(n)(Fp.order()).multiplicative_order()
    foo2 = T.tate_pairing(H, n, k)
    print 'n:',n
    print 'k:',k
    print 'e(T*P,H*Q):',foo2

    H = PoC['h(s)*P']*CRS['P']
    T = CRS['t(s)P']*CRS['Q']
    n = H.order()
    k = GF(n)(Fp.order()).multiplicative_order()
    foo2 = H.tate_pairing(T, n, k)
    print 'n:',n
    print 'k:',k
    print 'e(H*P,T*Q):',foo2

    Y = PoC['Y*P']*CRS['P']
    Q = CRS['Q']
    n = Y.order()
    k = GF(n)(Fp.order()).multiplicative_order()
    foo3 = Y.tate_pairing(Q, n, k)
    print 'n:',n
    print 'k:',k
    print 'e(Y*P,Q):',foo3

    P = CRS['P']
    n = P.order()
    k = GF(n)(Fp.order()).multiplicative_order()
    foo4 = P.tate_pairing(Q, n, k)
    print 'n:',n
    print 'k:',k
    print 'e(P,Q):',foo4
    V = int(PoC['V*P'])
    W = int(PoC['W*Q'])
    Y = int(PoC['Y*P'])
    T = int(CRS['t(s)P'])
    H = int(PoC['h(s)*P'])
    print 'V:',V,'W:',W,'Y:',Y,'T:',T,'H:',T
    dummy1 = V*W
    dummy2 = Y+T*H
    print 'e(P,Q)^(dummy1):',pow(foo4,dummy1)
    print 'e(P,Q)^(dummy2):',pow(foo4,dummy2)
    print 'dummy1:',dummy1,'dummy2:',dummy2
    dummy1 %= Fp.order()
    dummy2 %= Fp.order()
    print 'e(P,Q)^(dummy1):',pow(foo4,dummy1)
    print 'e(P,Q)^(dummy2):',pow(foo4,dummy2)
    print 'dummy1:',dummy1,'dummy2:',dummy2
    print 'e(P,Q)^(V*W):',pow(foo4,PoC['V*P']*PoC['W*Q'])
    print 'e(P,Q)^(Y+T*H):',(pow(foo4,PoC['Y*P']+PoC['h(s)*P']*CRS['t(s)P']))
    print 'e(P,Q)^(Y)  :',pow(foo4,PoC['Y*P'])
    print 'e(P,Q)^(T*H):',pow(foo4,PoC['h(s)*P']*CRS['t(s)P'])

    print 'e(P,Q)^(V*W-Y)',pow(foo4,(PoC['V*P']*PoC['W*Q']-PoC['Y*P']))
    
    return True  

####################################
# c: fInput
# h: polynomial s.t. p(x) = h(x)t(x)

def generatePOC(E, x5, x6, x7, c, c5, c6, c7, CRS):
    Fp = E.base_field()
    R  = PolynomialRing(Fp, 'x')
    p  = Fp.order()
    xlist = []
    ylist = []
    ###print 'R:',R
    for i in range(4):
    	foo = False
	while foo == False:
		x_i = int(Fp.random_element())
		if fHardCoded == true:
		   if   i==0: x_i = 8
		   elif i==1: x_i = 1
		   elif i==2: x_i = 2
		   else     : x_i = 3
		if (x_i!=0) and (x_i!=x5) and (x_i!=x6) and (x_i!=x7): foo = True
	dummy1 = getTarget(x_i, x5, x6, x7, p)
	dummy1 = inverse_mod(dummy1, p)
	dummy2 = getP(x_i, x5, x6, x7, p, c, c5, c6, c7)
	y_i  = dummy2 * dummy1
	y_i %= Fp.order()
	xlist.append(x_i)
	ylist.append(y_i) 
    h = R.lagrange_polynomial([(xlist[i], ylist[i]) for i in range(4)])
    h_coeff = h.coefficients()
    print '   h(x)=',h
    ###print '   h_i:',h_coeff[0], ' ', h_coeff[1]
    PoCK = {}
    PoCK['V*P']    = c[0]*CRS['V0(s)P']+c[2]*CRS['V2(s)P']+c5*CRS['V5(s)P']
    PoCK['W*Q']    = c[1]*CRS['W1(s)Q']+c[3]*CRS['W3(s)Q']+c[4]*CRS['W4(s)Q']+c6*CRS['W6(s)Q']
    PoCK['Y*P']    = c5*CRS['Y5(s)P']+c6*CRS['Y6(s)P']+c7*CRS['Y7(s)P']
    PoCK['h(s)*P'] = int(h_coeff[0])*CRS['P']+int(h_coeff[1])*CRS['s^1*P']   
    PoCK['h(s)*Q'] = int(h_coeff[0])*CRS['Q']+int(h_coeff[1])*CRS['s^1*Q']   

    print '   PUBLIC'
    print '      PoC:',PoCK

    s    = secret
    foo1 = c[0]*getV0(s, x5, x6, x7, Fp.order())+c[2]*getV2(s, x5, x6, x7, Fp.order())+c5*getV5(s, x5, x6, x7, Fp.order())
    foo2 = c[1]*getW1(s, x5, x6, x7, Fp.order())+c[3]*getW3(s, x5, x6, x7, Fp.order())+c[4]*getW4(s, x5, x6, x7, Fp.order())+c6*getW6(s, x5, x6, x7, Fp.order())
    foo3 = c5*getY5(s, x5, x6, x7, Fp.order())+c6*getY6(s, x5, x6, x7, Fp.order())+c7*getY7(s, x5, x6, x7, Fp.order())
    foo4 = getTarget(s, x5, x6, x7, Fp.order())
    foo5 = int(h_coeff[0])+int(h_coeff[1])*s  
    print '   Sanity check:'
    print '        v(s)          =',foo1
    print '        w(s)          =',foo2
    print '        y(s)          =',foo3
    print '        t(s)          =',foo4
    print '        h(s)          =',foo5
    print '        v(s)*w(s)     =',(foo1*foo2),(foo1*foo2)%Fp.order()
    print '        y(s)+t(s)*h(s)=',(foo3+foo4*foo5),(foo3+foo4*foo5)%Fp.order()
   

    return PoCK

#############################################
# c5, c6, c7: values of the circuit wires 

def evaluateQAP(x5, x6, x7, E, c):
    Fp = E.base_field()
    p  = Fp.order()
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
    Fp = E.base_field()
    PubK = {}

    PubK['P']      = P
    PubK['Q']      = Q
    PubK['V0(s)P'] = getV0(s, x5, x6, x7, Fp.order())*P
    PubK['V2(s)P'] = getV2(s, x5, x6, x7, Fp.order())*P
    PubK['V5(s)P'] = getV5(s, x5, x6, x7, Fp.order())*P
    PubK['W1(s)Q'] = getW1(s, x5, x6, x7, Fp.order())*Q
    PubK['W3(s)Q'] = getW3(s, x5, x6, x7, Fp.order())*Q
    PubK['W4(s)Q'] = getW4(s, x5, x6, x7, Fp.order())*Q
    PubK['W6(s)Q'] = getW6(s, x5, x6, x7, Fp.order())*Q
    PubK['Y5(s)P'] = getY5(s, x5, x6, x7, Fp.order())*P
    PubK['Y6(s)P'] = getY6(s, x5, x6, x7, Fp.order())*P
    PubK['Y7(s)P'] = getY7(s, x5, x6, x7, Fp.order())*P
    PubK['s^1*P']  = s*P
    PubK['s^2*P']  = pow(s,2)*P
    PubK['s^3*P']  = pow(s,3)*P
    PubK['s^1*Q']  = s*Q
    PubK['s^2*Q']  = pow(s,2)*Q
    PubK['s^3*Q']  = pow(s,3)*Q
    PubK['t(s)P']  = getTarget(s, x5, x6, x7, Fp.order())*P
    PubK['t(s)Q']  = getTarget(s, x5, x6, x7, Fp.order())*Q

    ####
    #PubK['V0(s)P'] = getV0(s, x5, x6, x7, Fp.order())
    #PubK['V2(s)P'] = getV2(s, x5, x6, x7, Fp.order())
    #PubK['V5(s)P'] = getV5(s, x5, x6, x7, Fp.order())
    #PubK['W1(s)Q'] = getW1(s, x5, x6, x7, Fp.order())
    #PubK['W3(s)Q'] = getW3(s, x5, x6, x7, Fp.order())
    #PubK['W4(s)Q'] = getW4(s, x5, x6, x7, Fp.order())
    #PubK['W6(s)Q'] = getW6(s, x5, x6, x7, Fp.order())
    #PubK['Y5(s)P'] = getY5(s, x5, x6, x7, Fp.order())
    #PubK['Y6(s)P'] = getY6(s, x5, x6, x7, Fp.order())
    #PubK['Y7(s)P'] = getY7(s, x5, x6, x7, Fp.order())
    #PubK['s^1*P']  = s
    #PubK['s^2*P']  = pow(s,2)
    #PubK['s^3*P']  = pow(s,3)
    #PubK['s^1*Q']  = s
    #PubK['s^2*Q']  = pow(s,2)
    #PubK['s^3*Q']  = pow(s,3)
    #PubK['t(s)P']  = getTarget(s, x5, x6, x7, Fp.order())
    ###

    print '      CRS:',PubK
    return (PubK)

###################################################
# Fp: base field
# (x5, x6, x7): set of multiplication gates labeled
#      	   	by the index of their output values

def generateQAP(E, s):
    Fp = E.base_field()
    x5 = getLabel(Fp)
    x6 = getLabel(Fp)
    x7 = getLabel(Fp)

    if fHardCoded == True:
       x5 = 6
       x6 = 7
       x7 = 4

    t = getTarget(s, x5, x6, x7, Fp.order())

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
    return dummy

######################
# t: target polynomial

def getTarget(x, x5, x6, x7, p):
    foo = (x-x5)*(x-x6)*(x-x7)
    return foo%p

################
# Fp: base field

def setInputData(E):
    mylist = []
    Fp = E.base_field()
    mylist.append(1)
    for i in range(4):
    	foo = False
	while foo == False:
	      dummy = int(Fp.random_element())
	      if dummy != 0:
	      	 foo = True
		 mylist.append(dummy)
    if fHardCoded == True:
       mylist[:] = []
       mylist = [1,2,6,5,3]
       mylist = [1,2,2,7,3]
    return (mylist)