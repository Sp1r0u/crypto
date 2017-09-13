###############
# page 3-5 
def g(P, Q, R, IE):
    (xP, yP) = P.xy()
    (xQ, yQ) = Q.xy()
    (xR, yR) = R.xy()
    if (P!=IE) and (Q!=IE) and (P+Q!=IE):
       if (P!=Q): slope = (yP-yQ)/(xP-xQ)
       else:      slope = (3*xP^2)/(2*yP)
       y_intercept = yP - slope * xP 
       xPQ = slope^2 - xP - xQ
       ###print 'P:',P
       ###print 'Q:',Q
       ###print 'slope:',slope
       ###print 'y_intercept:',y_intercept
       ###print 'xPQ:',xPQ
       return (yR-slope*xR-y_intercept)/(xR-xPQ)
    elif (P==-Q) and (P!=IE):
       return (xR-xP)
    elif (P==IE) or (Q==IE):
       return 1
    else: print 'error in g'

################
# Algorithm 3.1 
def miller(n, P, Q, IE):
    n = bin(n)[2:]
    I = len(n)-1
    ###print 'n:',n 
    ###print 'I:',I
    T = P
    f = 1
    for i in range(I-1, -1, -1):
	f = (f^2)*g(T, T, Q, IE)
	T = 2*T
	if n[i]==1:
	   f = g(T, P, Q, IE)
	   T = T + P
    return f	   

################
#

def my_tate_pairing(E, P, Q):
    foo = False
    IE = E(0, 1, 0) #point at infinity on E
    while foo == False:
    	  R = E.random_element()
	  if (R!=P) and (R!=IE) and (R+Q!=P) and (R+Q!=IE):
	     foo = True
    ###print 'P :', P
    ###print 'Q :', Q
    ###print 'R :', R
    ###print 'IE:', IE

    Fp = E.base_field()
    Fp = Fp.base()
    p  = Fp.characteristic()
    n  = P.order()
    k  = GF(n)(p).multiplicative_order()

    foo = ( p^k - 1 )/n

    ###print 'foo:', foo

    ###dummy = miller(n, P, Q+R, IE)/miller(n, P, R, IE)
    dummy = miller(n, P, Q, IE)

    ###print 'dummy:', dummy

    return dummy^foo
    ###print 'ans3:', P.tate_pairing(Q, n, k)
    