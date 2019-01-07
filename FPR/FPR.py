# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import sys
import string
symbols = string.digits

##########################################

def DivisionMethod (number, new_base):
    # source: https://en.wikipedia.org/wiki/Division_algorithm
    sign = -1 if number < 0 else 1
    number *= sign
    ans = ''
    while number:
        ans += symbols[number % new_base]
        number //= new_base
    if sign == -1:
        ans += '-'
    return ans[::-1]

##########################################

def Map2Z (X, B, f):
    # Given an element X of Z_<k> returns X in Z
    return X * B ** - f

##########################################
    
def FPR (X, B, k, f, original_base=10):
    # source: https://codereview.stackexchange.com/questions/44553/converting-float-from-arbitrary-base-to-arbitrary-base
    # Given a rational number X in base original_base, return the fixed-point representation 
    #       of X in base B of precision f and the corresponding signed integer
    # B: new base
    # k: Length of the ouput representation
    # f: length of the fractional part
    # k-f: length of the integer part 
    #      (Note that if (k-f) is small, the integer part may be truncated and 
    #      the data may be mis-represented.)

    if k-f < 0 or f < 0:
        sys.exit ('error with output length k or resolution factor f')

    #from original_base
    integer, point, fractional = X.strip().partition('.')
    num = int (integer + fractional, original_base) * original_base ** -len (fractional)
    
    #to new_base
    f = len (fractional) if f is None else f
    s = DivisionMethod (int (round (num / B ** -f)), B)
    
    if num>0:
        s = '+' + s.zfill(k)
    else:
        s = s.zfill(k+1)
    
    return s #s[:-f] + '.' + s[-f:]

##########################################

def SignedInt (s, k, B):
    # Returns a signed-integer in Z_k of a given a fixed-point datatype s in Q_<k,f>
    t = 0
    for i in range (0, k):
        t += int(s[k-i]) * B ** (i)
    
    if int(s) < 0:
        return t * (-1)
    else:
        return t
 
##########################################
   
def Map2Field (X, p):
    # Maps a signed integer X in Z_k to the field Fp
    return X % p

##########################################

def Scale (X, f1, f2, B):
    # converts a signed integer X in Z-k with resolution f1 to
    # a signed integer in Z_k with resolution f2 by scaling-up if f2-f1 >= 0
    # or by scaling-down (i.e. truncate x) if f2-f2 < 0.
    m = f2 - f1
    if m > 0:
        return X * B ** m
    else:
        return int ( round (X * B ** m ) ) 
    
##########################################

def AddSI (X1_SI, X2_SI, f1, f2, B):
    # returns the sum of two signed integers and the length of the output 
    # fractional part (i.e. the resolution factor of the output signal fOut)   
    # the function tests the resolution f1 and f2, and rescale, if necessary, 
    # the signed integer with the best resolution factor.
    # fOut = min (f1, f2)  
    if f1 == f2:
        return X1_SI + X2_SI, f1
    elif f2>f1:
        return X1_SI + Scale (X2_SI, f2, f1, B), f1
    else:
        return Scale (X1_SI, f1,f2, B) + X2_SI, f2

##########################################

def AddF (X1_F, X2_F, p):
    # returns the sum of two elements from the same field Fp
    # Note: we assume that X1_F and X2_F belongs to Fp 
    #       and the signed integers X1_SI and X2_SI have the same resolution 
    #       factor f
    return (X1_F + X2_F) % p

######################################################
######################################################

# p: characteristic of the base field
#    it is a large (i.e. 128bits or more) prime number
p  = 231279409

# x: rational number in base10
X1 = '-3.14159'
X2 = '-2.21'

# B: new base, can be the same as the original base (i.e. base10)
B  = 10

# k: size of of the fixed-point data type (i.e. number of digits)
k  = 20

# f: resolution factor (i.e. the fixed-point data type has resolution B^(-f))
f1  = 5
f2  = 3

# X_FPR: fixed-point datatyoe in Q_<k,f>
# X_SI: signed integers in Z_k
# X_F: elements of the field F_p
# X_Z: output of the reverse mapping from Z_k -> Z
#      the absolute error err = |X_Z - X| 

X1_FPR = FPR (X1, B, k, f1)
X1_SI  = SignedInt (X1_FPR, k ,B)
X1_F   = Map2Field (X1_SI, p)
X1_Z   = Map2Z (X1_SI, B, f1)

print ('X1     :', X1)
print ('X1_FPR :', X1_FPR)
print ('X1_SI  :', X1_SI)
print ('X1_F   :', X1_F)
print ('X1_Z   :', X1_Z)
print ('')

###f1_tmp = f1 - 3

###X1_FPR_tmp = FPR (X1, B, k, f1_tmp)
###X1_SI_tmp  = SignedInt (X1_FPR_tmp, k ,B)
###X1_F_tmp   = Map2Field (X1_SI_tmp, p)
###X1_Z_tmp   = Map2Z (X1_SI_tmp, B, f1_tmp)

###print ('X1_FPR :', X1_FPR_tmp)
###print ('X1_SI  :', X1_SI_tmp)
###print ('X1_Z   :', X1_Z_tmp)
###print ('')

###X1_SI_scaled = Scale (X1_SI, f1 ,f1_tmp, B)
###X1_Z_scaled  = Map2Z (X1_SI_scaled, B, f1_tmp)

###print ('X1_SI_scaled :', X1_SI_scaled)
###print ('X1_Z_scaled  :', X1_Z_scaled)

###AbsErr = abs (float(X1) - X1_Z)
###print ('Abs. Err.(',AbsErr,') should be smaller than Trunc. Err.(', B ** -f1,')')

###AbsErr = abs (float(X1) - X1_Z_scaled)
###print ('Abs. Err.(',AbsErr,') should be smaller than Trunc. Err.(', B ** -f1_tmp,')') 

X2_FPR = FPR (X2, B, k, f2)
X2_SI  = SignedInt (X2_FPR, k ,B)
X2_F   = Map2Field (X2_SI, p)
X2_Z   = Map2Z (X2_SI, B, f2)

print ('X2     :', X2)
print ('X2_FPR :', X2_FPR)
print ('X2_SI  :', X2_SI)
print ('X2_F   :', X2_F)
print ('X2_Z   :', X2_Z)
print ('')

X3_SI, f3 = AddSI (X1_SI, X2_SI, f1, f2, B)
X3_Z      = Map2Z (X3_SI, B, f3)

if f1 == f2:
    X3_F =  AddF (X1_F, X2_F, p)
elif f2>f1:
    X2_SI_scaled = Scale     (X2_SI, f2, f1, B)
    X2_F_scaled  = Map2Field (X2_SI_scaled, p)
    X3_F         = AddF (X1_F, X2_F_scaled, p)
else:
    X1_SI_scaled = Scale     (X1_SI, f1, f2, B)
    X1_F_scaled  = Map2Field (X1_SI_scaled, p)
    X3_F         = AddF (X1_F_scaled, X2_F, p)

print ('X3_SI  :', X3_SI, f3)
print ('X3_F   :', X3_F)
print ('X3_Z   :', X3_Z)

Sum = float(X1) + float(X2)
AbsErr = abs ( Sum - X3_Z)
print ('Abs. Err.(',AbsErr,') should be smaller than Trunc. Err.(', B ** -f3,')') 

if f1>= f2:
    Sum_FPR = FPR (str (Sum), B, k, f2)
else:
    Sum_FPR = FPR (Sum, B, k, f1)
Sum_SI  = SignedInt (Sum_FPR, k ,B)
Sum_F   = Map2Field (Sum_SI, p)

print ('Sum_F :', Sum_F)

###X1_FPR_tmp = FPR (X1, B, k, f1_tmp)
###X1_SI_tmp  = SignedInt (X1_FPR_tmp, k ,B)
###X1_F_tmp   = Map2Field (X1_SI_tmp, p)
###X1_Z_tmp   = Map2Z (X1_SI_tmp, B, f1_tmp)

#print ('XB2 :', XB2)
#print ('XF2 :', XF2)
###print ('   :', xbar2 * new_base ** - precision)
###print ('   ', int(s) * (new_base ** precision))

    ###add = (x1f + x2f) % p
    ###print ('add :', add)
    ###print ('    :', (add) * new_base ** - precision)
    ###print ('    :', (add-p) * new_base ** - precision)
    
    ###mul = (x1f * x2f) % p
    ###print ('mul :', mul)
    ###print ('    :', (mul) * new_base ** - (2*precision))
    ###print ('    :', (mul-p) * new_base ** - (2*precision))


