reset()
import sys
from PIL import Image

load("PBCLib.sage")
load("QAP.sage")
load("tate_pairing.sage")

fDebug = False 

####################################
# E: EEC defined over the base field

E = setECC()

if fDebug==True: displayPointsOnECC(E)

########################################
# P: generator of the n-torsion group G1
# k: embedding degree of the subgroup G1
# G1: torsion subgroup over E

P, k, G1 = setTorsionFieldGenerator(E, False, -99)

########################################################
# FE: extension of the base field
# EK: elliptic curve defined over the extension field FE

FE, EK = setFieldExtension(E, k)

if fDebug == True: displayPointsOnECC(EK)

##################################################
# Q: generator of the n-torsion group G2
# G2: torsion subgroup of EK
# G1_EK: torsion subroup over EK (equivalent to G1)

Q, G2, G1_EK = setTorsionFieldGenerator(EK, True, EK(P))

#############################
# hashing in torsion subgroup

###foo = hashToG1('tintin et les picaros', E, G1, fFlag)
###if foo in G1: print 'H1:',foo, foo in G1
###foo = hashToG1('tintin et les cigares du pharaon', E, G1, fFlag)
###if foo in G1: print 'H1:',foo, foo in G1

###raw_input('\nPress enter to continue\n')

################################################
# Phase 1: pre-processing phase
#          build QAP(t(x), V(x), W(x), Y(x)
#          construct Common Reference String CRS

#imgAdd = '/home/elmagnifico/crypto/sage/circuit.png'
imgAdd = 'circuit.png'
img = Image.open(imgAdd)
img.show()

fInput = setInputData(E)
secret = 0
while( secret==0 ): secret = int(E.base_field().random_element())

x5, x6, x7 = generateQAP(E, secret)

EvalK, VerifK = generatePublicParams(secret, EK(P), Q, x5, x6, x7, E)

raw_input('\nPress enter to continue\n')

###########################################################
# Phase 2: circuit evaluation on input fInput and PoC,
#          worker evaluates QAP(fInput) to obtain
#          y<-QAP(fInput)
#          worker knows c5, c6, c7=y of the circuit's wires
#          worker constructs Proof-of-Correctness (PoC)

c5, c6, c7 = evaluateQAP(x5, x6, x7, E, fInput)

Proof = generatePOC(E, x5, x6, x7, fInput, c5, c6, c7, EvalK)

#############################################################
# Phase 3: verification
#          anyone with access to the CRS and PoC can use the
#          pairing function e to check that:
#          1. e(V(s)*P,W(s)*Q)/e(Y(s)*P,Q) = e(t(s)*P,h(s)*Q)

###print 'type(fInput):',type(fInput),fInput

fIO = fInput
fIO.append(c7)

###print 'fIO:',fIO

verif(EK, VerifK, Proof, fIO)