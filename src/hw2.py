from DFT import *
from Lagrange import *

f_convolution = DFT([-1,0,1,6])

f_coefficent = Inverse_DFT(f_convolution)

print("DFT of f(x)")
for ele in f_convolution:
    print(str(int(round(ele.real,0))) + " + " +str(int(round(ele.imag,0)))+"j")

g_convolution = DFT([3,-4,0,2])

g_cofficent = Inverse_DFT(g_convolution)

print("\nDFT of g(x)")
for ele in g_convolution:
    print(str(int(round(ele.real,0))) + " + " +str(int(round(ele.imag,0))) +"j")

fg_convolution = []

for i in range(len(f_convolution)):
    fg_convolution.append(f_convolution[i]*g_convolution[i])

print("\nComponentwise product of f(x) * g(x)")
for ele in fg_convolution:
    print(str(int(round(ele.real,0))) + " + " +str(int(round(ele.imag,0))) +"j")    

fg_false_coff = Inverse_DFT(fg_convolution)

print("\nFalse fg Coefficent Vector")
for ele in fg_false_coff:
    print(str(int(round(ele.real,0))) + " + " +str(int(round(ele.imag,0))) +"j")   

h_coff = poly_mult([-1,0,1,6], [3,-4,0,2])

print(h_coff)

print(32768 % 1025)