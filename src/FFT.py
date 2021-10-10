import cmath
import math
import numpy

def FFT(coefficents):
    n = len(coefficents)
    
    if n == 1: #base
        return coefficents
    else: #recursive
        wn= cmath.rect(1,(2*math.pi/n))
        w = 1

        temp = splitList(coefficents)
        coff_even = temp[0]
        coff_odd = temp[1]
        a = FFT(coff_even)
        b = FFT(coff_odd)

        y = []
        while len(y) < n:
            y.append(complex(0,0))
        x = int(n/2)
        for k in range(x):
            y[k] = a[k] + w*b[k]
            y[k+x] = a[k] - w*b[k]
            w = w*wn
        return y



    return 


def padListLenToPower2(coefficents):
    
    x = 1
    while x < len(coefficents):
        x = x*2
    
    while len(coefficents) < x:
        coefficents.append(0)

    return coefficents

def splitList(coefficents):
    even = []
    odd =[]
    for i in range(len(coefficents)):
        if i%2 == 1:
            odd.append(coefficents[i])
        else:
            even.append(coefficents[i])
    
    return [even,odd]