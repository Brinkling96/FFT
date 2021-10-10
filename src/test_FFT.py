import pytest
import cmath
import math

from naive_poly_mult import *
from Lagrange import *
from horners import *
from DFT import *
from FFT import *




def test_poly_mult():
    convolution  = poly_mult([9,-10,7,6], [-5,4,0,-2])
    assert convolution[0] == -45
    assert convolution[1] == 86
    assert convolution[2] == -75
    assert convolution[3] == -20
    assert convolution[4] == 44
    assert convolution[5] == -14
    assert convolution[6] == -12

def test_lagrange():
    coefficents = [-5,4,0,-2]
    convolution = lagrange([[0, horners(0, coefficents)],
                            [1, horners(1, coefficents)],
                            [2, horners(2, coefficents)],
                            [3, horners(3, coefficents)]])

    print(convolution)
    for i,p in enumerate (convolution):
        assert  coefficents[i] - 0.001 <= p <= coefficents[i] + 0.001

def test_DFT():
    coefficents = [9,-10,7,6]
    yk = DFT(coefficents)

    assert len(yk) == 4

    assert yk[0].real - 0.001 <= 12 <= yk[0].real + 0.001
    assert yk[0].imag - 0.001 <= 0 <= yk[0].imag + 0.001

    assert yk[1].real - 0.001 <= 2 <= yk[1].real + 0.001
    assert yk[1].imag - 0.001 <= -16 <= yk[1].imag + 0.001

    assert yk[2].real - 0.001 <= 20 <= yk[2].real + 0.001
    assert yk[2].imag - 0.001 <= 0 <= yk[2].imag + 0.001

    assert yk[3].real - 0.001 <= 2 <= yk[3].real + 0.001
    assert yk[3].imag - 0.001 <= 16 <= yk[3].imag + 0.001


def test_IDFT():
    convolution = [complex(12,0), complex(2,-16), complex(20,0), complex(2,16)]
    
    aj = Inverse_DFT(convolution)

    assert len(aj) == 4

    assert aj[0].real - 0.001 <= 9 <= aj[0].real + 0.001
    assert aj[0].imag - 0.001 <= 0 <= aj[0].imag + 0.001

    assert aj[1].real - 0.001 <= -10 <= aj[1].real + 0.001
    assert aj[1].imag - 0.001 <= 0 <= aj[1].imag + 0.001

    assert aj[2].real - 0.001 <= 7 <= aj[2].real + 0.001
    assert aj[2].imag - 0.001 <= 0 <= aj[2].imag + 0.001

    assert aj[3].real - 0.001 <= 6 <= aj[3].real + 0.001
    assert aj[3].imag - 0.001 <= 0 <= aj[3].imag + 0.001

def test_DFT_Wrap():
    a_coefficents = [9,-10,7,6,0]
    a_convolution = DFT(a_coefficents)
    a_test = Inverse_DFT(a_convolution)

    for i,p in enumerate (a_test):
        assert  a_coefficents[i] - 0.001 <= p.real <= a_coefficents[i] + 0.001



def test_DFT_Mult():
    a_coefficents = [9,-10,7,6]
    b_coefficents = [-5,4,0,-2]
    
    degA = len(a_coefficents)
    degB = len(b_coefficents)
    degC = degA + degB

    while len(a_coefficents) < degC:
        a_coefficents.append(0)

    while len(b_coefficents) < degC:
        b_coefficents.append(0)

    a_convolution = DFT(a_coefficents)
    b_convolution = DFT(b_coefficents)

    a_test = Inverse_DFT(a_convolution)
    
    c_convolution = []

    for i in range(len(a_convolution)):
        c_convolution.append(a_convolution[i]*b_convolution[i])

    c_coefficents = Inverse_DFT(c_convolution)

    assert len(c_coefficents) == 8

    assert c_coefficents[0].real - 0.001 <= -45 <= c_coefficents[0].real + 0.001
    assert c_coefficents[0].imag - 0.001 <= 0 <= c_coefficents[0].imag + 0.001

    assert c_coefficents[1].real - 0.001 <= 86 <= c_coefficents[1].real + 0.001
    assert c_coefficents[1].imag - 0.001 <= 0 <= c_coefficents[1].imag + 0.001

    assert c_coefficents[2].real - 0.001 <= -75 <= c_coefficents[2].real + 0.001
    assert c_coefficents[2].imag - 0.001 <= 0 <= c_coefficents[2].imag + 0.001

    assert c_coefficents[3].real - 0.001 <= -20 <= c_coefficents[3].real + 0.001
    assert c_coefficents[3].imag - 0.001 <= 0 <= c_coefficents[3].imag + 0.001

    assert c_coefficents[4].real - 0.001 <= 44 <= c_coefficents[4].real + 0.001
    assert c_coefficents[4].imag - 0.001 <= 0 <= c_coefficents[4].imag + 0.001

    assert c_coefficents[5].real - 0.001 <= -14 <= c_coefficents[5].real + 0.001
    assert c_coefficents[5].imag - 0.001 <= 0 <= c_coefficents[5].imag + 0.001

    assert c_coefficents[6].real - 0.001 <= -12 <= c_coefficents[6].real + 0.001
    assert c_coefficents[6].imag - 0.001 <= 0 <= c_coefficents[6].imag + 0.001


def test_FFT():
    coefficents = [9,-10,7,6]
    yk = FFT(coefficents)

    assert len(yk) == 4

    assert yk[0].real - 0.001 <= 12 <= yk[0].real + 0.001
    assert yk[0].imag - 0.001 <= 0 <= yk[0].imag + 0.001

    assert yk[1].real - 0.001 <= 2 <= yk[1].real + 0.001
    assert yk[1].imag - 0.001 <= -16 <= yk[1].imag + 0.001

    assert yk[2].real - 0.001 <= 20 <= yk[2].real + 0.001
    assert yk[2].imag - 0.001 <= 0 <= yk[2].imag + 0.001

    assert yk[3].real - 0.001 <= 2 <= yk[3].real + 0.001
    assert yk[3].imag - 0.001 <= 16 <= yk[3].imag + 0.001

def test_padCoefficentVector():
    coefficents = [9,-10,7]
    coefficents = padListLenToPower2(coefficents)

    assert len(coefficents) == 4

    coefficents.append(6)

    coefficents = padListLenToPower2(coefficents)

    assert len(coefficents) == 8

    coefficents = padListLenToPower2(coefficents)

    assert len(coefficents) == 8

