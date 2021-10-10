import numpy

def poly_mult(a_coefficents, b_coefficients):
    convolution = numpy.zeros(len(a_coefficents)+len(b_coefficients))
    for i,a in enumerate(a_coefficents):
        for j,b in enumerate(b_coefficients):
            convolution[i+j] += a*b
    return convolution

def point_vector_mult(a_point_vector, b_point_vector):
    if len(a_point_vector) != len(b_point_vector):
        return
    for i in range(len(a_point_vector)):
       a_point_vector[i][1] *= b_point_vector[i][1] 
    return a_point_vector

    
