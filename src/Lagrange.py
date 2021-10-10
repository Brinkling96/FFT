import numpy

def lagrange(point_vector):
    convolution = numpy.zeros(len(point_vector))
    for k in range(len(point_vector)):
        yk = float(point_vector[k][1])
        xk = point_vector[k][0]
        curr_poly =[1]
        curr_denominator =1
        for j in range(len(point_vector)):
            if j != k:
               curr_poly = poly_mult(curr_poly,[-point_vector[j][0],1])
               curr_denominator = curr_denominator*(xk - point_vector[j][0])
        scalar = yk/curr_denominator
        for i in range(len(point_vector)):
            convolution[i] += curr_poly[i]*scalar
    return convolution

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

#print(lagrange(point_vector_mult([[0,9],[1,10],[2,13],[3,18],[4,25]], [[0,1],[1,6],[2,15],[3,28],[4,45]])))

print(lagrange([[0,1024],[1,2048],[2,4096],[3,8192],[4,16384],[5,32768]]))