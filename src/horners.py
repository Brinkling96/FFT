

def horners(input, coefficents):
    if len(coefficents) >1:
        return coefficents[0] + input*horners(input, coefficents[1:])
    else:
        return coefficents[0]



for i in range(5):
    print(horners(i,[1,3,2]))