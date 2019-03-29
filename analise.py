import numpy as np
import math
from inn import readMecSol


##### RECEIVES THE INPUTS IN ORDER TO CREATE OUR ELEMENT'S KNOTS #######
##### RETURNS A LIST OF (X, Y) COORDINATES LISTS 

def nos_sistema(dataa):
    list = []
    for i in range(1,len(dataa['COORDINATES'])):
        x = float(dataa['COORDINATES'][i][1])
        y = float(dataa['COORDINATES'][i][2])
        list.append([x,y])
    return list


########### CREATE OUR ELEMENTS ##########
########### RETURNS A LIST WITH [ELEMENT_MATRIX, KNOT1, KNOT2]

def elemento(lista_de_nos, no1, no2, area, elast, comprimento):
    prefixo = elast * area/comprimento
    x1 = lista_de_nos[no1 - 1][0]
    y1 = lista_de_nos[no1 - 1][1]
    x2 = lista_de_nos[no2 - 1][0]
    y2 = lista_de_nos[no2 - 1][1]
    dist = (x2 - x1)**2 + (y2-y1)**2
    cos = (x2 - x1)/np.sqrt(dist)
    sen = (y2 - y1)/np.sqrt(dist)
    mak = np.matrix([[cos**2, cos*sen, -(cos**2), -(cos*sen)],[cos*sen, sen**2, -(cos*sen), -(sen**2)],[-(cos**2), -(cos*sen), cos**2, cos*sen],[-(cos*sen), -(sen**2), cos*sen, sen**2]])
    k_e = prefixo*mak
    return [k_e/(10**8), no1, no2]


######### RETURNS THE BIG MATRIX WITH OVERLAPPED ELEMENTS

def make_Kg(maks, n_nos):
    big_matrix = np.zeros((n_nos*2, n_nos*2))
    for mak in maks:
        first = (mak[1] * 2) - 1
        second = mak[1] * 2
        third = (mak[2] * 2) - 1
        fourth = mak[2] * 2
        if third > second:
            for line in range(first, fourth + 1):
                for column in range(first, fourth + 1):
                    big_matrix[line - 1][column - 1] += mak[0][(line - first, column - first)]
        else:
            for line in range(first, second + 1):
                for column in range(first, second + 1):
                    big_matrix[line - 1][column - 1] += mak[0][(line - first, column - first)]
                for column in range(third, fourth + 1):
                    big_matrix[line - 1][column - 1] -= mak[0][(line - first, column - first)]

            for line in range(third, fourth + 1):
                for column in range(first, second + 1):
                    big_matrix[line - 1][column - 1] -= mak[0][(line - first, column - first)]
                for column in range(third, fourth + 1):
                    big_matrix[line - 1][column - 1] += mak[0][(line - first, column - first)]
    return big_matrix


###### IMPLEMENTED THE GAUSS NUMERIC METHOD - RECEIVES TWO MATRIX TO FIND A THIRD ONE
###### RETURN THE U MATRIX

# | a b c |   | u1 |   | F1 |
# | d e f | x | u2 | = | F2 |
# | g h i |   | u3 |   | F3 |


def gauss_method(matriz_nos, forcas):

    lista_u = [0] * len(matriz_nos) # create the U array with zeros

    preview = lista_u[0] - 1 # start the first part of the comparisson as zero
    tolerance = 0.00001 # convergence parameter
    iterations  = 0 # start a counter of iteratios (how many times until reach the real value)
    lista_u[0] = 1

    while (np.abs(preview - lista_u[0]) > tolerance):
        iterations += 1
        preview = lista_u[0]

        for i in range(len(matriz_nos)):

            divisor = 0 #start our divisor (wich is matriz_nos[i][j])
            soma = 0

            for j in range(len(matriz_nos[i])):
                if i==j:
                    divisor = matriz_nos[i][j]
                else:
                    soma+=matriz_nos[i][j]*lista_u[j]

            if divisor == 0:
                lista_u[i] = 0
            else:
                lista_u[i] = (forcas[i]-soma)/divisor

    lista_u = np.array(lista_u)*(10**-8) #TODO: fix
    print(lista_u)
    return lista_u


##### 

def findStrains(elast, length, lista_deslocamentos, no1, no2):
    ones = np.array([-1, 1])
    lowerUKnot = (lista_deslocamentos[((no1*2)-2)] + lista_deslocamentos[((no1*2)-1)])/2
    upperUKnot = (lista_deslocamentos[((no2*2)-2)] + lista_deslocamentos[((no2*2)-1)])/2

    lista_deslocs = np.array([lowerUKnot, upperUKnot])

    strain = (elast/length) * np.dot(ones, lista_deslocs)
    return strain


def findStress(length, lista_deslocamentos, no1, no2):
    ones = np.array([-1, 1])
    lowerUKnot = (lista_deslocamentos[((no1*2)-2)] + lista_deslocamentos[((no1*2)-1)])/2
    upperUKnot = (lista_deslocamentos[((no2*2)-2)] + lista_deslocamentos[((no2*2)-1)])/2

    lista_deslocs = np.array([lowerUKnot, upperUKnot])

    stress = (1/length) * np.dot(ones, lista_deslocs)
    return stress

####### DEFININDO CONDICAO DE CONTORNO

def boundaryConditions(matrix1, matrix2):
    indexesToGet = []
    for i in range(len(matrix1)):
        if matrix1[i] != 0:
            indexesToGet.append(i)

    new_matrix = []
    
    for k in range(len(matrix2)):
        reducedLine = []
        if k in indexesToGet:
            for j in range(len(matrix2[i])):
                if j in indexesToGet:
                    reducedLine.append(matrix2[i][j])
        new_matrix.append(reducedLine)
    return np.array(new_matrix)



# lista_u = gauss_method(matrize_nos, [0,150,-100])

resultInput = readMecSol("input.txt")

forcas = resultInput
incidences = resultInput['INCIDENCES']
area = resultInput['GEOMETRIC_PROPERTIES']
materials = resultInput['MATERIALS']

a = nos_sistema(resultInput)

maks = []

for i in range(len(incidences)):
    maks.append(elemento(a, incidences[i][1], incidences[i][2], area[5], materials[i][0], materials[i][2]))

# maks.append(elemento(a,1,2,0.0002,210000000000,0.4))
# maks.append(elemento(a,2,3,0.0002,210000000000,0.3))
# maks.append(elemento(a,3,1,0.0002,210000000000,0.5))

# matriz_nos_elementos = make_Kg(maks, 3)
# print( "AAAAAAAA", matriz_nos)


# lista_u = gauss_method(matrize_nos, forcas)

# livres = np.zeros(2*int(resultInput['COORDINATES'][0][0]))

# for i in range(1,len(resultInput['BCNODES'])):
#     livres[int(((resultInput['BCNODES'][i][0]-1)*2) + (resultInput['BCNODES'][i][1]-1))]=-1
# print(livres[0])
# j=0
# for i in range(len(livres)):
#     if livres[i]==0:
#         livres[i]=lista_u[j]
#         j+=1
#     elif livres[i]<0:
#         livres[i]=0