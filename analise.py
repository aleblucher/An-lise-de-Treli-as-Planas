import numpy as np
import math
from inn import readMecSol

def nos_sistema(quantidade):
    i = 0
    list = []
    for i in range (quantidade):
        x = float(input("Escreva a coordena x desse ponto: "))
        y = float(input("Escreva a coordena y desse ponto: "))
        list.append([x,y])
    return list

a = nos_sistema(3)
print(a)

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

maks = []

maks.append(elemento(a,1,2,0.0002,210000000000,0.4))
maks.append(elemento(a,2,3,0.0002,210000000000,0.3))
maks.append(elemento(a,3,1,0.0002,210000000000,0.5))

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

    print("maks", maks)
    #print("maks0", maks[0])
    #print("maks00", maks[0][0][(0,0)])
    print("big", big_matrix)

    return big_matrix


matriz_nos = make_Kg(maks, 3)


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

def findStrains(elast, length, lista_deslocamentos, no1, no2):
    ones = np.array([-1, 1])
    lowerUKnot = (lista_deslocamentos[((no1*2)-2)] + lista_deslocamentos[((no1*2)-1)])/2
    upperUKnot = (lista_deslocamentos[((no2*2)-2)] + lista_deslocamentos[((no2*2)-1)])/2

    lista_deslocs = np.array([lowerUKnot, upperUKnot])

    strain = (elast/length) * np.dot(ones, lista_deslocs)


def findStress(length, lista_deslocamentos):
    ones = np.array([-1, 1])
    lowerUKnot = (lista_deslocamentos[((no1*2)-2)] + lista_deslocamentos[((no1*2)-1)])/2
    upperUKnot = (lista_deslocamentos[((no2*2)-2)] + lista_deslocamentos[((no2*2)-1)])/2

    lista_deslocs = np.array([lowerUKnot, upperUKnot])

    stress = (1/length) * np.dot(ones, lista_deslocs)


matrize_nos= [[1.59, -0.4, -0.54], [-0.4, 1.7, 0.4], [-0.54, 0.4, 0.54]]
lista_u = gauss_method(matrize_nos, [0,150,-100])
resultInput = readMecSol("input.txt")


livres = np.zeros(2*int(resultInput['COORDINATES'][0][0]))
for i in range(1,len(resultInput['BCNODES'])):
    livres[int(((resultInput['BCNODES'][i][0]-1)*2) + (resultInput['BCNODES'][i][1]-1))]=-1
print(livres[0])
j=0
for i in range(len(livres)):
    if livres[i]==0:
        livres[i]=lista_u[j]
        j+=1
    elif livres[i]<0:
        livres[i]=0

print(livres)

print(np.transpose(livres)*matriz_nos)
