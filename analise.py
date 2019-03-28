import numpy as np

def nos_sistema(quantidade):
    i = 0
    list = []
    for i in range (quantidade):
        x = int(input("Escreva a coordena x desse ponto: "))
        y = int(input("Escreva a coordena y desse ponto: "))
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
    return [k_e, no1, no2]

maks = []

maks.append(elemento(a,1,2,1.2,0.5,10))
maks.append(elemento(a,2,3,1.6,0.7,10))
maks.append(elemento(a,3,1,1.6,0.7,10))

def make_Kg(maks, n_nos):
    big_matrix = np.zeros((n_nos*2, n_nos*2))
    for mak in maks:
        first = (mak[1] * 2) - 1
        second = mak[1] * 2
        third = (mak[2] * 2) - 1
        fourth = mak[2] * 2
        for line in range(first, fourth + 1):
            for column in range(first, fourth + 1):
                big_matrix[line - 1][column - 1] += mak[0][(line - first, column - first)]

    print("maks", maks)
    #print("maks0", maks[0])
    #print("maks00", maks[0][0][(0,0)])
    print("big", big_matrix)

    return big_matrix


matriz_nos = make_Kg(maks, 3)


def gauss_method(matriz_nos, forcas):

    lista_u = [0] * len(matriz_nos) # create the U array with zeros
    
    preview = lista_u[0] # start the first part of the comparisson as zero
    tolerance = 1 # convergence parameter
    iterations  = 0 # start a counter of iteratios (how many times until reach the real value)

    while (np.abs(preview - lista_u[0]) < tolerance):
        
        iterations += 1
        
        for i in range(len(matriz_nos)):
            
            divisor = 0 #start our divisor (wich is matriz_nos[i][j])
            soma = 0
            
            for j in range(len(matriz_nos[i])):
                if i==j:
                    divisor = matriz_nos[i][j]
                else:
                    soma+=matriz_nos[i][j]*lista_u[j]
            preview = lista_u[0]
            if divisor == 0:
                lista_u[i] = 0
            else:
                lista_u[i] = (forcas[i]-soma)/divisor
        print(iterations)
        print(preview)
        print(np.abs(preview - lista_u[0]))
        print("AAAAAA", lista_u[0])
    #print(lista_u)


gauss_method(matriz_nos, [1,2,3,4,5,6])