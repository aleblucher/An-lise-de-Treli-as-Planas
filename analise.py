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

def calc_iteration(b1, x2, x3, a12, a13, a11):
    x1 = (b1-a12*x2-a13*x3)/a11
    return x1


def gauss_method(matriz_nos, forcas, graus_de_liberdade):
    lista_deslocs = [0] * graus_de_liberdade
    for i in range(graus_de_liberdade):
        #calc_iteration(forcas[i], )               #fazer entrada de forças ser uma lista
