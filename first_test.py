import numpy as np


def nos_sistema(quantidade):
    i = 0
    list = []
    for i in range (quantidade):
        x = int(input("Escreva a coordena x desse ponto: "))
        y = int(input("Escreva a coordena y desse ponto: "))
        list.append([x,y])
    return list

a = nos_sistema(2)
print(a)

def elemento(lista_de_nos, no1, no2):
    x1 = lista_de_nos[no1][0]
    y1 = lista_de_nos[no1][1]
    x2 = lista_de_nos[no2][0]
    y2 = lista_de_nos[no2][1]
    dist = (x2 - x1)**2 + (y2-y1)**2
    cos = (x2 - x1)/np.sqrt(dist)
    sen = (y2 - y1)/np.sqrt(dist)
    return sen, cos

print(elemento(a, 0, 1))

def making_matrix(sen, cos):
    mak = np.matrix([[cos**2, cos*sen, -(cos**2), -(cos*sen)],[cos*sen, sen**2, -(cos*sen), -(sen**2)],[-(cos**2), -(cos*sen), cos**2, cos*sen],[-(cos*sen), -(sen**2), cos*sen, sen**2]])
    return mak

x,y = elemento(a,0,1)
print(making_matrix(x,y))
