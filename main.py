from math import *

from numpy import linspace
from scipy.misc import derivative


def f(x):
    return x**2 + log(x)

def diff_func(x):
    return 2*x + 1/x


def finite_differences_table(x):
    table = []
    for i in range(len(x)):
        table.append([])
        for j in range(len(x)):
            table[i].append(0)

    for i in range(len(x)):
        table[i][0] = f(x[i])

    for j in range(1, len(x)):
        for i in range(len(x) - j):
            table[i][j] = table[i + 1][j - 1] - table[i][j - 1]

    return table



def firstNewton(nodes, differenceTable, t):
    result = 0
    for i in range(nodes):
        resultMult = 1
        for j in range(i - 1):
            resultMult *= (t - j)
        result += (differenceTable[i][0] / factorial(i + 1)) * resultMult
    return result


def secondNewton(nodes, differenceTable, t):
    result = 0
    for i in range(nodes):
        resultMult = 1
        for j in range(1, i):
            resultMult *= ((t + (j - 1)) / j)
        result += differenceTable[i][nodes - i - 1] * resultMult
    return result


def firstGaus(nodes, differenceTable, t):
    result = 0
    for i in range(nodes):
        resultMult = 1
        for j in range(i - 1):
            if (j % 2 != 0):
                resultMult *= (t - (j + 1) / 2)
            else:
                resultMult *= (t + (j - 1) / 2)
        result += (differenceTable[nodes // 2 - 1][0] * (1 / factorial(i + 1))) * resultMult
    return result


def secondGaus(nodes, differenceTable, t):
    result = 0
    for i in range(nodes):
        resultMult = 1
        for j in range(i - 1):
            if (j % 2 != 0):
                resultMult *= (t + (j + 1) / 2)
            else:
                resultMult *= (t - (j + 1) / 2)
        result += (differenceTable[nodes // 2 - 1][nodes - i - 1] * (1 / factorial(i + 1))) * resultMult
    return result


def Rn(nodes, x0, x_values, x_star):
    result = abs(derivative(diff_func, x0, n=nodes, order=nodes+1))
    result /= factorial(nodes+1)
    for x in x_values:
        result *= (x_star - x)
    return result


a = 0.4
b = 0.9
x_arr_stars = [0.42, 0.87, 0.67]

x_arr = linspace(a, b, num=10)
rnx_arr = []
rnx_values = []
difTable = finite_differences_table(x_arr)
h = x_arr[1] - x_arr[0]
for i in range(len(x_arr_stars)):
    distance = 1000000
    x0 = -1
    for j in range(len(x_arr)):
        if abs(x_arr[j] - x_arr_stars[i]) < distance:
            distance = x_arr[j] - x_arr_stars[i]
            x0 = x_arr[j]
    rnx_values.append(Rn(10, x0, x_arr, x_arr_stars[i]))
    t = abs(x_arr_stars[i] - x0) / h
    Lnx = 0
    if (x_arr_stars[i] <= x_arr[3]):
        Lnx = firstNewton(10, difTable, t)
    if (x_arr_stars[i] >= x_arr[len(x_arr)-4]):
        Lnx = secondNewton(10, difTable, t)
    if (x_arr_stars[i] > x_arr[3] and x_arr_stars[i] < x_arr[int(len(x_arr)/2)]):
        Lnx = firstGaus(10, difTable, t)
    if (x_arr_stars[i] < x_arr[len(x_arr)-4] and x_arr_stars[i] > x_arr[int(len(x_arr)/2)]):
        Lnx = secondGaus(10, difTable, t)

    rnx = Lnx - f(x_arr_stars[i])
    rnx_arr.append(rnx)

for i in rnx_arr:
    print(min(rnx_values) < i < max(rnx_values))