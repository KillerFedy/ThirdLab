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
    result = abs(derivative(diff_func, x0, n=nodes, order=nodes+2))
    result /= factorial(nodes+1)
    for x in x_values:
        result *= (x_star - x)
    return result


a = 0.4
b = 0.9
x1 = 0.52
x2 = 0.42
x3 = 0.87
x4 = 0.67
x_arr_stars = [0.42, 0.87, 0.67]

x_arr = linspace(a, b, num=15)
difTable = finite_differences_table(x_arr_stars)
h = x_arr[1] - x_arr[0]


