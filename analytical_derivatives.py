import numpy as np
import copy

# `ExpressionDer` dictionary represents a generalised term in the full expression of
#   an analytical derivative of the nuclear-nuclear repulsion energy.
#   It consists of a polynomial part and distance R power.
#   E.g., a term -(x_0 - x_1) * Z_0*Z_1 * R^3 has (x_0 - x_1) polynomial part and R^3 part.
#   These two parts will be differentiated separately and the results will be added.

# ExpressionDer dictionary with keys 'sign', 'power', 'atoms_pair', 'factoredPoly,' 'coefficient'
#     'sign' - keeps track of the sign of the expression
#     'power' - R**(-power/2) term
#     'atoms_pair' - pair of atoms that are considered in the given term; R is a distance between pair[0] and pair[1]
#     'factoredPoly' - polynomial in a factored form; a dictionary that keeps track of variables in the polynomial;
#                   key - a tuple of 2 variables (tuples), value - number of such terms
#     'coefficient' - coefficient that is accumulated from the power derivatives
#  e.g., for term -3*(x_0-x_1)**2 * Z_0*Z_1 * R**(-5/2) :
#     'sign' = -1, power = 5, 'atoms_pair' = (0, 1), 'factoredPoly' = {((0, 0), (0, 1)): 2}, 'coefficient' = 3
#
# A variable is given by a tuple of 2 numbers, first is an index of atom (numbering from 0),
#                                              second is an index of cartesian coordinate
# x_0 = (0, 0), y_5 = (1, 5), z_3 = (2, 3)

# template for the ExpressionDer dictionary
# ExpressionDer = {'sign': 1, 'power': 1, 'atoms_pair': pair, 'factoredPoly': {}, 'coefficient': 1}


# Updating keys of ExpressionDer dictionary
def updExpression(expression, new):
    for i in new:
        ind = ['sign', 'power', 'atoms_pair', 'factoredPoly', 'coefficient'].index(i)
        expression[ind] = new[i]


# ['factoredPoly'] dictionary is headed with a tuple of (cart_A, atom_A, cart_B, atom_B).
# E.g. (0, 0, 0, 1) represents (x_0 - x_1)

# This function adds one to the value for a key of the ['factoredPoly'] dictionary
def updPolynomialDict_Add(expr, variable):
    # (cart_A, atom_A, cart_B, atom_B)
    # the key is given by pair of atoms and the cartesian component of the `variable`
    keyVar = (variable[0])

    # if the `expr` has already such term, add +1 to the value (which represents the power):
    #   (x_0 - x_1)**value
    if keyVar in expr['factoredPoly']:
        expr['factoredPoly'][keyVar] += 1

    # if not in the `expr`, then create a key:value pair
    else:
        expr['factoredPoly'][keyVar] = 1


# Differentiation of the polynomial part of the term (expression)
def derPolynom(expression, curvar):
    newX = copy.deepcopy(expression)

    # no operations on zero terms
    if newX['coefficient'] == 0.:
        return newX

    # check if the expression involves curvar atom
    if curvar[1] not in newX['atoms_pair']:
        newX['coefficient'] = 0.
        return newX

    # if no polynomial in expression
    if not newX['factoredPoly']:
        newX['coefficient'] = 0.
        return newX

    # if not yet in the expression, a term will be added
    key = (curvar[0])

    if key not in newX['factoredPoly']:

        newX['coefficient'] = 0.

    else:

        newX['coefficient'] = newX['coefficient'] * newX['factoredPoly'][key]
        newX['factoredPoly'][key] -= 1

        # e.g., for (x_a - x_b), if curvar is x_b then change sign
        #   curvar[1] - atom number of curvar
        newX['sign'] = newX['sign'] * (-1) if curvar[1] == newX['atoms_pair'][1] else newX['sign']

        if newX['factoredPoly'][key] < 1:
            del newX['factoredPoly'][key]

    return newX


# Differentiation of the 1/R**n part of the term (expression)
def derPower(expression, curvar):
    newX = copy.deepcopy(expression)

    # no operations on zero terms
    if newX['coefficient'] == 0.:
        return newX

    # check if the expression involves curvar atom
    if curvar[1] not in newX['atoms_pair']:
        newX['coefficient'] = 0.
        return newX

    # sign from R**(-n)
    nsign = newX['sign'] * (-1)

    # sign from the (x0-x1) term
    if curvar[1] == newX['atoms_pair'][1]:
        nsign = nsign * (-1)

    newX.update({'power': newX['power'] + 2,
                 'sign': nsign,
                 'coefficient': newX['coefficient'] * newX['power']})

    updPolynomialDict_Add(newX, curvar)

    return newX


# Sum of derivatives of polynomial and 1/R**n parts of the expression for all terms
def general_one_order_derivative(expressionSum, var):
    total = []
    for expr in expressionSum:
        aa, bb = copy.deepcopy(expr), copy.deepcopy(expr)

        # two new expressions
        a, b = derPolynom(aa, var), derPower(bb, var)

        if a['coefficient'] != 0.:
            total.append(a)
        if b['coefficient'] != 0.:
            total.append(b)

    return total


def one_expression_Evaluation(expression, atoms, charges):
    if expression['coefficient'] == 0.:
        return 0.

    # polynomials such as (x_0 - x_1)
    polynomial = 1.

    for k in expression['factoredPoly']:

        # e.g. (x_0 - x_1) != 0
        # print(atoms[expression['atoms_pair'][0], k], atoms[expression['atoms_pair'][1], k])
        polyn_val = atoms[expression['atoms_pair'][0], k] - atoms[expression['atoms_pair'][1], k]
        if polyn_val != 0.:

            for n in range(expression['factoredPoly'][k]):
                polynomial *= polyn_val

        else:
            return 0.

    denominator = np.linalg.norm(atoms[expression['atoms_pair'][0]] - atoms[expression['atoms_pair'][1]]) ** expression[
        'power']
    charges = charges[expression['atoms_pair'][0]] * charges[expression['atoms_pair'][1]]

    return charges * expression['sign'] * expression['coefficient'] * polynomial / denominator


def general_derivative_Expression(start, vars):
    """
    Returns the final expression
    :param start: starting list of terms (sum of expressions)
    :param vars: cartesian variables - list of tuples of tuples (one variable - 1 tuple)
    :return: new list of terms after all differentiations
    """
    current_order = start

    for j in vars:
        t1 = general_one_order_derivative(current_order, j)
        # updating starting point as next order derivative expression sum
        current_order = t1

    return current_order


def general_derivative_Evaluation(listExrp, atoms, charges):
    total = 0.

    for expression in listExrp:
        t = one_expression_Evaluation(expression, atoms, charges)
        total += t

    return total


def startingExp(molecule):
    from itertools import combinations
    pairss = list(combinations(np.arange(len(molecule)), 2))

    # list of expressions to sum up
    start = []

    for pair in pairss:
        start.append({'sign': 1, 'power': 1, 'atoms_pair': pair, 'factoredPoly': {}, 'coefficient': 1})

    return start
