# testing against numerical derivatives order by order

import numpy as np
from itertools import combinations


def nuc_repulsionZeroOrder(molecule, charges):
    pairsAll = list(combinations(np.arange(len(molecule)), 2))
    total = 0
    for i in pairsAll:
        total += charges[i[0]] * charges[i[1]] / np.linalg.norm(molecule[i[0]] - molecule[i[1]])
    return total


def upd1varmolUP(molecule, variable, h):
    import copy
    newmol = copy.deepcopy(molecule)
    newmol[variable[1], variable[0]] += h

    return newmol


def upd1varmolDOWN(molecule, variable, h):
    import copy
    newmol = copy.deepcopy(molecule)
    newmol[variable[1], variable[0]] -= h

    return newmol


def mol_difference(molecule, variable, h, sign):
    import copy
    newmol = copy.deepcopy(molecule)
    newmol[variable[1], variable[0]] += sign * h

    return newmol


def numericalDerGeneral(molecule, charges, h, variables):
    import analytical_derivatives as nnr

    molUp = mol_difference(molecule, variables[-1], h, sign=1)
    molDown = mol_difference(molecule, variables[-1], h, sign=-1)

    if len(variables) == 1:
        return (nuc_repulsionZeroOrder(molUp, charges) - nuc_repulsionZeroOrder(molDown, charges)) / (2 * h)

    else:
        start = nnr.startingExp(molecule)
        norderder = nnr.general_derivative_Expression(start, variables[:-1])

        return (nnr.general_derivative_Evaluation(norderder, molUp, charges) - nnr.general_derivative_Evaluation(norderder, molDown, charges)) / (2 * h)

