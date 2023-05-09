import analytical_derivatives as nnr
import numerical_derivatives as num
import expression2Latex as latexstr
import numpy as np
import itertools
import os


# parsing XYZ files of test molecules
def get_molecule_xyz(file):
    atomsDict = {'C': 12, 'H': 1, 'O': 8, 'F': 9, 'N': 7, 'S': 16}

    if file[-4:] != '.xyz':
        return "Not an xyz file"

    coords = []
    charges = []

    with open(file) as xyzfile:
        for line in xyzfile:
            clean = line.strip().split()

            if clean:
                if clean[0] in list(atomsDict.keys()):
                    charges.append(atomsDict[clean[0]])
                    coords.append(np.array([float(i) for i in clean[1:]]))

    return np.array(coords), charges


# -----------------------------------------------------------

# a function for testing all given variables in a list 'curtest' (defined below)
def test_SingleMolecule():
    zerothOrderExpression = nnr.startingExp(molecule)
    a, n = {}, {}
    # print('start \n', latexstr.exprSum2Latex(zerothOrderExpression))

    for i in curtest:
        # list of expressions (terms of the final sum)
        alldersExpr = nnr.general_derivative_Expression(zerothOrderExpression, list(i))

        # evaluation of the sum above
        analytical = nnr.general_derivative_Evaluation(alldersExpr, molecule, charges)
        numerical = num.numericalDerGeneral(molecule, charges, 10 ** (-5), list(i))

        a[i] = round(analytical, 4)
        n[i] = round(numerical, 4)

        if a[i] != n[i]:
            print(f"Analytical is {a[i]} but numerical is {n[i]} for variables {latexstr.vars2Latex(i)} \\\\ \n" +
                  'Analytical expression: ' + latexstr.exprSum2Latex(alldersExpr) + '\n')

    assert a == n
    print(f'Success for {len(curtest[0])} order derivatives for molecule {mol_name}')


# Test molecules from molecules folder
path = "./molecules/"
dir_list = os.listdir(path)

molecules = {}

for index, file in enumerate(dir_list):
    molecules[file[:-4]] = (get_molecule_xyz(path + file))

# print(molecules.keys())
# dict_keys([''HCOOH', 'HF'])

# choose a test molecule
mol_name = 'CH4'
molecule, charges = molecules[mol_name]

# getting lists of combinations of all variables for four orders of derivatives (1-4)
all_possible_variables = [(i, j) for i in range(3) for j in range(len(molecule))]

# list of first order variables; e.g., for HF:
# [(0, 0), (0, 1), (1, 0), (1, 1), (2, 0), (2, 1)] - x_0, x_1, y_0, y_1, z_0, z_1
firstOrd = [(i,) for i in all_possible_variables]

# test first order
curtest = firstOrd
test_SingleMolecule()

secondOrd = [(x, y) for x, y in itertools.combinations_with_replacement(all_possible_variables, 2)]

# test second order
curtest = secondOrd
test_SingleMolecule()

thirdOrd = [(x, y, z) for x, y, z in itertools.combinations_with_replacement(all_possible_variables, 3)]

# test third order
curtest = thirdOrd
test_SingleMolecule()

# test fourth order
fourthOrd = [(x, y, z, w) for x, y, z, w in itertools.combinations_with_replacement(all_possible_variables, 4)]
curtest = fourthOrd
test_SingleMolecule()

# -----------------------------------------------------------
# single derivative
curtest_oneder = [(0, 1), (0, 1), (1, 1)]
zerothOrderExpression = nnr.startingExp(molecule)
alldersExpr = nnr.general_derivative_Expression(zerothOrderExpression, curtest_oneder)

# evaluation of the sum above
analytical = nnr.general_derivative_Evaluation(alldersExpr, molecule, charges)
numerical = num.numericalDerGeneral(molecule, charges, 10 ** (-5), curtest_oneder)

# print(latexstr.exprSum2Latex(alldersExpr))
# print(latexstr.exprSum2Latex(zerothOrderExpression))