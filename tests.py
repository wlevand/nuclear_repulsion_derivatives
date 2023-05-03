import analytical_derivatives as nnr
import numerical_derivatives as num
import expression2Latex as latexstr
import numpy as np
import itertools
import os


# parsing XYZ files of test molecules
def get_molecule_xyz(file):
    atomsDict = {'C': 12, 'H': 1, 'O': 8, 'N': 7, 'S': 16}

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


# Test molecules from molecules folder
path = "./molecules/"
dir_list = os.listdir(path)

molecules = {}

for index, file in enumerate(dir_list):
    molecules[file[:-4]] = (get_molecule_xyz(path + file))

# print(molecules.keys())
# dict_keys(['acetone', 'AcOH', 'CH4', 'H2O2', 'H2S', 'HCOOH', 'HF', 'HONO'])

# choose a test molecule
molecule, charges = molecules['HCOOH']


# -----------------------------------------------------------


def test_allmolecules_allvars():
    errors = []

    for index, k in molecules.items():
        print(index)
        a, n = {}, {}
        molecule_cur, charges_cur = (np.array(j) for j in k)
        start = nnr.startingExp(molecule_cur)
        allvars = [(i, j) for i in range(3) for j in range(len(molecule_cur))]
        first_alvars = [(i,) for i in allvars]
        second_alvars = [(x, y) for x, y in itertools.combinations_with_replacement(allvars, 2)]
        third_alvars = [(x, y, z) for x, y, z in itertools.combinations_with_replacement(allvars, 3)]

        # list of combinations of variables
        curtest = third_alvars

        for j in curtest:
            allders = nnr.general_derivative_Expression(start, list(j))
            analytical = nnr.general_derivative_Evaluation(allders, molecule_cur, charges_cur)
            numerical = num.numericalDerGeneral(molecule_cur, charges_cur, 10 ** (-7), list(j))

            a[j] = round(analytical, 4)
            n[j] = round(numerical, 4)

        for key in a:
            if a[key] != n[key]:
                errors.append(f"Expected {n[key]} but got {a[key]} for key {key}, molecule {index}")

        if errors:
            raise AssertionError(f"{len(errors)} assertions failed:\n" + "\n".join(errors))
        assert a == n


def test_SingleMolecule_AnyOrderAllVars():
    # collecting results in the list
    info = []

    allvars = [(i, j) for i in range(3) for j in range(len(molecule))]
    first_alvars = [(i,) for i in allvars]
    second_alvars = [(x, y) for x, y in itertools.combinations_with_replacement(allvars, 2)]
    third_alvars = [(x, y, z) for x, y, z in itertools.combinations_with_replacement(allvars, 3)]
    fourth_alvars = [(x, y, z, w) for x, y, z, w in itertools.combinations_with_replacement(allvars, 4) if x <= y]

    # list of combinations of variables
    curtest = fourth_alvars

    zerothOrderExpression = nnr.startingExp(molecule)
    a, n = {}, {}
    print('start \n', latexstr.exprSum2Latex(zerothOrderExpression))

    for i in curtest:
        # list of expressions (terms of the final sum)
        alldersExpr = nnr.general_derivative_Expression(zerothOrderExpression, list(i))

        # evaluation of the sum above
        analytical = nnr.general_derivative_Evaluation(alldersExpr, molecule, charges)
        numerical = num.numericalDerGeneral(molecule, charges, 10**(-5), list(i))

        a[i] = round(analytical, 5)
        n[i] = round(numerical, 5)

        if a[i] != n[i]:
            info.append(f"Analytical is {a[i]} but numerical is {n[i]} for variables {latexstr.vars2Latex(i)} \\\\ \n" +
                        'Analytical expression: ' + latexstr.exprSum2Latex(alldersExpr) + '\n')

        else:
            info.append(f"Analytical is {a[i]} and numerical is {n[i]} for variables {latexstr.vars2Latex(i)} \\\\ \n" +
                        'Analytical expression: ' + latexstr.exprSum2Latex(alldersExpr) + '\n')

    for i in info: print(i)

    assert a == n


# curtest=[((0, 0), (0, 0), (2, 2))]
# curtest=[((2, 2), (0, 0), (0, 0))]
curtest = [((0, 3), (0, 2), (1, 2))]
curtest = [((0, 3),)]


def test_one_der():
    # print('\n', molecule)
    for i in curtest:
        print(latexstr.vars2Latex(i))

    info = []

    zerothOrderExpression = nnr.startingExp(molecule)
    a, n = {}, {}
    print('start \n', latexstr.exprSum2Latex(zerothOrderExpression))

    for i in curtest:
        print(i)
        # list of expressions (terms of the final sum)
        alldersExpr = nnr.general_derivative_Expression(zerothOrderExpression, list(i))

        # evaluation of the sum above
        analytical = nnr.general_derivative_Evaluation(alldersExpr, molecule, charges)
        numerical = num.numericalDerGeneral(molecule, charges, 0.0001, list(i))

        a[i] = round(analytical, 4)
        n[i] = round(numerical, 4)

        if a[i] != n[i]:
            info.append(f"Analytical is {a[i]} but numerical is {n[i]} for variables {latexstr.vars2Latex(i)} \\\\ \n" +
                        'Analytical expression: ' + latexstr.exprSum2Latex(alldersExpr) + '\n')

        else:
            info.append(f"Analytical is {a[i]} and numerical is {n[i]} for variables {latexstr.vars2Latex(i)} \\\\ \n" +
                        'Analytical expression: ' + latexstr.exprSum2Latex(alldersExpr) + '\n')

    for i in info: print(i)

    assert a == n


## example of use

# molecule
# [[ 0.        0.419376  0.      ]
#  [-1.025551 -0.441399  0.      ]
#  [ 1.153173  0.112123  0.      ]
#  [-0.374635  1.445556  0.      ]
#  [-0.646343 -1.327598  0.      ]]

# charges
# [12, 8, 8, 1, 1]

# x_1, x_0, y_0, z_2
variables = [((1, 3), (0, 1)), ((0, 0), (2, 2))]

# starting expression (zeroth order)
zerothOrderExpression = nnr.startingExp(molecule)

# list of expressions (terms of the final sum)
alldersExpr = nnr.general_derivative_Expression(zerothOrderExpression, variables)

# evaluation of the sum above
analytical = nnr.general_derivative_Evaluation(alldersExpr, molecule, charges)
