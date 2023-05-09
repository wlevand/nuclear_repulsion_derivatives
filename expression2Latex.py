# Nuclear repulsion derivatives expressions to latex string
# For usage examples, see README

# construction of a single expression string
def exr2Latex(expression):

    # if expression is multiplied by zero
    if expression['coefficient'] == 0:
        return '0'

    multiplier = '' if expression['coefficient'] == 1 else f'{expression["coefficient"]}'

    sign = '' if expression['sign'] == 1 else '-'
    charges = f'\\frac{{ Z_{expression["atoms_pair"][0]} Z_{expression["atoms_pair"][1]} }} '

    denominator = f'{{ R_{{ {expression["atoms_pair"][0]}{expression["atoms_pair"][1]} }}^{{ {expression["power"]} }} }}'

    vardict = {0: 'x', 1: 'y', 2: 'z'}
    polynomial = ''

    if expression["factoredPoly"]:

        for cart_comp in expression['factoredPoly']:

            polynomial += f'({vardict[cart_comp]}_{expression["atoms_pair"][0]} - {vardict[cart_comp]}_{expression["atoms_pair"][1]})'

            if expression['factoredPoly'][cart_comp] > 1:
                polynomial += f'^{expression["factoredPoly"][cart_comp]}'

    return sign + multiplier + polynomial + charges + denominator


# returns a LaTeX string for a list of variables
# variable: tuple/list of tuples - ((0, 0), (1, 5))
def vars2Latex(variable):

    vardict = {0: 'x', 1: 'y', 2: 'z'}
    res = ''

    for i in variable:
        res += f'{vardict[i[0]]}_{i[1]} '

    return '$' + res + '$'


# construction of a string for a sum of expressions
def exprSum2Latex(exprSum):
    res = ''
    count = len(exprSum)

    for i in exprSum:
        addition = exr2Latex(i)
        res = res + addition if addition[0] != '-' else res[:-1] + addition
        count -= 1

        if count > 0:
            res += '+'

    return '$$' + res + '$$'