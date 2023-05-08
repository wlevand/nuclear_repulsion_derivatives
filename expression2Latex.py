# nuclear repulsion derivatives expressions to latex string
def exrt2Latex(expression):

    # if expression is multiplied by zero
    if expression.coefficient == 0:
        return '0'

    multiplier = '' if expression.coefficient == 1 else f'{expression.coefficient}'

    sign = '' if expression.sign == 1 else '-'
    result = f'\\frac{{ Z_{expression.atoms_pair[0]} Z_{expression.atoms_pair[1]} }} '

    denominator = f'{{ R_{{ {expression.atoms_pair[0]}{expression.atoms_pair[1]} }}^{{ {expression.power} }} }}'

    vardict = {0: 'x', 1: 'y', 2: 'z'}
    # polynomial part
    uder = ''

    if expression.factoredPoly:

        for i in expression.factoredPoly.keys():

            uder += f'({vardict[i[0]]}_{i[1]} - {vardict[i[2]]}_{i[3]})'

            if expression.factoredPoly[i] > 1:
                uder += f'^{expression.factoredPoly[i]}'

    return sign + multiplier + uder + result + denominator


# returns a LaTeX string for a list of variables
def vars2Latex(variable):
    """
    :param variable: tuple/list of tuples - ((0, 0), (1, 5))
    :return:
    """
    vardict = {0: 'x', 1: 'y', 2: 'z'}
    res = ''

    for i in variable:
        res += f'{vardict[i[0]]}_{i[1]} '

    return '$' + res + '$'


def exprSum2Latex(exprSum):
    res = ''
    count = len(exprSum)

    for i in exprSum:
        addition = exrt2Latex(i)
        res = res + addition if addition[0] != '-' else res[:-1] + addition
        count -= 1

        if count > 0:
            res += '+'

    return '$$' + res + '$$'