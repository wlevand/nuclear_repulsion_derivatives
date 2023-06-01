# analytical_derivatives.py
import analytical_derivatives as nnr
from tests import get_molecule_xyz

molecule, charges = get_molecule_xyz('./molecules/HCOOH.xyz')

# molecule HCOOH
# [[ 0.        0.419376  0.      ]
#  [-1.025551 -0.441399  0.      ]
#  [ 1.153173  0.112123  0.      ]
#  [-0.374635  1.445556  0.      ]
#  [-0.646343 -1.327598  0.      ]]

# charges
# [12, 8, 8, 1, 1]

# fourth order derivative w.r.t. x_1, x_0, y_0, z_2
variables = [(0, 1), (0, 0), (1, 0)]

# starting expression (zeroth order)
zerothOrderExpression = nnr.startingExp(molecule)

# list of expressions (terms of the final sum)
alldersExpr = nnr.general_derivative_Expression(zerothOrderExpression, variables)

# evaluation of the sum above
analytical = nnr.general_derivative_Evaluation(alldersExpr, molecule, charges)

# -------------------------------------------------------
# expression2Latex.py
import expression2Latex as latexstr

# for the final expression above
latex_string = latexstr.exprSum2Latex(alldersExpr)
print(latex_string)
# $$-3(y_0 - y_1)\frac{ Z_0 Z_1 } { R_{ 01 }^{ 5 } }+15(x_0 - x_1)^2(y_0 - y_1)\frac{ Z_0 Z_1 } { R_{ 01 }^{ 7 } }$$
