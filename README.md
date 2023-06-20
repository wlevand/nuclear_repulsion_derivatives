# Analytical geometric derivatives of nuclear-nuclear repulsion energy

Nuclear-nuclear repulsion energy contribution to the total energy:
$$E = \sum^N_{A>B} \frac{Z_AZ_B}{R_AB},$$

with

$$R^2_{AB} = (X_A-X_b)^2+(Y_A-Y_B)^2+(Z_A-Z_B)^2$$

Python script `analytical_derivatives.py` contains an implemented functionality for the calculation of **analytical** derivatives for a given molecule. 

How to use it:
``` python
import analytical_derivatives as nnr
from tests import get_molecule_xyz

molecule, charges = get_molecule_xyz('./molecules/HCOOH.xyz')

# molecule HCOOH - type: numpy.ndarray
# [[ 0.        0.419376  0.      ]
#  [-1.025551 -0.441399  0.      ]
#  [ 1.153173  0.112123  0.      ]
#  [-0.374635  1.445556  0.      ]
#  [-0.646343 -1.327598  0.      ]]

# charges - type: list or numpy.array
# [12, 8, 8, 1, 1]

# fourth order derivative w.r.t. x_1, x_0, y_0, z_2 
variables = [(0, 1), (0, 0), (1, 0)]

# starting expression (zeroth order)
zerothOrderExpression = nnr.startingExp(molecule)

# list of expressions (terms of the final sum)
alldersExpr = nnr.general_derivative_Expression(zerothOrderExpression, variables)

# evaluation of the sum above
analytical = nnr.general_derivative_Evaluation(alldersExpr, molecule, charges)
```

For **numerical** evaluation of these derivatives, use `numerical_derivatives.py`:

```python
# evaluate this derivative numerically
import numerical_derivatives as num
numerical = num.numericalDerGeneral(molecule, charges, 10 ** (-5), variables)
```

To print $\LaTeX$ equations with `expression2Latex.py`:

``` python
import expression2Latex as latexstr

# for the final expression above
latex_string = latexstr.exprSum2Latex(alldersExpr)
print(latex_string)
# $$-3(y_0 - y_1)\frac{ Z_0 Z_1 } { R_{ 01 }^{ 5 } }+15(x_0 - x_1)^2(y_0 - y_1)\frac{ Z_0 Z_1 } { R_{ 01 }^{ 7 } }$$
```
The string will have `$$` symbols in the beginning and in the end of it:

$$-3(y_0 - y_1)\frac{ Z_0 Z_1 } { R_{ 01 }^{ 5 } }+15(x_0 - x_1)^2(y_0 - y_1)\frac{ Z_0 Z_1 } { R_{ 01 }^{ 7 } }$$


## References
Yamaguchi, Y., Osamura, Y., Goddard, J. D., & Schaefer III, H. F. A New Dimension to Quantum Chemistry Analytic Derivative Methods in _Ab Initio_ Molecular Electronic Structure Theory (pp. 51–52). Oxford University Press, Inc. (1994)
