Quadratic Programming
=================================

Quadratic programs represent a relatively straightforward extension to linear programs, wherein we now allow the objective to be a quadratic function. More precisely, we define this as problem of the following form,

.. math::

    \begin{aligned}
    & \underset{x}{\text{minimize}} &\quad& x^T Q x + c^T x, \\
    & \text{subject to} && A_{eq} x = b_{eq}, \\
    & && A_{ub} x \le b_{ub}, \\
    & && \underline{x} \le x \le \overline{x}.
    \end{aligned}

Notice that the only additional parameter is :math:`Q` which is reflected in the structure you need to define in order to construct such problems.

Constructing a quadratic program
-----------------------------

As in the case of linear programs, we provide a useful structure in which you can define all relevant properties of a quadratic program to pass on to the solvers.

.. doxygenstruct:: stuka::QP::QuadraticProgram
    :members:

Solving a quadratic program
---------------------------

Similar to the linear programming case, we provide a convenience method that emulates both MATLAB and SciPy in order to solve such problems called :code:`quadprog` which also is included in the header file :code:`<stuka/util/functions.h>`.

.. doxygenfunction:: stuka::util::quadprog
    :outline:


Available solvers and their options
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

- Gurobi (stuka::Solver::Gurobi)

Unfortunately, we only have an implementation of quadratic programs for :code:`Gurobi`. We hope to remedy this in the near future. If you are interested in adding additional solvers, please read the section titled `Extending Stuka` (to be written shortly).
