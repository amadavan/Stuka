Linear Programming
=================================

Linear programs are perhaps the simplest form of convex optimization that we can easily solve. These programs are defined as the following

.. math::

    \begin{aligned}
    & \underset{x}{\text{minimize}} &\quad& c^T x, \\
    & \text{subject to} && A_{eq} x = b_{eq}, \\
    & && A_{ub} x \le b_{ub}, \\
    & && \underline{x} \le x \le \overline{x}.
    \end{aligned}

Constructing a linear program
-----------------------------

We provide a structure that contains all problem information to pass to the solvers.

.. doxygenstruct:: stuka::LP::LinearProgram
    :members:

Available solvers
-----------------

- Mehrotra's Predictor-Corrector (stuka::Solver::MPC)
- Gurobi (stuka::Solver::Gurobi)
