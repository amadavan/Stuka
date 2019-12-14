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

We provide a structure that contains all problem information to pass to the solvers. The members shown here should be all that are required to define a linear program. Note that we use a :code:`shared_ptr` for each the attributes. Leaving an attribute empty amounts to setting the :code:`shared_ptr` to be the :code:`nullptr`. This is the default! So all you need to do is set the values you care about in the :code:`struct` and you are ready to move on the next step.

.. doxygenstruct:: stuka::LP::LinearProgram
    :members:

Solving linear programs
-----------------------

In the previous step, we constructed a :code:`LinearProgram` and now we will explore how to solve this. Well, it is in fact very simple. Emulating routine from MATLAB and SciPy, we provide a utility function called :code:`linprog` that you can use to solve your problem. This requires the header file :code:`<stuka/util/functions.h>` which provides the following useful function,

.. doxygenfunction:: stuka::util::linprog(const LP::LinearProgram&, const Options&)
    :outline:


Available solvers and their options
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

- Mehrotra's Predictor-Corrector (stuka::Solver::MPC)
- Gurobi (stuka::Solver::Gurobi)

Details about the solver options will be added shortly.
