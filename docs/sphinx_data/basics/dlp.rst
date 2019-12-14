Decomposed Linear Programming
=================================

While we traditionally assume linear programs to be "easy", this section focuses on the case where this is no longer true due to the large problem size. The approach that we consider takes advantage of the particular structure that is often found in such problems to provide algorithms that are more scalable. Stated mathematically, we are still looking at linear programs which are defined as,

.. math::

    \begin{aligned}
    & \underset{x}{\text{minimize}} &\quad& c^T x, \\
    & \text{subject to} && A_{eq} x = b_{eq}, \\
    & && A_{ub} x \le b_{ub}, \\
    & && \underline{x} \le x \le \overline{x}.
    \end{aligned}

The distinction lies in the particular structure of the matrices :math:`A_{eq}` and :math:`A_{ub}`. Specifically, we have :math:`A` matrices that look like,

IMAGE

This allows us to separate the original problem into smaller problems of the form,

.. math::

    \begin{aligned}
    & \underset{x}{\text{minimize}} &\quad& c_0^T x_0 + \sum_{i = 1}^N J_i^*(x_0), \\
    & \text{subject to} && A_{0,{eq}} x_0 = b_{0,{eq}}, \\
    & && A_{0,{ub}} x_0 \le b_{0,{ub}}, \\
    & && \underline{x}_0 \le x_0 \le \overline{x}_0,
    \end{aligned}

where, we have

.. math::
   \begin{aligned}
    J_i^*(x_0) := \
    & \underset{x}{\text{minimize}} &\quad& c_i^T x_i, \\
    & \text{subject to} && C_{i,{eq}} x_0 + A_{i,{eq}} x_i = b_{i,{eq}}, \\
    & && C_{i,{ub}} x_0 + A_{i,{ub}} x_i \le b_{i,{ub}}, \\
    & && \underline{x}_i \le x_i \le \overline{x}_i.
    \end{aligned}

Notice that we decoupled the original problem into smaller subproblems that we hope are more tractable to solve that the original large-scale linear program.

Constructing a decomposed linear program
----------------------------------------

This can done very similarly to a linear program, where the attributes are now arrays corresponding to the subproblems. Note that we index the subproblems from 0 to N-1 and the final attribute corresponds to the base problem in :math:`x_0`.

.. doxygenstruct:: stuka::dLP::DecomposedLinearProgram
    :members:

Solving decomposed linear programs
----------------------------------

Once again we provide a utility function to solve decomposed linear programs. As these problems are still technically linear programs, we overload the original :code:`linprog` function for these problems, which can be called from the following function that is available in the :code:`<stuka/util/functions.h>` header file.

.. doxygenfunction:: stuka::util::linprog(const dLP::DecomposedLinearProgram&, const Options&)
    :outline:


Available solvers and their options
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

- Bender's Decomposition (stuka::Solver::Benders)
- Critical Region Exploration (stuka::Solver::CRE)
