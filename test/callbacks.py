import numpy as np
import scipy as sp
import scipy.sparse

import stukapy as st

c = np.array([2, 1])
A_ub = sp.sparse.bmat([[-1, 1],
                       [-1, -1],
                       [0, -1],
                       [1, -2]], 'csc')
b_ub = np.array([1, -2, 0, 4])

lp = st.LinearProgram(c=c, A_ub=A_ub, b_ub=b_ub)
opts = st.Options()
# opts.lp_solver = st.solver.MPC

print("No Callback")
res = st.linprog(lp, opts)
print(res)

print("SaveHDF5 Callback")
opts.callback = st.callback.SaveHDF5("test.h5", 10)
res = st.linprog(lp, opts)
print(res)

print("Functional Callback")


def fcallback(state):
    print(state['nit'])


opts.callback = st.callback.Function(fcallback)
res = st.linprog(lp, opts)
print(res)
