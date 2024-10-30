# Project 1 COE352
# name: Alana Gaughan
# eid: arg5345

# (30 pts)  Code a callable function that will calculate the SVD for a general, NXM matrix.
# The routine should return
# (1) each matrix of the SVD decomposition,
# (2) the matrix condition number using its singular/eigenvalues and
# (3) the matrix inverse.
# If the matrix has no inverse, the routine should return with an error.
# note: You will need to use a blackbox eigenvalue/vector call for this (np.linalg.eig in numpy for example).
import numpy as np
import numpy.linalg as la
import math
import scipy.linalg as sla


def find_SVD(matrix):
    """Calculates the SVD decomposition, condition number, and matrix inverse of the given matrix

    Inputs:
    matrix: mxn numpy array

    Outputs:
    U = mxm numpy array
    V = nxn numpy array
    Sigma = mxn numpy array
    inverse = nxm numpy array
    condition_num = float
    """
    m, n = matrix.shape

    left_vals, left_vecs = la.eig(np.matmul(matrix, matrix.T))
    left_dict = {}
    for i in range(len(left_vals)):
        left_dict[left_vals[i]] = left_vecs[:,i]
    l_keys = list(left_dict.keys())
    l_keys.sort(reverse=True)
    # Sorted Dictionary
    left_dict_sorted = {i: left_dict[i] for i in l_keys}

    right_vals, right_vecs = la.eig(np.matmul(matrix.T, matrix))
    right_dict = {}
    for i in range(len(right_vals)):
        right_dict[right_vals[i]] = right_vecs[:,i]
    r_keys = list(right_dict.keys())
    r_keys.sort(reverse=True)
    # Sorted Dictionary
    right_dict_sorted = {i: right_dict[i] for i in r_keys}

    # print(f"right sorted = {right_dict_sorted}\n")
    # print(f'left sorted = {left_dict_sorted}\n')

    left_vecs = np.array(list(left_dict_sorted.values()))
    left_vals = np.array(list(left_dict_sorted.keys()))

    right_vecs = np.array(list(right_dict_sorted.values()))
    right_vals = np.array(list(right_dict_sorted.keys()))

    if m > n:
        # case if there are more rows than cols
        singular_values = np.sqrt(right_vals)
        long_sing_val = np.sqrt(left_vals)

        # create U matrix from left eigenvectors
        U = np.zeros((m, m))
        for col in range(len(left_vecs)):
            U[:, col] = left_vecs[col]

        # create V from vi = 1/si * [A].T * ui
        V = np.zeros((n,n))
        for col in range(n):
            V[:, col] = (1/singular_values[col]) * np.matmul(matrix.T, left_vecs[col])
    else:
        # case if there are more cols than rows
        singular_values = np.sqrt(left_vals)
        long_sing_val = np.sqrt(right_vals)

        # create V from right eigenvectors
        V = np.zeros((n, n))
        for col in range(len(right_vecs)):
            V[:, col] = right_vecs[col]

        # create U from graham schmidt
        U = np.zeros((m,m))
        for col in range(m):
            U[:, col] = (1/singular_values[col]) * np.matmul(matrix, right_vecs[col])

    # # create U matrix from left eigenvectors
    # U = np.zeros((m, m))
    # for col in range(len(left_vecs)):
    #     U[:,col] = left_vecs[col]
    #
    # # create V from right eigenvectors
    # V = np.zeros((n, n))
    # for col in range(len(right_vecs)):
    #     V[:,col] = right_vecs[col]

    # create Sigma from singular values
    Sigma = np.zeros((m,n))
    # print(f'len sing vals = {len(long_sing_val)}')
    for val in range(len(singular_values)):
        Sigma[val, val] = singular_values[val]


    # Find the inverse
    for s in long_sing_val:
        if math.isclose(s, 0):
            inverse = 'Error: Matrix not invertible, eigenvalue of zero exists'
        else:
            sigma_inv = np.zeros((m, n))
            for val in range(len(singular_values)):
                sigma_inv[val, val] = (1 / singular_values[val])
            # print("sigma_inv = ")
            inverse = np.matmul(V, np.matmul(sigma_inv.T, U.T))

    # calculate condition number
    max_sigma = max(singular_values)
    min_sigma = min(singular_values)
    condition_number = max_sigma/min_sigma

    return U, V, Sigma, inverse, condition_number

# (10 pts)  Compare the results from your routine to that of MATLAB or python (or any c, c++, etc SVD blackbox)
# for an invertible matrix.

# Note: For invertible matrices, these routines return a “pseudo-inverse”,
# which utilize the SVD decomposition to solve Ax=b even when A is not invertible.
# When it is invertible, it is equivalent to solving the problem exactly.


A = np.array([[3, 3, 2], [2, 3, -2]])
B = np.array([[3, 3],[2, 3], [-2, 2]])
C = np.array([[3, 4], [0, 5]])
matrix_list = [A, B, C]

print("Comparing my Function vs SciPy Routines\n")
for matrix in matrix_list:
    my_U, my_V, my_S, my_inv, my_cond = find_SVD(matrix)
    scipy_U, scipy_S, scipy_V = sla.svd(matrix)
    scipy_inv = sla.pinv(matrix)
    np_cond = la.cond(matrix)
    print(f"Matrix:\n{matrix}")
    print(f"My SVD:\nU={my_U}\nS={my_S}\nV={my_V}\n")
    print(f"SciPy SVD:\nU={scipy_U}\nS={scipy_S}\nV={scipy_V}\n")
    print(f'My Inverse:\n{my_inv}\n')
    print(f'Scipy Inverse:\n{scipy_inv}\n')
    print(f'My Condition Number: {my_cond}\n')
    print(f'Numpy Condition Number: {np_cond}\n')

