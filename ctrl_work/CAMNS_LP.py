import numpy as np
from numpy.core.defchararray import index
import scipy
import scipy.stats as sps
import cvxpy as cp

def is_ext_pt(C, d, alpha, tol):
    L, D = np.shape(C)
    T = C[(C@alpha + d) < tol, :]
    sin_value = np.abs(np.linalg.eig(T))
    temp = sin_value / np.sum(sin_value)
    ran = np.sum(temp > tol)
    return ran == D


def CAMNS_LP(X, N):
    TOL_LP = 1e-3
    TOL_EXT = 1e-6
    TOL_ZEROS = 1e-6

    L, M = np.shape(X)
    
    Xn = X[np.abs(np.sum(X, axis=1)) > TOL_ZEROS]
    LL = np.sum(np.abs(np.sum(X, axis=1)) > TOL_ZEROS)

    d = Xn @ np.ones(M) / M

    C, Sigma, V = scipy.sparse.linalg.svds(Xn - (d @ np.ones(M)),N - 1)

    el = 0
    Q1 = np.zeros(LL)
    hS = []
    lp_cnt = 0
    while el < N:
        w = sps.norm.rvs(size=LL)
        w2 = np.dot(Q1, w)
        r = w - w2 * Q1

        b = -C.T @ r
        A = -C.T.copy()
        c = d.copy()

        x = cp.Variable(len(b))
        prob = cp.Problem(cp.Minimize(c.T@x),
                 [A @ x <= b])
        alpha1 = prob.solve()
        x1 = x.value.copy()


        x = cp.Variable(len(b))
        prob = cp.Problem(cp.Minimize(c.T@x),
                 [A @ x <= -b])
        alpha2 = prob.solve()
        x2 = x.value

        if el == 0:
            if is_ext_pt(C, d, alpha1, TOL_EXT):
                hS.append(C @ alpha1 + d)

            if is_ext_pt(C, d, alpha2, TOL_EXT):
                hS.append(C @ alpha2 + d)
        else:
            p_star = abs(r.T @ (C @ alpha1 + d))
            q_star = abs(r.T @ (C @ alpha2 + d))

            if p_star / np.linalg.norm(r) * np.linalg.norm(C@alpha1 + d) >= TOL_LP:
                if is_ext_pt(C, d, alpha1, TOL_EXT):
                    hS.append(C @ alpha1 + d)

            if q_star / np.linalg.norm(r) * np.linalg.norm(C@alpha2 + d) >= TOL_LP:
                if is_ext_pt(C, d, alpha2, TOL_EXT):
                    hS.append(C @ alpha2 + d)
        
        el = len(hS) # len(hS[0])
        if el > 0:
            Q1, R = np.linalg.qr(np.array(hS))
    
    #if el > N:
    #    hS=(p_star > q_star)*hS[:, 1:1:N] + (p_star < q_star)*hS[:, [1:N-1, N+1]]

    Y = np.zeros(L, N)
    Y[index, :] = hS.copy()
    hS = Y.copy()
    return hS
