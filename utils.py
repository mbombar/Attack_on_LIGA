from linear_rank_metric import from_matrix_representation
from skew_polynomials import reduction
from sage.matrix.constructor import matrix
from sage.matrix.special import block_matrix
from sage.modules.free_module_element import vector
from copy import copy

from math import log

def decode_supercode(Ts, c, g, k, d, Fqm, q=None, R=None):
    """
    Tk_pub = trace of the public key over Fqm. It is a vector of length n
    and coordinates in Fqm.

    c is the ciphertext

    g is the support of the gabidulin, k it's dimension

    d is the public rank of the small error in the ciphertext

    R is a ring of linear polynomial, because Sage crashes due to https://trac.sagemath.org/ticket/28617.
    """
    n = len(g)
    if not q:
        q = Fqm.characteristic()

    if R:
        s = R.twisting_morphism()
    else:
        s = Fqm.frobenius_endomorphism()
        R = Fqm["y", s]  # skew polynomials

    Tinterp = []
    for Tk_pub in Ts:
        pointsT = list(zip(g, Tk_pub))
        T = reduction(R.lagrange_polynomial(pointsT))
        Tinterp.append(T)
    G = matrix(Fqm, n, k + d, lambda i, j: (s ** (j))(g[i]))
    S = [G]
    for T in Tinterp:
        Tm = matrix(Fqm, n, d + 1, lambda i, j: (T(g[i])) ** (q ** j))
        S.append(Tm)
    Cm = matrix(Fqm, n, d + 1, lambda i, j: (s ** j)(c[i]))
    S.append(Cm)
    S = block_matrix(Fqm, 1, 2 + len(Ts), S)
    S = S.right_kernel()
    S = S.basis_matrix().row(0)
    N = R(S.list_from_positions(range(len(S) - (d + 1))))
    V = R(S.list_from_positions(range(len(S) - (d + 1), len(S))))
    f, rem = (-N).left_quo_rem(V)
    if rem == 0:
        return (f, N, V, S)
    else:
        print(S)
        raise

def trace(Fl,Fs,vec):
    # Trace operation from a vector vec over Fl to a projection over Fs
    val_mu = log(Fl.cardinality(),Fl.characteristic())
    val_m = log(Fs.cardinality(),Fs.characteristic())
    u = int(val_mu / val_m)
    qm = Fs.cardinality()
    output = matrix(Fl,1,vec.ncols())
    for ii in range(vec.ncols()):
        output[0,ii] = sum([vec[0,ii]**(qm**j) for j in range(u)])
    return output


def rand_vec(Fl,Fs,rk,leng):
    # Generate a random vector of length leng over Fl with rank rk over Fs
    u = Fl.degree() // Fs.degree()
    # 1. Compute echelonizable matrix with rank rk
    rank_false = True
    while rank_false:
        A_mat = matrix.random(Fs,u,rk)
        if A_mat.rank()==rk:
            rank_false = False

    rank_false = True
    while rank_false:
        B_mat = matrix.random(Fs,rk,leng)
        if B_mat.rank()==rk:
            rank_false = False

    out_mat = A_mat*B_mat
    out_vec = from_matrix_representation(out_mat, Fl)
    return out_vec

def dual_basis(Fqm, Fq, basis):
    # Generate a dual basis
    basis = vector(basis[0,:])
    m = Fqm.degree() // Fq.degree()
    entries = matrix(Fqm, [xi*xj for xi in basis for xj in basis])
    entries = (Fq**(m**2))(vector(trace(Fqm, Fq, entries)))
    B = matrix(Fq, len(basis), entries).inverse()
    db = [sum(x * y for x, y in zip(col, basis)) for col in B.columns()]
    return db
