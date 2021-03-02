from copy import copy
from linear_rank_metric import from_matrix_representation, to_matrix_representation
from sage.matrix.constructor import matrix, vector

def reduction(p, m=None, verbose=False):
    R = p.parent()
    F = R.base_ring()
    if not m:
        m = F.degree()
    coeff = [F(0)] * m
    if verbose:
        print(f"Ici, m = {m}")
        print(f"Coeff de p sont {p.coefficients(sparse=False)}")
    for (i, u) in enumerate(p.coefficients(sparse=False)):
        if verbose:
            print(f"i = {i}, on ajoute u={u} à {coeff[i%m]}")
            print(f"on obtient {coeff[i%m]+u}")
        coeff[i % m] += u
    return R(coeff)


def adjoint(p, m=None, verbose=False):
    """
    Take a skew polynomial p and return its adjoint.
    m is the degree of the relative extension. If not specified,
    the extension is taken over the prime subfield of the base_field
    """
    if verbose:
        print(f"p = {p}.")
        print("On réduit :")
    R = p.parent()
    t = p.degree()
    if not m:
        m = R.base_ring().degree()
    np = reduction(copy(p), m)
    if verbose:
        print(f"np = {np}")
    # coeff = np.coefficients(sparse=True)
    coeff = np.padded_list(m + 1)
    if verbose:
        print(f"p est de q-degré {t} et ses coeff reduits sont {coeff}")
    # if t < m:
    #     coeff += [0]*(m-t) # you may want to use method padded_list instead !
    # print(f"On complète par des 0 : {coeff}")
    # res = []
    if verbose:
        for i, u in enumerate(coeff):
            print(f"i={i}, u={u} --> m-i = {m-i} et le coeff est {R.twist_map(m-i)(u)}")
    new_coeff = list(reversed([R.twist_map(m - i)(u) for (i, u) in enumerate(coeff)]))
    if verbose:
        print(f"Les nouveaux coeffs sont {new_coeff}\n")
    new_p = R(new_coeff)
    pstar = reduction(new_p, m)
    if verbose:
        print(f"p* = {new_p}")
        print("On réduit")
        print(f"p* = {pstar}")
    return pstar


def rand_vec(Fl, Fs, rk, leng):
    # Generate a random vector of length leng over Fl with rank rk over Fs
    u = Fl.degree() // Fs.degree()
    # 1. Compute echelonizable matrix with rank rk
    rank_false = True
    while rank_false:
        A_mat = matrix.random(Fs, u, rk)
        if A_mat.rank() == rk:
            rank_false = False

    rank_false = True
    while rank_false:
        B_mat = matrix.random(Fs, rk, leng)
        if B_mat.rank() == rk:
            rank_false = False

    out_mat = A_mat * B_mat
    out_vec = from_matrix_representation(out_mat, Fl)
    return out_vec


def vanishing_space(p, Fqm=None, Fq=None):
    """
    Compute the vanishing space of the Ore polynomial p.
    """
    if not Fqm:
        Fqm = p.base_ring()
    if not Fq:
        Fq = Fqm.prime_subfield()
    m = Fqm.degree() // Fq.degree()
    d = p.degree()
    extension, to_big_field, from_big_field = Fqm.vector_space(Fq, None, map=True)
    basis = extension.basis_matrix()
    ev = from_matrix_representation(basis, Fqm)
    ev = vector(p.multi_point_evaluation(ev))
    ev = to_matrix_representation(ev, Fq)
    ev = ev.right_kernel()
    return ev

