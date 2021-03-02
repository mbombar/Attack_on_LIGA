from utils import decode_supercode, trace, rand_vec, dual_basis
from sage.modules.free_module_element import vector
from gabidulin_code import GabidulinCode
from linear_rank_metric import from_matrix_representation, to_matrix_representation, rank_weight
from skew_polynomials import vanishing_space

import time


# SkewPolynomials are still an experimental feature.
# Comment out to display FutureWarnings.
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)

repair = True
zeta = 2

#set_random_seed(1234)

p = 2
s = 1
q = p^s

# m, k, w, u, version =  92, 53, 27, 5, 'LIGA-128'
# m, k, w, u, version = 120, 69, 35, 5, 'LIGA-192'
# m, k, w, u, version = 148, 85, 43, 5, 'LIGA-256'

m, k, w, u, version = 10, 4, 4, 3, 'LIGA, toy example'  # toy example

n = m

if repair:
    print(f"# ----------- Repaired Version ({version}, zeta={zeta}) -----------------")
else:
    print("# ----------  Original Version -----------------")


print(f"\nm={m}\nk={k}\nw={w}\nu={u}")

assert w < n-k
assert w > (n-k)//2
assert (n-k-w)//2 > 0


# Define all fields:
Fq = GF(q)
Fqm = GF(q^m)
Fqmu = GF(q^(m*u))

twist = Fqm.frobenius_endomorphism(s)


print("\n# ----------- Key Generation -------------------")
# 1. Choose g at random with rank(g)= n
g = rand_vec(Fqm,Fq,rk=n,leng=n)


# 2. Choose x at random such that the last u positions form a basis of Fqmu over Fqm
x_vec = block_matrix(Fqmu,[matrix.random(Fqmu,1,k-u),matrix(rand_vec(Fqmu,Fqm,rk=u,leng=u))],nrows=1,subdivide=False)


# LIGA implementation of the error
G_A = matrix.random(Fqm, zeta, w)
while G_A.rank()<zeta:
    G_A = matrix.random(Fqm, zeta, w)
good = [False*u]
while not all(good):
    Sprime = matrix.random(Fqm, u, zeta)*G_A
    while Sprime.rank()<zeta:
        Sprime = matrix.random(Fqm, u, zeta)*G_A
    good = list(map(lambda s: rank_weight(vector(s)), Sprime))
s = from_matrix_representation(Sprime)

# 4. Choose an invertible matrix P at random
P = matrix.random(Fq,n,n)
while not (P.rank() == n):
    P = matrix.random(Fq,n,n)

# 5. Build Gabidulin code
G = GabidulinCode(Fqm, n, k, Fq, twist, evaluation_points=g)
G_gab = G.generator_matrix()

# 6. Generate vector z
z_vec = block_matrix(Fqmu, [matrix(s), matrix([0]*(n-w))], nrows=1)*P.inverse()


# 7. Generate k_pub
k_pub = x_vec*G_gab + z_vec

# 8. Compute t_pub
t_pub = (n-w-k)//2

print("Key Generation done")

print("\n# ------------------ Encryption --------------------------")
# 0. Construct a random message:
m_vec = vector(Fqm, [Fqm.random_element() for ii in range(k-u)] + [0]*u)

# 1. Choose alpha at random
alpha = Fqmu.random_element()

# 2. Choose e such that rank_q(e) = t_pub
e_vec = rand_vec(Fqm,Fq,t_pub,n)

# 3. Build generatohar matrix of Gab
G = GabidulinCode(Fqm, n, k, Fq, twist, evaluation_points=g)
G_gab = G.generator_matrix()
E = G.encoder("VectorEvaluation")

# 4. Calculate ciphertext
trace_alpha_kpub = (Fqm^n)(vector(trace(Fqmu,Fqm,alpha * k_pub)))

c_vec = E.encode(m_vec) + trace_alpha_kpub + e_vec

print("Encryption done")

print("\n# ------------------- Decryption ---------------------------------")
# 1. Compute cP
c_prime = (c_vec*P)[w:]

# 2. Build code G'
g_prime = (g*P)[w:]
Gab_prime = GabidulinCode(Fqm, n-w, k, Fq, evaluation_points=g_prime)
G_prime = Gab_prime.generator_matrix()

# 3. Decode c' ind Gab' to get m'
c_dec = Gab_prime.decode_to_code(c_prime, decoder_name="BerlekampWelch")
m_prime = Gab_prime.unencode(c_dec)

# 4. Retrieve alpha
x_dual = dual_basis(Fqmu, Fqm, x_vec[0,(k-u):])
alpha_hat = sum([m_prime[k-u+ii]*x_dual[ii] for ii in range(u)])

# 5. Calculate m
m_hat = m_prime - vector(trace(Fqmu,Fqm,alpha_hat*x_vec))

print(f"Encryted message equals the decrypted message ? --> {m_hat == m_vec}")
assert m_hat == m_vec

print("\n# ------------------- Attack ---------------------------------")
start = time.time()



## Secret: P, s, s1, z_vec, x_vec
## Public: k_pub, g, G_gab, t_pub, k, n, w, u, Fq, Fqm, Fqmu


# Ring of linearized polynomials
f = Fqm.frobenius_endomorphism(Fq.degree())
R.<y> = Fqm['y', f]

d = t_pub

# Compute a basis of Fqmu over Fqm
gamma = matrix(rand_vec(Fqmu, Fqm, rk=u, leng=u))
print("\n# ------------------- Step 1: Get rid of the small error ---------------------------------")

# 1. Compute the supercode, then apply modified Berlekamp-Welch algorithm
Ts = []
for s in range(zeta):
    # beta = Fqmu.random_element()
    beta = gamma[0][s]
    Ts.append((Fqm^n)(vector(trace(Fqmu, Fqm, beta*k_pub))))

assert len(Ts) == zeta

# V is the annulator polynomial of the additional error vector
(f, N, V, S) = decode_supercode(Ts, c_vec, g, k, d, Fqm, q, R)

# 2. Compute the support of the error
error_support = vanishing_space(V, Fqm, Fq)

# Check that we computed the support of the error in the ciphertext
E_matrix = to_matrix_representation(e_vec, Fq)
Priv_supp = E_matrix.column_space()

assert error_support == Priv_supp

# 3. Syndrome decoding

# Now we will try to perform a syndrom decoding algorithm on the supercode Cpub = Gab + \sum <Trace(\gamma_i*Kpub)>
# First, we need to compute a parity-check matrix of this code.
# A generator matrix of the sum Cpub is the superpositionf of a generator matrix Mg and Tk_pub.
# Then we can just perform gaussian elimination to find a parity-check matrix.

Mg = G.generator_matrix()  # generator matrix of the Gabidulin code
blocks = [[Mg]]
for s in range(zeta):
    blocks.append([matrix(Ts[s])])
Sg = block_matrix(blocks)  # generator matrix of Cpub
Sg = Sg.echelon_form()  # standard form of the previous generator matrix of Cpub
H = Sg.submatrix(0, k+zeta, k+zeta, n-k-zeta)
H = block_matrix([[-H.transpose(), matrix.identity(n-k-zeta)]])  # parity-check matrix of Cpub


synd = c_vec*H.transpose()  # syndrome

F = from_matrix_representation(error_support.basis_matrix().transpose())  # d elements in Fqm such that <e1,..., en> = <F1,...,Fd>

# From now on we can express every ei in terms of the Fj :
#
#  ei = eps_i,1*F1 + ... + eps_i,d*Fd where eps_i,j in Fq are to be found --> d*n unknowns in Fq
#
# Equations are the syndrom equations
#
# Let H^T = [H1, ..., H_(n-k-1)]  where Hj is the j-th column of H^T.
#
# Then Hj is a column of n elements in Fqm
#
#
# We then have (n-k-1) equations in Fqm :
#
# e1*H1,1 + ... + en*H1,n = s1
# .....
# e1*H_(n-k-1),1 + ... + en*H_(n-k-1),n = s_(n-k-1)
#
#
# Let's focus on the first equation :
# We may decompose it using the Fi :
#
# (eps_1,1*F1+...+eps_1,d*Fd)*H1,1 + ... + (eps_1,1*F1+...+eps_1,d*Fd)*H1,n = s1
#
# We may expand all the Fl*H1,j onto Fq. Each equation on Fqm gives m equations over Fq
#
# And we finally have (n-k-1)*m equations over Fq with n*d unknowns, which we can solve with linear algebra.
#
#

# Because I wrote t on the board ...
t = d

# We compute all the products
S = [[F[j]*H[i][l] for l in range(n) for j in range(t)] for i in range(n-k-zeta)]
S = [vector(v) for v in S]

# We expand everything over Fq
S = [to_matrix_representation(v, Fq) for v in S]

# We compute the matrix of the system
M = block_matrix(n-k-zeta, 1, S)

# We expand the syndrome over Fq
synd = to_matrix_representation(synd, Fq)

# And we compute the right member of the system
synd = matrix([[synd.transpose()[i][d]] for i in range(n-k-zeta) for d in range(m)])

# We solve it
eps = M.solve_right(synd).transpose()

# We can now recover the error_vector
e_rec = vector([Fqm(0)]*n)
for i in range(n):
    e_rec[i] = sum([eps[0][i*t+l]*F[l] for l in range(t)])

# We check that we recovered the right error vector
assert e_rec == e_vec

# 4. Compute the new noisy codeword
new_c = c_vec - e_rec  # = mG + Tr(alpha*k_pub)

print("\n# ------------------- Step 2: Remove the z dependency-----------------------------")

# 5. Find the set of all alpha's such that new_c - trace(alpha*k_pub) is in G.

# Compute a basis of Fqmu over Fqm
gamma = matrix(rand_vec(Fqmu, Fqm, rk=u, leng=u))

# dual basis
gamma_dual = dual_basis(Fqmu, Fqm, gamma)

gamma = vector(gamma)
gamma_dual = vector(gamma_dual)

# Compute the traces of the gamma_i*k_pub (components of the public key in the basis gamma_dual)
gk = [(Fqm^n)(vector(trace(Fqmu, Fqm, gamma[i]*k_pub))) for i in range(u)]

Hgab = G.parity_check_matrix()
Lambda = [gk[i]*Hgab.transpose() for i in range(u)]
new_synd = new_c*Hgab.transpose()

new_M = matrix(
    Fqm,
    n-k,
    u,
    lambda i, j: Lambda[j][i]
)

Em = new_M.right_kernel()

new_alpha_coord = new_M.solve_right(new_synd)

def get_new_alpha():
    random_ksi_dual = Em.random_element()
    new_alpha = random_ksi_dual + new_alpha_coord
    new_alpha = sum(new_alpha[i]*gamma[i] for i in range(u))
    return new_alpha

def get_new_trace():
    a = get_new_alpha()
    trace_a_kpub = (Fqm^n)(vector(trace(Fqmu, Fqm, a*k_pub)))
    return trace_a_kpub


def random_m_affine():
    """
    Return a random element of the affine space m + F
    """
    return E.unencode(new_c - get_new_trace())


print("\n# ------------------- Step 3: Recover the plaintext-----------------------------")

# Now we want to recover the message. We need to find a basis of E_x = {Tr(bx) | b in \cap <mu_i>^dual}
# E_x is of dimension at most u-1
# So we find u-1 vectors in E_x : For that, we just need to get 2*(u-1) random_m_affine() and then take the u-1 consecutive differences.

basis_E_x = []
for i in range(u-1):
    m1, m2 = random_m_affine(), random_m_affine()
    basis_E_x.append(matrix(m1 - m2).transpose())

basis_E_x = block_matrix(Fqm, 1, u-1, basis_E_x)

# Define Identity and zero matrices
Ik = matrix.identity(Fqm, k)
Iu = matrix.identity(Fqm, u)

Zu = matrix.zero(Fqm, nrows=u, ncols=u-1)
Zu1 = matrix.zero(Fqm, nrows=u, ncols=1)
Zk = matrix.zero(Fqm, nrows=u, ncols=k-u)

top = block_matrix([[Ik, basis_E_x]])
bot = block_matrix([[Zk, Iu, Zu]])

S = block_matrix(2, 1, [top, bot])

b = random_m_affine()
B = matrix(b).transpose()

r_member = block_matrix(2, 1, [B, Zu1])

m_rec = vector(S.solve_right(r_member).transpose())[:k]

end = time.time()

if m_rec == m_vec:
    print("I Won !!")
    print(f"Attack lasted {round(end - start, 2)} seconds !")
else:
    print("Attack failed")
