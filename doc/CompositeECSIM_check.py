"""Toy composite-grid ECSIM: check the identities claimed in doc/CompositeECSIM.tex.
Periodic coarse grid (h=2 on the fine index lattice), one refined box (h=1),
ratio 2, CIC particles, FLEKS corner-averaged curl/div/grad stencils.
usage: python3 CompositeECSIM_check.py 2|3 [uniform]"""
import sys, itertools
import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spl

d = int(sys.argv[1]) if len(sys.argv) > 1 else 2
rng = np.random.default_rng(1)
Nc = 8 if d == 2 else 6                     # coarse cells per dim
L = 2 * Nc                                  # fine-lattice period
lo, hi = (4, 12) if d == 2 else (4, 8)      # refined box in fine index (even)
refine = len(sys.argv) < 3 or sys.argv[2] != "uniform"

def inbox_cell(c):          # fine cell index c (lower corner) inside box
    return all(lo <= ci < hi for ci in c)

# leaf cells: (corner, h)
cells = []
for c in itertools.product(range(0, L, 2), repeat=d):
    if refine and all(lo <= ci < hi for ci in c):
        for o in itertools.product(range(2), repeat=d):
            cells.append((tuple(ci + oi for ci, oi in zip(c, o)), 1))
    else:
        cells.append((c, 2))
offs = list(itertools.product(range(2), repeat=d))
def corners(c, h):
    return [tuple((ci + h * oi) % L for ci, oi in zip(c, o)) for o in offs]

# expanded nodes, hanging nodes, unknowns
exp = sorted({n for c, h in cells for n in corners(c, h)})
iexp = {n: i for i, n in enumerate(exp)}
def on_box_bnd(n):
    return refine and all(lo <= ni <= hi for ni in n) and any(ni in (lo, hi) for ni in n)
hang = [n for n in exp if on_box_bnd(n) and any(ni % 2 for ni in n)]
hs = set(hang)
unk = [n for n in exp if n not in hs]
iunk = {n: i for i, n in enumerate(unk)}
ne, nu, nc = len(exp), len(unk), len(cells)

# P: expanded <- unknowns (multilinear interpolation along the interface)
rows, cols, vals = [], [], []
for n in exp:
    if n in hs:
        odd = [k for k in range(d) if n[k] % 2]
        for sgn in itertools.product((-1, 1), repeat=len(odd)):
            m = list(n)
            for k, s in zip(odd, sgn):
                m[k] = (m[k] + s) % L
            rows.append(iexp[n]); cols.append(iunk[tuple(m)]); vals.append(0.5 ** len(odd))
    else:
        rows.append(iexp[n]); cols.append(iunk[n]); vals.append(1.0)
P1 = sp.csr_matrix((vals, (rows, cols)), shape=(ne, nu))
P = sp.kron(sp.identity(3), P1).tocsr()          # 3 components, component-major

# per-cell derivative rows: dk[c] = corner weights of d/dx_k (FLEKS averaging)
def deriv(k, h):
    w = np.zeros(len(offs))
    for j, o in enumerate(offs):
        w[j] = (1 if o[k] else -1) / (2 ** (d - 1) * h)
    return w
def scalar_op(k):       # cell <- expanded nodes, d/dx_k
    r, cc, v = [], [], []
    for ic, (c, h) in enumerate(cells):
        w = deriv(k, h)
        for j, n in enumerate(corners(c, h)):
            r.append(ic); cc.append(iexp[n]); v.append(w[j])
    return sp.csr_matrix((v, (r, cc)), shape=(nc, ne))
Dk = [scalar_op(k) for k in range(d)] + [sp.csr_matrix((nc, ne))] * (3 - d)
Z = sp.csr_matrix((nc, ne))
# curl: cell (Bx,By,Bz) <- node (Ex,Ey,Ez)
C = sp.bmat([[Z, -Dk[2], Dk[1]], [Dk[2], Z, -Dk[0]], [-Dk[1], Dk[0], Z]]).tocsr()
Dv = sp.hstack(Dk[:3]).tocsr()                       # cell div of nodal E
Gr1 = sp.vstack(Dk[:3]).tocsr()                      # cell grad of nodal scalar
vol = np.array([h ** d for c, h in cells], float)
WB = sp.diags(np.tile(vol, 3)); WB1 = sp.diags(vol)
wt = np.zeros(ne)
for c, h in cells:
    for n in corners(c, h):
        wt[iexp[n]] += h ** d / 2 ** d
Wt = sp.diags(np.tile(wt, 3))
CP = (C @ P).tocsr()
G = (P.T @ Wt @ P).tocsr()                           # consistent field mass matrix
GL = sp.diags(np.asarray(G.sum(axis=1)).ravel())     # lumped
print(f"d={d} refine={refine}: unknown nodes {nu}, hanging {len(hang)}, leaf cells {nc}")

def nrm(A):
    A = sp.csr_matrix(A); return np.abs(A.data).max() if A.nnz else 0.0

# 1. weak div B: (Gr1 P1)^T W_B C P == 0 ?
Gam = sp.vstack([Dk[k] @ P1 for k in range(3)]).tocsr()   # 3nc x nu
print("weak divB  max|Gam^T W_B C P|            =", nrm(Gam.T @ WB @ CP))
# 2. Gauss law: div(P G^-1 (CP)^T W_B) == 0 ?
DP = (Dv @ P).tocsr()
for name, M in (("consistent G", G), ("lumped G", GL)):
    X = spl.spsolve(M.tocsc(), (CP.T @ WB).tocsc())
    print(f"Gauss      max|D P {name}^-1 (CP)^T W_B| =", nrm(DP @ X))

# 3. one ECSIM step with particles, theta = 1/2: total energy change
theta, dt, fourpi = 0.5, 0.7, 4 * np.pi
npart = 40 * nc
def find_leaf(x):
    for ic, (c, h) in enumerate(cells):
        if all(0 <= (x[k] - c[k]) % L < h for k in range(d)):
            return ic
S = sp.lil_matrix((npart, ne)); spec = []
xs = rng.uniform(0, L, (npart, d))
for p in range(npart):
    ic = find_leaf(xs[p]); c, h = cells[ic]
    f = [((xs[p][k] - c[k]) % L) / h for k in range(d)]
    for j, n in enumerate(corners(c, h)):
        S[p, iexp[n]] = np.prod([f[k] if offs[j][k] else 1 - f[k] for k in range(d)])
S = S.tocsr(); SP = (S @ P1).tocsr()                 # gather weights on unknowns
q = np.where(rng.random(npart) < 0.5, 1.0, -1.0) * 0.01
m = np.where(q > 0, 1.0, 0.04) * 0.01
v0 = rng.normal(0, 0.1, (npart, 3)); Bp = rng.normal(0, 1.0, (npart, 3))
def cross_mat(b):
    return np.array([[0, -b[2], b[1]], [b[2], 0, -b[0]], [-b[1], b[0], 0]])
a = q * dt / (2 * m)
Rinv = np.array([np.linalg.inv(np.eye(3) + a[p] * cross_mat(Bp[p])) for p in range(npart)])
# v_bar - a v_bar x B = v0 + a E  ->  (I + a [B x]) v_bar ... (v x B = -[B x] v)
vhat = np.einsum('pij,pj->pi', Rinv, v0)
beta = a[:, None, None] * Rinv
# assemble jhat and Mtilde on unknowns (component-major)
SPd = SP.toarray()
jh = np.zeros(3 * nu); Mt = np.zeros((3 * nu, 3 * nu))
for p in range(npart):
    s = SPd[p]; nzi = np.nonzero(s)[0]
    for a_ in range(3):
        jh[a_ * nu + nzi] += q[p] * s[nzi] * vhat[p, a_]
        for b_ in range(3):
            Mt[np.ix_(a_ * nu + nzi, b_ * nu + nzi)] += q[p] * beta[p, a_, b_] * np.outer(s[nzi], s[nzi])
En = rng.normal(0, 1, 3 * nu); Bn = rng.normal(0, 1, 3 * nc)
for name, M in (("consistent G", G), ("lumped G", GL)):
    Md = M.toarray(); CPd = CP.toarray(); WBd = WB.toarray()
    A = Md + theta ** 2 * dt ** 2 * CPd.T @ WBd @ CPd + fourpi * theta * dt * Mt
    rhs = Md @ En + theta * dt * (CPd.T @ WBd @ Bn - fourpi * jh)
    Eth = np.linalg.solve(A, rhs)
    E1 = (Eth - (1 - theta) * En) / theta
    B1 = Bn - dt * CPd @ Eth
    Ep = (SPd @ Eth.reshape(3, nu).T)                 # gathered E at particles
    vbar = vhat + np.einsum('pij,pj->pi', beta, Ep)
    v1 = 2 * vbar - v0
    K = lambda v: 0.5 * np.sum(m[:, None] * v ** 2)
    WE = lambda E: E @ Md @ E / (2 * fourpi)
    # the energy measured with the consistent G (what calc_E_field_energy computes)
    WEc = lambda E: E @ G.toarray() @ E / (2 * fourpi)
    WBf = lambda B: B @ WBd @ B / (2 * fourpi)
    tot0 = WE(En) + WBf(Bn) + K(v0); tot1 = WE(E1) + WBf(B1) + K(v1)
    totc0 = WEc(En) + WBf(Bn) + K(v0); totc1 = WEc(E1) + WBf(B1) + K(v1)
    print(f"energy ({name}): rel change in own norm {abs(tot1 - tot0) / tot0:.2e},"
          f" in diagnostic (consistent) norm {abs(totc1 - totc0) / totc0:.2e}")

# where is the weak-divB residual (rows = test-function nodes)?
R = sp.csr_matrix(Gam.T @ WB @ CP)
bad = sorted({unk[i] for i in range(nu) if R[i].nnz and np.abs(R[i].data).max() > 1e-12})
print("weak divB residual rows at nodes:", bad[:12], "... total", len(bad))
print("  of which on the box boundary:", sum(on_box_bnd(n) for n in bad))
Rg = sp.csr_matrix(DP @ spl.spsolve(G.tocsc(), (CP.T @ WB).tocsc()))
badc = [cells[i] for i in range(nc) if Rg[i].nnz and np.abs(Rg[i].data).max() > 1e-12]
print("Gauss residual cells (consistent G):", len(badc), "fine:", sum(h == 1 for c, h in badc),
      "coarse:", sum(h == 2 for c, h in badc))
