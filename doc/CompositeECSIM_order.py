"""Order of accuracy of the composite-grid curls of src/PicComposite.cpp at a
resolution change. Periodic unit box, N coarse cells per dimension, the
middle half refined by 2. Operators as in CompositeECSIM_check.py:
  Faraday  C P        (nodes -> leaf cell centers)
  Ampere   V^-1 (C P)^T W_B   (leaf cell centers -> unknown nodes), V = P^T w~
  curl curl = Ampere(Faraday)
applied to smooth periodic fields and compared with the analytic result.
usage: python3 CompositeECSIM_order.py 2|3 N1 N2 ... [uniform]"""
import sys, itertools
import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spl

d = int(sys.argv[1])
refine = sys.argv[-1] != "uniform"
Ns = [int(a) for a in sys.argv[2:] if a != "uniform"]
tp = 2 * np.pi
# E_a = sin(2 pi k_a.x + ph_a)
K = np.array([[1, 2, 1], [2, -1, 1], [1, 1, -2]], float)
PH = np.array([0.3, 1.1, 2.0])
if d == 2:
    K[:, 2] = 0

def field(x):                     # x: (n, 3) -> (n, 3)
    return np.sin(tp * x @ K.T + PH)
def curl(x):
    c = tp * np.cos(tp * x @ K.T + PH)          # c[:, a] * K[a, j] = d_j E_a
    dE = lambda a, j: c[:, a] * K[a, j]
    return np.stack([dE(2, 1) - dE(1, 2), dE(0, 2) - dE(2, 0), dE(1, 0) - dE(0, 1)], 1)
def curlcurl(x):                  # grad div - laplace
    E = field(x)
    out = np.zeros_like(E)
    for a in range(3):
        out[:, a] = tp ** 2 * (K[a] @ K[a]) * E[:, a]
        for b in range(3):
            out[:, a] -= tp ** 2 * K[b, a] * K[b, b] * E[:, b]
    return out

def build(Nc):
    L = 2 * Nc; lo, hi = L // 4, 3 * L // 4
    lo -= lo % 2; hi -= hi % 2
    offs = list(itertools.product(range(2), repeat=d))
    cells = []
    for c in itertools.product(range(0, L, 2), repeat=d):
        if refine and all(lo <= ci < hi for ci in c):
            cells += [(tuple(ci + oi for ci, oi in zip(c, o)), 1) for o in offs]
        else:
            cells.append((c, 2))
    corners = lambda c, h: [tuple((ci + h * oi) % L for ci, oi in zip(c, o)) for o in offs]
    exp = sorted({n for c, h in cells for n in corners(c, h)})
    iexp = {n: i for i, n in enumerate(exp)}
    onb = lambda n: refine and all(lo <= ni <= hi for ni in n) and any(ni in (lo, hi) for ni in n)
    hs = {n for n in exp if onb(n) and any(ni % 2 for ni in n)}
    unk = [n for n in exp if n not in hs]
    iunk = {n: i for i, n in enumerate(unk)}
    ne, nu, nc = len(exp), len(unk), len(cells)
    r, cc, v = [], [], []
    for n in exp:
        if n in hs:
            odd = [k for k in range(d) if n[k] % 2]
            for sgn in itertools.product((-1, 1), repeat=len(odd)):
                m = list(n)
                for k, s in zip(odd, sgn):
                    m[k] = (m[k] + s) % L
                r.append(iexp[n]); cc.append(iunk[tuple(m)]); v.append(0.5 ** len(odd))
        else:
            r.append(iexp[n]); cc.append(iunk[n]); v.append(1.0)
    P1 = sp.csr_matrix((v, (r, cc)), shape=(ne, nu))
    P = sp.kron(sp.identity(3), P1).tocsr()
    dx = 1.0 / L
    Dk = []
    for k in range(d):
        r, cc, v = [], [], []
        for ic, (c, h) in enumerate(cells):
            for o, n in zip(offs, corners(c, h)):
                r.append(ic); cc.append(iexp[n]); v.append((1 if o[k] else -1) / (2 ** (d - 1) * h * dx))
        Dk.append(sp.csr_matrix((v, (r, cc)), shape=(nc, ne)))
    Dk += [sp.csr_matrix((nc, ne))] * (3 - d)
    Z = sp.csr_matrix((nc, ne))
    C = sp.bmat([[Z, -Dk[2], Dk[1]], [Dk[2], Z, -Dk[0]], [-Dk[1], Dk[0], Z]]).tocsr()
    vol = np.array([(h * dx) ** d for c, h in cells])
    wt = np.zeros(ne)
    for c, h in cells:
        for n in corners(c, h):
            wt[iexp[n]] += (h * dx) ** d / 2 ** d
    V = P1.T @ wt
    CP = (C @ P).tocsr()
    pad = lambda a: np.c_[a, np.zeros((len(a), 3 - d))]
    xn = pad(np.array(unk, float) * dx)
    xc = pad(np.array([[ci + h / 2 for ci in c] for c, h in cells]) * dx)
    hc = np.array([h for c, h in cells])
    # distance (in coarse cells) of nodes and cells from the box surface
    def dist(x):       # Chebyshev distance to the surface of the box, fine index
        x = np.asarray(x, float)
        ins = np.all((x >= lo) & (x <= hi), axis=1)
        dout = np.max(np.maximum(np.maximum(lo - x, x - hi), 0), axis=1)
        din = np.min(np.minimum(x - lo, hi - x), axis=1)
        return np.where(ins, din, dout)
    dn = dist(np.array(unk)) if refine else np.full(nu, 99.)
    dc = dist(np.array([[ci + h / 2 for ci in c] for c, h in cells])) if refine else np.full(nc, 99.)
    ua = np.array(unk)
    nsurf = np.sum((ua == lo) | (ua == hi), axis=1) * (dn < 0.5)
    return CP, vol, V, xn, xc, dn, dc, hc, nu, nc, nsurf

def norms(err, w, mask):
    e = np.abs(err).max(axis=1)
    if not mask.any():
        return np.nan, np.nan
    return e[mask].max(), np.sum((w * e)[mask]) / np.sum(w)   # max, L1 share of the whole domain

rows = []; extra = []
for Nc in Ns:
    CP, vol, V, xn, xc, dn, dc, hc, nu, nc, nsurf = build(Nc)
    comp = lambda f, n: f.T.reshape(-1)                  # component-major
    un = lambda a, n: a.reshape(3, n).T
    En = field(xn)
    # Faraday: error at leaf cell centers
    eF = un(CP @ comp(En, nu), nc) - curl(xc)
    # Ampere: B = analytic field at cell centers
    eA = un(CP.T @ (np.tile(vol, 3) * comp(field(xc), nc)), nu) / V[:, None] - curl(xn)
    # curl curl of nodal E
    eCC = un(CP.T @ (np.tile(vol, 3) * (CP @ comp(En, nu))), nu) / V[:, None] - curlcurl(xn)
    nearc = dc < 2; nearn = dn < 0.5                     # cells touching the interface; interface nodes
    nearn2 = dn < 2.5                                     # interface nodes and their neighbours
    # solve (V + delta^2 (CP)^T W CP) E = V (E + delta^2 curl curl E)_exact
    delta = 0.1
    A = (sp.diags(np.tile(V, 3)) + delta ** 2 * (CP.T @ sp.diags(np.tile(vol, 3)) @ CP)).tocsc()
    if d == 2 or Nc <= 16:
        Es = spl.spsolve(A, np.tile(V, 3) * comp(En + delta ** 2 * curlcurl(xn), nu))
    else:
        Es, info = spl.cg(A, np.tile(V, 3) * comp(En + delta ** 2 * curlcurl(xn), nu), rtol=1e-12, maxiter=20000)
    eS = un(Es, nu) - En
    face = nearn & (nsurf == 1); edge = nearn & (nsurf >= 2)
    extra.append((Nc, norms(eA, V, face)[0], norms(eA, V, edge)[0], norms(eCC, V, face)[0], norms(eCC, V, edge)[0],
                  norms(eS, V, ~nearn2)[0], norms(eS, V, nearn2)[0], norms(eS, V, V > 0)[1]))
    rows.append((Nc, nu,
                 norms(eF, vol, ~nearc), norms(eF, vol, nearc), norms(eF, vol, vol > 0),
                 norms(eA, V, ~nearn2), norms(eA, V, nearn), norms(eA, V, V > 0),
                 norms(eCC, V, ~nearn2), norms(eCC, V, nearn2), norms(eCC, V, V > 0)))

def table(title, i, la, lb):
    print(f"\n{title}  (absolute errors)")
    print(f"{'N':>4s} {'unknowns':>9s} | {'max '+la:>14s} ord | {'max '+lb:>14s} ord | {'max all':>9s} ord | {'L1 all':>9s} ord")
    prev = None
    for r in rows:
        a, b, c = r[i], r[i + 1], r[i + 2]
        vals = [a[0], b[0], c[0], c[1]]
        if prev:
            o = [np.log2(p / q) if q > 0 and p == p and q == q else np.nan for p, q in zip(prev, vals)]
        else:
            o = [np.nan] * 4
        print(f"{r[0]:4d} {r[1]:9d} | {vals[0]:14.3e} {o[0]:4.2f} | {vals[1]:14.3e} {o[1]:4.2f} | {vals[2]:9.3e} {o[2]:4.2f} | {vals[3]:9.3e} {o[3]:4.2f}")
        prev = vals

print(f"d = {d}, refine = {refine}; coarse N cells per dimension, middle half refined by 2")
table("Faraday curl C P E at leaf cell centers", 2, "away", "at interface")
table("Ampere curl V^-1 (CP)^T W B at nodes", 5, "away", "interface nodes")
table("curl curl E at nodes", 8, "away", "within 1 coarse cell")

print("\nInterface nodes split: on a face of the box / on an edge or corner; and the solve")
print("(V + delta^2 S) E = V (E + delta^2 curl curl E)_exact, delta = 0.1: error of E (max|E| = 1)")
hdr = ["Amp face", "Amp edge", "cc face", "cc edge", "E max away", "E max near", "E L1"]
print(f"{'N':>4s} | " + " | ".join(f"{h:>10s}  ord" for h in hdr))
prev = None
for e in extra:
    o = [np.log2(p / q) if prev and q > 0 else np.nan for p, q in zip(prev or e[1:], e[1:])]
    print(f"{e[0]:4d} | " + " | ".join(f"{v:10.3e} {oo:4.2f}" for v, oo in zip(e[1:], o)))
    prev = e[1:]
