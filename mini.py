# mini.py -- the relations of the big proof at small weight, with the whole matrix printed (needs sympy)
import itertools, sys
import sympy as sp
from collections import Counter

Z = 4                                   # letters: 0,1,2,3 <-> dt/(t - i^k)  (roots 1, i, -1, -i);  4 <-> dt/t
L, pi, G, z3 = sp.symbols('L pi G z3', positive=True)
I = sp.I
zeta = {2: pi**2/6, 3: z3}
beta = {1: pi/4, 2: G, 3: pi**3/32}

def Li1(s, k):                          # Li_s(i^k): the only "known" values (depth one)
    k %= 4
    if k == 0: return zeta[s]
    m1 = -L if s == 1 else -(1 - sp.Rational(1, 2**(s-1)))*zeta[s]          # Li_s(-1)
    return m1 if k == 2 else m1/2**s + (I if k == 1 else -I)*beta[s]

conv  = lambda w: w[0] != Z and w[-1] != 0              # first letter not 0, last letter not 1
words = lambda n: [w for w in itertools.product(range(5), repeat=n) if conv(w)]
depth = lambda w: sum(a != Z for a in w)
name  = lambda w: ','.join({0: '1', 1: 'i', 2: '-1', 3: '-i', 4: '0'}[a] for a in w)
fl    = lambda f: 'Li(' + ','.join(map(str, f[0])) + ';' + ','.join({0: '1', 1: 'i', 2: '-1', 3: '-i'}[z] for z in f[1]) + ')'

def to_li(w):                           # word -> (s, z):  I(w) = (-1)^depth * Li_s(i^z)
    nz = [(p, a) for p, a in enumerate(w) if a != Z]; k = len(nz)
    c, s = [0]*(k+1), [0]*(k+1)
    for j in range(1, k+1):
        p, a = nz[j-1]; q = nz[j][0] if j < k else len(w)
        c[k-j+1] = -a % 4; s[k-j+1] = q - p
    return tuple(s[1:]), tuple((c[m] - c[m-1]) % 4 for m in range(1, k+1))
def to_word(ss, zs):
    w = []
    for j in range(1, len(ss)+1): w.append(-sum(zs[:len(ss)-j+1]) % 4); w += [Z]*(ss[len(ss)-j]-1)
    return tuple(w)
def lis(n):
    def comps(n):
        if n == 0: yield (); return
        for f in range(1, n+1):
            for r in comps(n-f): yield (f,) + r
    return [(s, z) for s in comps(n) for z in itertools.product(range(4), repeat=len(s)) if s[0] >= 2 or z[0]]
def shuffle(u, v):
    res = Counter(); n = len(u)
    for pos in itertools.combinations(range(n+len(v)), n):
        iu, iv, ps = iter(u), iter(v), set(pos)
        res[tuple(next(iu) if i in ps else next(iv) for i in range(n+len(v)))] += 1
    return res
def qshuffle(a, b):                     # quasi-shuffle of two nested sums (merged terms add exponents, multiply arguments)
    (sa, za), (sb, zb) = a, b
    if not sa: yield b; return
    if not sb: yield a; return
    for s, z in qshuffle((sa[1:], za[1:]), b): yield (sa[:1]+s, za[:1]+z)
    for s, z in qshuffle(a, (sb[1:], zb[1:])): yield (sb[:1]+s, zb[:1]+z)
    for s, z in qshuffle((sa[1:], za[1:]), (sb[1:], zb[1:])): yield ((sa[0]+sb[0],)+s, ((za[0]+zb[0]) % 4,)+z)
PHI = {0: 4, 4: 0, 1: 3, 3: 1}          # the map t -> (1-t)/(1+t) permutes the letters 1<->0, i<->-i, and sends -1 to infinity
def dual(w):
    res = Counter()
    for ch in itertools.product(*[([(PHI[a], 1)] if a in PHI else []) + [(2, -1)] for a in w]):
        sg = (-1)**len(w)
        for _, s in ch: sg *= s
        res[tuple(c for c, _ in reversed(ch))] += sg
    return res

VAL = {}
def build(n):
    W = words(n); known = {w: -Li1(n, -w[0]) for w in W if depth(w) == 1}
    unk = [w for w in W if w not in known]; X = {w: sp.Symbol('x[' + name(w) + ']') for w in unk}
    rows = []                           # (label, {word: coefficient}, right-hand side)
    def rel(label, lhs, rhs):
        r, rhs = {}, sp.expand(rhs)
        for w, c in lhs.items():
            if c and w in known: rhs -= c*known[w]
            elif c: r[w] = r.get(w, 0) + c
        rows.append((label, r, sp.expand(rhs)))
    for p in range(1, n//2 + 1):
        for u in words(p):
            for v in words(n-p):
                if p == n-p and u > v: continue
                rel(f'sh {name(u)} | {name(v)}', shuffle(u, v), VAL[u]*VAL[v])
        for f in lis(p):
            for g in lis(n-p):
                if p == n-p and f > g: continue
                lhs = Counter()
                for s, z in qshuffle(f, g): lhs[to_word(s, z)] += (-1)**len(s)
                sg = (-1)**(len(f[0]) + len(g[0]))
                rel(f'st {fl(f)} | {fl(g)}', lhs, sg*VAL[to_word(*f)]*VAL[to_word(*g)])
    for s, z in lis(n):
        if all(x in (0, 2) for x in z):
            ts = [(s, tuple((x//2 + e) % 4 for x, e in zip(z, es))) for es in itertools.product((0, 2), repeat=len(s))]
            if all(t[0][0] >= 2 or t[1][0] for t in ts):
                lhs = Counter({to_word(s, z): (-1)**len(s)})
                for t in ts: lhs[to_word(*t)] -= 2**(n-len(s)) * (-1)**len(s)
                rel(f'di {fl((s, z))}', lhs, 0)
    for w in W:
        lhs = Counter({w: 1})
        for c, m in dual(w).items(): lhs[c] -= m
        rel(f'du {name(w)}', lhs, 0)
    cols = list(X)
    eqs = [sum(c*X[w] for w, c in r.items()) - rhs for _, r, rhs in rows]
    A, b = sp.linear_eq_to_matrix(eqs, [X[w] for w in cols]) if cols else (sp.zeros(0, 0), sp.zeros(0, 1))
    sol = sp.solve(eqs, [X[w] for w in cols], dict=True)[0] if cols else {}
    VAL.update(known)
    for w in cols: VAL[w] = sp.expand(sol.get(X[w], X[w]))
    return dict(words=W, known=known, cols=cols, rows=rows, A=A, b=b)

if __name__ == '__main__':
    for n in range(1, int(sys.argv[1]) + 1):
        d = build(n)
        rk = d['A'].rank() if d['cols'] else 0
        print(f'weight {n}: {len(d["words"])} words, {len(d["cols"])} unknown, {len(d["rows"])} relations, rank {rk}, '
              f'free {len(d["cols"]) - rk}')
        if n == 2:
            for w in d['cols']: print(f'   I({name(w)}) =', VAL[w])
