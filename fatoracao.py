"""
Fatoracao de semiprimos impares N = p*q (p < q) por soma e resto com sinal.
Notas de trabalho: ver metodologia.md. Tudo aqui e autocontido (python3, sem bibliotecas).

Modelo de custo: contam-se somas/subtracoes. Multiplicar por 10 ou dobrar conta como soma.
A raiz isqrt(N) e usada so como LIMITE de busca, nunca para fatorar.
"""
from math import isqrt, gcd
import random, time

# ---------------------------------------------------------------- 1. compressao (roda)
# Os candidatos coprimos com M se repetem a cada M: o padrao e a "constante comprimida".
M = 30030                                   # 2*3*5*7*11*13
RES = [r for r in range(1, M) if gcd(r, M) == 1]          # 5760 posicoes por periodo

def cands_from(lo):
    """Candidatos impares > lo coprimos com 3..13 (e >= 17), gerados pela regra, sem tabela."""
    k = lo // M
    while True:
        for r in RES:
            q = k * M + r
            if q > lo and q >= 17:
                yield q
        k += 1

def todos():
    for q in (3, 5, 7, 11, 13):
        yield q
    yield from cands_from(13)

# ---------------------------------------------------------------- 2. aritmetica por soma
def mul(a, b, st):
    """a*b por dobra-e-soma. st[1] acumula somas."""
    r = 0
    while a:
        if a & 1:
            r += b; st[1] += 1
        b += b; a >>= 1; st[1] += 1
    return r

# ---------------------------------------------------------------- 3. percurso de baixo
def walk(it, N, st):
    """Para cada candidato Q (crescente) mantem S = Q*P por soma, com P comecando em N.
    O passo de P cresce (D+D) ate passar de N e desce; se S==N, P=N/Q sai do contador."""
    P = N; tw = P + P; Q = None
    for q in it:
        if Q is None:
            Q = q; S = mul(q, N, st)
        else:
            if q > P:
                return None
            for _ in range((q - Q) // 2):
                S += tw; st[1] += 1
            Q = q
        D = 2 * Q; U = 2; stk = []
        while S - D >= N:
            stk.append((D, U)); D += D; U += U; st[1] += 2
        for D, U in reversed(stk):
            if S - D >= N:
                S -= D; P -= U; st[1] += 2
        if S == N:
            return (Q, P)
        tw = P + P
    return None

# ---------------------------------------------------------------- 4. crivo a partir de p^2
def plano(n):
    s = bytearray([1]) * (n + 1)
    for i in range(3, isqrt(n) + 1, 2):
        if s[i]:
            for j in range(i * i, n + 1, 2 * i):
                s[j] = 0
    return [i for i in range(3, n + 1, 2) if s[i]]

def crivo_com(base, limite):
    """Primos impares ate 'limite' marcando so a partir de p^2 com os primos de 'base'."""
    s = bytearray([1]) * (limite + 1); marcas = 0
    for p in base:
        if p * p > limite:
            break
        for j in range(p * p, limite + 1, 2 * p):
            s[j] = 0; marcas += 1
    return [i for i in range(3, limite + 1, 2) if s[i]], marcas

def torre(N):
    """Guarda so os primos ate N^(1/8); o nivel seguinte (ate N^(1/4)) sai do crivo."""
    r2 = isqrt(N); r4 = isqrt(r2) + 1
    base0 = plano(isqrt(r4) + 1)
    return crivo_com(base0, r4)[0]

def primos_seg(limit, W, peq, cnt):
    """Crivo segmentado: so a janela de W impares existe (o buffer); candidatos = primos."""
    lo = 3
    while lo <= limit:
        hi = min(lo + 2 * W, limit + 1); marks = bytearray(W)
        for p in peq:
            if p * p >= hi:
                break
            m = ((lo + p - 1) // p) * p
            if m % 2 == 0:
                m += p
            m = max(m, p * p)
            for k in range((m - lo) // 2, W, p):
                marks[k] = 1; cnt[0] += 1
        for i in range(W):
            v = lo + 2 * i
            if v > limit:
                return
            if not marks[i]:
                yield v
        lo += 2 * W

def fatora_baixo(N, W=16384):
    """De baixo para cima: candidatos primos por crivo em janelas (memoria ~W bytes) + walk."""
    cnt = [0]; st = [0, 0]
    r = walk(primos_seg(isqrt(N), W, torre(N), cnt), N, st)
    return r, st[1] + cnt[0]

def fatora_baixo_roda(N):
    """De baixo para cima so com a roda (sem crivo): candidatos ~19% dos impares."""
    st = [0, 0]
    return walk(todos(), N, st), st[1]

# ---------------------------------------------------------------- 5. percurso do centro
def centro(N, limiar=8, lim=10**8):
    """Do centro para fora: n = (P+Q)/2, y = (Q-P)/2, N = n^2 - y^2. y^2 acumula impares.
    Se a inclinacao (passos de y por passo de n) for >= limiar, salta (y+m)^2 por dobra-e-soma
    e corrige com o resto positivo ou o oposto. Bom quando p/sqrt(N) >= ~0.5."""
    n = isqrt(N) + 1; t = n * n - N; y = 0; ysq = 0; m = 0; ops = 0; first = True
    while ops < lim:
        y0 = y
        if (not first) and m >= limiar:
            a = y + y + m; ops += 2; pr = 0; kk = m; aa = a
            while kk:
                if kk & 1:
                    pr += aa; ops += 1
                aa += aa; kk >>= 1; ops += 1
            ysq += pr; y += m; ops += 1
        while ysq < t:
            ysq += y + y + 1; y += 1; ops += 2
        while y > 0 and ysq - (y + y - 1) >= t:
            ysq -= y + y - 1; y -= 1; ops += 2
        if ysq == t:
            return (n - y, n + y), ops
        m = y - y0; first = False
        t += n + n + 1; n += 1; ops += 2
    return None, ops

# ---------------------------------------------------------------- 6. testes
def _primos(a, b):
    s = bytearray([1]) * (b + 1); s[0] = s[1] = 0
    for i in range(2, isqrt(b) + 1):
        if s[i]:
            s[i * i::i] = bytearray(len(s[i * i::i]))
    return [i for i in range(max(3, a), b + 1) if s[i]]

if __name__ == "__main__":
    random.seed(1)
    pr = _primos(3, 20000)
    print("fatora_baixo (crivo em janelas) e fatora_baixo_roda:")
    ok1 = ok2 = tot = 0; s1 = s2 = 0
    for _ in range(60):
        p, q = sorted(random.sample(pr, 2)); N = p * q
        r1, t1 = fatora_baixo(N); r2, t2 = fatora_baixo_roda(N)
        ok1 += (r1 == (p, q)); ok2 += (r2 == (p, q)); tot += 1; s1 += t1; s2 += t2
    print(f"  corretos {ok1}/{tot} e {ok2}/{tot} | somas medias {s1/tot:,.0f} e {s2/tot:,.0f}")
    print("centro (fatores proximos, p/sqrt(N) >= 0.6):")
    ok = tot = 0; sc = 0
    for _ in range(60):
        p = random.choice(pr); q = random.choice([x for x in pr if p < x <= p * 2])
        N = p * q; r, o = centro(N); ok += (r == (p, q)); tot += 1; sc += o
    print(f"  corretos {ok}/{tot} | somas medias {sc/tot:,.0f}")
    N = 999983 * 1000003
    print("N=999983*1000003:", fatora_baixo(N), "| torre guarda so os primos ate N^(1/8):", plano(isqrt(isqrt(isqrt(N))) + 1))
