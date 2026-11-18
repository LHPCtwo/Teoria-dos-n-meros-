from math import isqrt
rows=[(0.8,800011,1249999,558710),(0.6,600011,1666657,2133279),(0.4,400009,2499947,4949817),(0.2,200003,4999933,13599765)]
def mulsoma(a,b):
    r=0; o=0
    while a:
        if a&1: r+=b; o+=1
        b+=b; a>>=1; o+=1
    return r,o
def centro_mod(N,mod,limiar=8,lim=10**8):
    sq={x*x%mod for x in range(mod)}                       # residuos quadraticos: simetricos (palindromo)
    n0=isqrt(N)+1
    offs=sorted(o for o in range(mod) if ((n0+o)**2-N)%mod in sq)   # n que podem dar y^2 (em um periodo)
    ops=mod+85; frac=len(offs)/mod
    base=0; idx=0; n=n0+offs[0]; t=n*n-N
    r=isqrt(t); y=r if r*r>=t else r+1; ysq=y*y; m=0; first=True
    while ops<lim:
        y0=y
        if (not first) and m>=limiar:
            a=y+y+m; ops+=2; pr,o=mulsoma(m,a); ops+=o; ysq+=pr; y+=m; ops+=1
        while True:
            ops+=1
            if not ysq<t: break
            ysq+=y+y+1; y+=1; ops+=2
        while True:
            ops+=3
            if not (y>0 and ysq-(y+y-1)>=t): break
            ysq-=y+y-1; y-=1; ops+=2
        ops+=1
        if ysq==t: return (n-y,n+y),ops,frac
        m=y-y0; first=False
        idx+=1
        if idx==len(offs): idx=0; base+=mod
        nn=n0+base+offs[idx]; g=nn-n                      # proximo n permitido
        a=n+n+g; ops+=2; pr,o=mulsoma(g,a); ops+=o+2      # t += g*(2n+g) por dobra-e-soma
        t+=pr; n=nn
    return None,ops,frac
print("p/sqrtN | referencia (melhor sem filtro, com comparacoes) | m=720 | m=5040 | m=55440  (somas; fracao de n que sobra)")
for f,p,q,ref in rows:
    N=p*q; cols=[]
    for mod in (720,5040,55440):
        r,o,fr=centro_mod(N,mod)
        assert r==(p,q), (f,mod,r)
        cols.append(f"{o:,} ({100*fr:.1f}%)")
    print(f"{f:7.2f} | {ref:,} | "+" | ".join(cols))
