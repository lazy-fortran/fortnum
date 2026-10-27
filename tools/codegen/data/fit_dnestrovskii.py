"""Reproduce real-axis Dnestrovskii fits and independent high-precision oracles.

Definition: Bornatici et al., Nucl. Fusion 23 (1983), Eqs. 2.3.68–2.3.71,
https://doi.org/10.1088/0029-5515/23/9/005. Requires mpmath.
Coefficients are rounded to 22 significant decimal digits before code generation.
"""
from pathlib import Path
import mpmath as m
m.mp.dps=90
root=Path(__file__).resolve().parents[3]
bounds=[0,1,2,4,8,16,32,64,128]
def fq(x,q):
    if x == 0: return 1/(q-1)
    f=-m.sqrt(m.pi)*m.exp(-x)*m.erfi(m.sqrt(x))/m.sqrt(x)
    r=m.mpf('.5')
    while r<q:
        f=(1+x*f)/r
        r+=1
    return f
def fq_positive(z,q):
    if z == 0: return 1/(q-1)
    f=m.sqrt(m.pi)*m.exp(z)*m.erfc(m.sqrt(z))/m.sqrt(z)
    order=m.mpf('.5')
    while order<q:
        f=(1-z*f)/order
        order+=1
    return f
lines=[]
maxerr=m.mpf(0)
for q in [m.mpf('2.5'),m.mpf('3.5')]:
    for a,b in zip(bounds,bounds[1:]):
        c,err=m.chebyfit(lambda s:fq((a+b)/2+m.mpf(b-a)/2*s,q),[-1,1],33,error=True)
        maxerr=max(maxerr,err)
        lines.append(' '.join(m.nstr(t,22) for t in c))
for q in [m.mpf('2.5'),m.mpf('3.5')]:
    for a,b in zip(bounds[1:],bounds[2:]):
        c,err=m.chebyfit(lambda s:fq_positive((a+b)/2+m.mpf(b-a)/2*s,q),[-1,1],33,error=True)
        maxerr=max(maxerr,err)
        lines.append(' '.join(m.nstr(t,22) for t in c))
(root/'tools/codegen/data/dnestrovskii.dat').write_text('\n'.join(lines)+'\n')
# Independent oracle values use 90-digit defining functions, not the fits.
refs=[]
points = {m.mpf(k)/2 for k in range(257)}
points.update(m.mpf(x) for x in ['1e-8', '.3', '3', '6', '9', '12', '15', '1000'])
for edge in bounds[1:]:
    points.update([m.mpf(edge)-m.mpf('1e-8'), m.mpf(edge)+m.mpf('1e-8')])
for z in sorted({-x for x in points} | points):
    row=[z]
    for q in [m.mpf('2.5'),m.mpf('3.5')]:
        if z <= 0:
            x=-z
            re=fq(x,q)
            im=-m.pi*x**(q-1)*m.exp(-x)/m.gamma(q)
            dre=(re-fq(x,q-1)) if x else -1/((q-1)*(q-2))
            dim=im*(1-(q-1)/x) if x else m.mpf(0)
        else:
            re=fq_positive(z,q)
            im=m.mpf(0)
            dre=re-fq_positive(z,q-1)
            dim=m.mpf(0)
        row += [re,im,dre,dim]
    refs.append(' '.join(m.nstr(t,22) for t in row))
(root/'test/oracle/data/dnestrovskii.dat').write_text('\n'.join(refs)+'\n')
print('maximum Chebyshev interpolation error estimate',m.nstr(maxerr,5))
