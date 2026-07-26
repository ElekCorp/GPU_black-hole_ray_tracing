#!/usr/bin/env python
# coding: utf-8

# In[1]:


from sympy.diffgeom.rn import R2
from sympy.diffgeom import metric_to_Christoffel_2nd, TensorProduct, metric_to_Riemann_components,twoform_to_matrix,metric_to_Christoffel_1st
TP = TensorProduct


# In[2]:


metric_to_Christoffel_2nd(TP(R2.dx, R2.dx) + TP(R2.dy, R2.dy))


# In[3]:


metric_to_Christoffel_2nd(R2.x*TP(R2.dx, R2.dx) + TP(R2.dy, R2.dy))


# In[4]:


from sympy import symbols, pi, sqrt, atan2, cos, sin
from sympy.diffgeom import Manifold, Patch, CoordSystem
m = Manifold('M', 4)
p = Patch('P', m)
tt, xt, yt, zt = symbols('t x y z', real=True)
a, Q, rs = symbols('a Q rs')
#a=0
#Q=0

coordinates=CoordSystem("kerr",p,symbols=(tt, xt, yt, zt))

t, x, y, z=coordinates.base_scalars()


# In[5]:


m=rs/2
r=x
delta=r**2-2*m*r+a**2+Q**2
rho2=r**2+(a**2)*(cos(y)**2)


# In[6]:


dxi=coordinates.base_oneforms()


# In[7]:


dx=dxi


# In[8]:


dx


# In[9]:


from sympy import Matrix
g =  Matrix([
    [0, 0, 0, 0],
    [0, 0, 0, 0],
    [0, 0, 0, 0],
    [0, 0, 0,0]
])


# In[10]:


g[0,0]=-a**2*sin(y)**2/rho2+delta/rho2
g[1,1]=-rho2/delta
g[2,2]=-rho2
g[3,3]=(delta*a**2*sin(y)**4-a**4*sin(y)**2-2*a**2*r**2*sin(y)**2-r**4*sin(y)**2)/rho2
g[0,3]=(a**3*sin(y)**2+a*r**2*sin(y)**2-delta*a*sin(y)**2)/rho2
g[3,0]=g[0,3]
g


v0, v1, v2, v3 = symbols('v[0] v[1] v[2] v[3]', real=True)





# In[16]:


def resub(expr):
    expr=expr.subs(t,tt)
    expr=expr.subs(x,xt)
    expr=expr.subs(y,yt)
    expr=expr.subs(z,zt)
    return expr


# Direct GR geodesic equation derivation (Lagrangian method)
from sympy import diff, Rational

# Convert metric components to use standard symbols xt, yt
delta_s = resub(delta)
rho2_s = resub(rho2)
g00_s = resub(g[0,0])
g11_s = resub(g[1,1])
g22_s = resub(g[2,2])
g33_s = resub(g[3,3])
g03_s = resub(g[0,3])

# Inverse metric blocks
det_block = g00_s * g33_s - g03_s**2
inv_g00 = g33_s / det_block
inv_g33 = g00_s / det_block
inv_g03 = -g03_s / det_block

inv_g11 = 1 / g11_s
inv_g22 = 1 / g22_s

# P0 and P3
P0 = g00_s * v0 + g03_s * v3
P3 = g03_s * v0 + g33_s * v3

# Directional derivatives along velocity: d_v = v1 * d/dxt + v2 * d/dyt
dp0_dxt = diff(g00_s, xt) * v0 + diff(g03_s, xt) * v3
dp0_dyt = diff(g00_s, yt) * v0 + diff(g03_s, yt) * v3

dp3_dxt = diff(g03_s, xt) * v0 + diff(g33_s, xt) * v3
dp3_dyt = diff(g03_s, yt) * v0 + diff(g33_s, yt) * v3

dv_P0 = v1 * dp0_dxt + v2 * dp0_dyt
dv_P3 = v1 * dp3_dxt + v2 * dp3_dyt

# Geodesic accelerations (dx_res)
dx_res = [0, 0, 0, 0]
dx_res[0] = - (inv_g00 * dv_P0 + inv_g03 * dv_P3)
dx_res[3] = - (inv_g03 * dv_P0 + inv_g33 * dv_P3)

S = g00_s * v0**2 + 2 * g03_s * v0 * v3 + g33_s * v3**2

dx_res[1] = - inv_g11 * (v1 * v2 * diff(g11_s, yt) + Rational(1, 2) * v1**2 * diff(g11_s, xt) - Rational(1, 2) * diff(S + g22_s * v2**2, xt))
dx_res[2] = - inv_g22 * (v1 * v2 * diff(g22_s, xt) + Rational(1, 2) * v2**2 * diff(g22_s, yt) - Rational(1, 2) * diff(S + g11_s * v1**2, yt))


from sympy import cse, ccode, fcode, pycode, MatrixSymbol
def c_code(expr):
    return ccode(cse(expr))


def py_code(expr):
    return pycode(cse(expr))



from sympy.codegen.ast import CodeBlock, Assignment, Variable, FloatType

x0, x1, x2, x3 = symbols('x[0] x[1] x[2] x[3]')

c = CodeBlock(Assignment(x0, dx_res[0]),Assignment(x1, dx_res[1]),Assignment(x2, dx_res[2]),Assignment(x3, dx_res[3]))
print(ccode(c.cse(optimizations='basic')))


print(ccode(c.cse(optimizations='basic'),user_functions={'Pow': [(lambda x,y: (y==-1), lambda x,y: '1.0/('+x+')' ),(lambda x,y: (y.is_integer and y<0 and y!=-1), lambda x,y: '1.0/('+'*'.join(['('+x+')']*(-int(y)))+')'),(lambda x, y: (y.is_integer and y>0), lambda x, y: '*'.join(['('+x+')']*int(y))),(lambda x, y: not y.is_integer, 'pow')]}))


print(ccode(c.cse(optimizations='basic'),user_functions={'Pow': [(lambda x, y: y.is_integer, 'pown'),(lambda x, y: not y.is_integer, 'pow')]}))





