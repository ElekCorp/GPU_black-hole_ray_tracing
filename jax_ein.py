from functools import partial
import jax 
import jax.numpy as jnp
from jax.scipy.linalg import solve

#jax.config.update("jax_enable_x64", True)



@partial(jax.jit, static_argnames=('model',))
def metric(coords, model=None):
    """
    Compute the metric tensor at a given coordinate.

    Parameters
    ----------
    coords : jnp.array
        Coordinates of the point at which to compute the metric.
    model :
        Function that returns the metric tensor at a given coordinate. If not specified, the Euclidean metric is used.
    
    Returns
    -------
    metric : jnp.array
        Metric tensor (two lower indices) at the given coordinate. 
    """

    if model is None:
        return jnp.eye(coords.shape[-1])
    return model(coords)

# partial derivatives of the metric
pd_metric = jax.jit(jax.jacfwd(metric), static_argnames=('model',))

@partial(jax.jit, static_argnames=('model',))
def christoffel(coords, model):
    """
    Compute the Christoffel symbols at a given coordinate.

    Parameters
    ----------
    coords : jnp.array
        Coordinates of the point at which to compute the Christoffel symbols.
    model :
        Function that returns the metric at a given coordinate.
    
    Returns
    -------
    christoffel : jnp.array
        Christoffel symbols at the given coordinate. The first index is upper, the rest two indices are lower.
    """

    met = metric(coords, model)
    pd = jnp.einsum('mns -> smn', pd_metric(coords, model))
    rhs = pd + jnp.einsum('nrm -> mnr', pd) - jnp.einsum('rmn -> mnr', pd)
    rhs = 0.5 * rhs  # shape (dim, dim, dim)
    rhs =  jnp.einsum('rmn -> nrm', rhs)
    # Solve g_{rs} X^s_{mn} = rhs_{r mn}
    # reshape to solve all RHS at once
    dim = met.shape[0]
    rhs_flat = rhs.reshape(dim, -1)

    sol_flat = solve(met, rhs_flat, assume_a='gen')
    christ = sol_flat.reshape(dim, dim, dim)

    return christ

# partial derivatives of the christoffel symbols
pd_christoffel = jax.jit(jax.jacfwd(christoffel), static_argnames=('model',))

@partial(jax.jit, static_argnames=('model',))
def riemann_curvature(coords, model):
    """
    Compute the Riemann curvature at a given coordinate.

    Parameters
    ----------
    coords : jnp.array
        Coordinates of the point at which to compute the Riemann curvature.
    model :
        Function that returns the metric at a given coordinate.
    
    Returns
    -------
    riemann_curvature : jnp.array
        Riemann curvature at the given coordinate. The first index is upper, the rest three indices are lower.
    """

    christ = christoffel(coords, model)
    pd_christ = jnp.einsum('rmns -> srmn', pd_christoffel(coords, model))
    return jnp.einsum('mrns -> rsmn', pd_christ) - jnp.einsum('nrms -> rsmn', pd_christ) + jnp.einsum('rml, lns -> rsmn', christ, christ) - jnp.einsum('rnl, lms -> rsmn', christ, christ)

@partial(jax.jit, static_argnames=('model',))
def ricci_tensor(coords, model):
    """
    Compute the Ricci tensor at a given coordinate.

    Parameters
    ----------
    coords : jnp.array
        Coordinates of the point at which to compute the Ricci tensor.
    model :
        Function that returns the metric at a given coordinate.
    
    Returns
    -------
    ricci_tensor : jnp.array
        Ricci tensor (two lower indices) at the given coordinate. 
    """

    riemann = riemann_curvature(coords, model)
    return jnp.einsum('rsru -> su', riemann)

@partial(jax.jit, static_argnames=('model',))
def ricci_scalar(coords, model):
    """
    Compute the Ricci scalar at a given coordinate.

    Parameters
    ----------
    coords : jnp.array
        Coordinates of the point at which to compute the Ricci scalar.
    model :
        Function that returns the metric at a given coordinate.
    
    Returns
    -------
    ricci_scalar : jnp.array
        Ricci scalar at the given coordinate. 
    """

    return jnp.einsum('mn, mn -> ', jnp.linalg.inv(metric(coords, model)), ricci_tensor(coords, model))

@partial(jax.jit, static_argnames=('model',))
def einstein_tensor(coords, model):
    """
    Compute the Einstein tensor at a given coordinate.

    Parameters
    ----------
    coords : jnp.array
        Coordinates of the point at which to compute the Einstein tensor.
    model :
        Function that returns the metric at a given coordinate.
    
    Returns
    -------
    einstein_tensor : jnp.array
        Einstein tensor (two lower indices) at the given coordinate. 
    """
    met = metric(coords, model)
    ricci_ts = ricci_tensor(coords, model)
    return ricci_ts - 0.5 * jnp.einsum('mn, mn -> ', jnp.linalg.inv(met), ricci_ts).reshape(1, 1) * met

def Kerr_metric(a,Q,rs):
    m=rs/2
    @jax.jit
    def Kerr(coords):
        t, x, y, z = coords[0],coords[1],coords[2],coords[3]
        r=x
        delta=r**2-2*m*r+a**2+Q**2
        rho2=r**2+(a**2)*(jnp.cos(y)**2)
        return jnp.array([
            [-a**2*jnp.sin(y)**2/rho2+delta/rho2,0,0,(a**3*jnp.sin(y)**2+a*r**2*jnp.sin(y)**2-delta*a*jnp.sin(y)**2)/rho2],
            [0,-rho2/delta,0,0],
            [0,0,-rho2,0],
            [(a**3*jnp.sin(y)**2+a*r**2*jnp.sin(y)**2-delta*a*jnp.sin(y)**2)/rho2,0,0,(delta*a**2*jnp.sin(y)**4-a**4*jnp.sin(y)**2-2*a**2*r**2*jnp.sin(y)**2-r**4*jnp.sin(y)**2)/rho2]
        ])
    return Kerr
