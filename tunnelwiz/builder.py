import MDAnalysis as mda
import numpy as np
import math
from scipy.spatial import KDTree


def get_axis(universe: mda.Universe, endpoints: list, type: str ="resid"):
    
    """
    This function creates an axis through a biomolecular tunnel with specified origins, it returns axis sample points, tangent, normal and binormal vectors.
    The search of the axis is based on a Voronoi diagram with weight of the ridges altered by the proximity to the nearest atoms. Dijsktras algorithm then finds the shortes/least resistant path.

    :param universe:  biomolecule loaded with mda
    :param endpoints: list of: resids/atomids/[cartesian coordinates] which will be the descriptors of then tunnel exits

    """
    
    from scipy.sparse import csr_matrix
    from scipy.sparse.csgraph import dijkstra
    from scipy.spatial import Voronoi
    from scipy.interpolate import make_splprep
    from itertools import combinations

    if type == "resid":
        endpoints = [universe.select_atoms(f"resid {endpoints[i]}").center_of_mass() for i in range(2)]
    if type == "atomid":
        endpoints = universe.atoms[endpoints].positions
    if type == "xyz":
        endpoints = endpoints

    # Reduce the number of atoms used in upcoming calculations
    radius=20
    n_samples=10
    tunnel_axis = endpoints[1]-endpoints[0]
    tunnel_samples = [endpoints[0]+tunnel_axis*n/(n_samples-1) for n in range(n_samples)]
    tunnel_samples = np.asarray(tunnel_samples, dtype=float).reshape(-1, 3)
    selections = []
    for n in range(len(tunnel_samples)):
        sel = 'point '+' '.join([str(tunnel_samples[n][i]) for i in range(3)])+f' {radius}'
        selections.append(universe.select_atoms(sel))
    
    index = [*set(np.concatenate([selections[i].atoms.indices for i in range(n_samples)]))]
    u_sel =  universe.atoms[index]

    # 
    tree = KDTree(u_sel.atoms.positions)
    vor = Voronoi(u_sel.atoms.positions)
    V = vor.vertices
    edge_w = {}
    
    # Cost function to alter weights og ridges from Voronoi diagram
    def cost_function(A: list,B: list,universe,alpha: float=10,sigma: float=1.2,r0: float=3):
        dist = np.linalg.norm(B-A)
        nn_d = [tree.query(pos,k=1)[0] for pos in [A,B,(A+B)/2]]
        delta = np.atleast_1d(nn_d).min() - r0
        cost = dist*(1+ alpha*math.exp(-delta/sigma))
        return cost
    
    for rv in vor.ridge_vertices:
            # Remove infinite vertices
            rv = [k for k in rv if k != -1]
            if len(rv) < 2:
                continue
    
            # Connect all pairs on a given ridge face
            for i, j in combinations(rv, 2):
                if i == j:
                    continue
                a, b = (i, j) if i < j else (j, i)
                w = cost_function(V[a],V[b],u_sel)
                prev = edge_w.get((a, b))
                if prev is None or w < prev:
                    edge_w[(a, b)] = w

    # Rearange data to a sparse matrix and find shortest/least resistant path
    tree_vor = KDTree(V)
    idx_start = tree_vor.query(endpoints[0])[1]
    idx_end   = tree_vor.query(endpoints[1])[1]    
    rows, cols, data = [], [], []
    for (a, b), w in edge_w.items():
        rows += [a]
        cols += [b]
        data += [w]
    A = csr_matrix((data, (rows, cols)), shape=(len(V), len(V)))
    
    dist, pred = dijkstra(A, directed=False, indices=idx_start, return_predecessors=True)
    
    # Reconstruct the path
    path = []
    cur = idx_end
    while cur != idx_start:
        path.append(cur)
        cur = int(pred[cur])
    path.append(idx_start)
    path.reverse()
    tunnel_path = vor.vertices[path]

    # Smooth and sample the axis
    n_points = 100
    smooth = 0.0
    path = np.linalg.norm(tunnel_path[1:] - tunnel_path[:-1], axis=1)
    path_cs = np.concatenate(([0.0], np.cumsum(path)))
    param = path_cs/path_cs[-1]
    xyz = np.linspace(0,1,n_points)
    
    # Fit parametric smoothing spline
    fit, u = make_splprep(tunnel_path.T,u=param, s=smooth)
    axis = np.vstack(fit(xyz).T)

    # Tangent vectors
    seg = axis[1:] - axis[:-1]
    seg_len = np.linalg.norm(seg, axis=1, keepdims=True)
    t = seg / np.maximum(seg_len, 1e-12)

    # Choose initial normal vector
    basis = [[1,0,0],[0,1,0],[0,0,1]]
    x = np.asarray([np.dot(v,t[0]) for v in basis])
    for i in range(3):
        if abs(x)[i] == min(abs(x)):
            a = i
    n0 = basis[a] - np.dot(t[0],basis[a])*t[0]
    b0 = np.cross(n0,t[0])
    n = np.zeros((len(t),3), float) 
    b = np.zeros((len(t),3), float)
    n[0] = n0
    b[0] = b0

    # Normal and binormal vectors using Rotation Minimizing Frames with Euler-Rodrigues' formula
    for i in range(1,len(t)):
        ti_p = t[i-1]
        ti = t[i]
        ni_p = n[i-1]
        
        v = np.cross(ti_p,ti)
        s = np.linalg.norm(v)
        c = np.clip(np.dot(ti,ti_p), -1.0, 1.0)
        
        k = v/s
        angle = np.arctan2(s,c)
        ni = (ni_p*np.cos(angle) + np.cross(k, ni_p)*np.sin(angle) +k*np.dot(k,ni_p)*(1-np.cos(angle)))
    
        ni = ni - np.dot(ni,ti)*ti
        ni = ni/np.linalg.norm(ni)
        bi = np.cross(ti,ni)
        
        n[i] = ni
        b[i] = bi

    return axis,t,n,b


def get_atoms(universe: mda.Universe, axes: np.array, beta: float=4.0):

    """
    This function realizes atom selection along the tunnel axis and returns a universe object
    
    :param universe: biomolecule loaded with mda
    :param axis: tunnel axis specified by sample points, tangent and normal vectors of the type [ [N,3], [N,3], [N,3] ], where N is the number of sample points
    """

    axis = axes[0]
    tree = KDTree(universe.atoms.positions)

    # Helper function - Fit spheres for a list of points from a given center
    def fit_sphere(points: list, center: list):
        P = np.asarray(points, float)
        c = np.asarray(center, float)
        d = np.linalg.norm(P - c[None, :], axis=1)
        r = np.sqrt((d * d).mean())
        return r
    sphere_radii = [fit_sphere(universe.atoms[tree.query(axis[i], k=4)[1]].positions, axis[i]) for i in range(len(axis))]

    # Select atoms around fitted spherese with an an extra radius of beta
    # A better selection scheme planned for future
    selections = []

    for n in range(len(axis)):
        sel = 'point '+' '.join([str(axis[n][i]) for i in range(3)])+' '+str(sphere_radii[n]+beta)
        selections.append(universe.select_atoms(sel))
    
    tunnel_atoms = universe.atoms[list(set(np.concatenate(
        [selections[i].atoms.ix for i in range(len(selections))]
    )))]
    
    return tunnel_atoms


def get_projection(points: list, axes: np.array):
    """
    This function creates a projection from cartesian coordinates to cylindrical coordinates relative to each axis frame, datapoints are assigned to a given frame based on distances
    Returns coordinates: [Z, Theta, R], (distance along axis, angle from a normal vector. radial distance from axis)
    
    :param points: list of coordinates of a shape (N,3)
    :param axis: tunnel axis specified by sample points, tangent and normal vectors of the shape ((N,3),4) where N is the number of sample points
    """

    axis = axes[0]
    t = axes[1]
    n = axes[2]
    b = axes[3]
    K = t.shape[0]
    N = points.shape[0]
    Axis = axis[:-1]
    tt = np.sum(t*t, axis=1)
    L = np.sqrt(tt)
    S = np.concatenate(([0,0],np.cumsum(L)))

    # Assign points to the best fitting segment of the axis
    best_d2 = np.full(N, np.inf)           #squared distance of the atom-tunnel path vertex
    best_seg = np.zeros(N, dtype=np.int32) #which segment best fits given atom
    best_u = np.zeros(N, dtype=float)      #length along found segment (0,1)
    best_Q = np.zeros((N,3), dtype=float)  #tunnel_path vertices

    for i in range(K-1):
        Ai=Axis[i]
        ti=t[i]
        denom = tt[i] if tt[i] > 0 else 1.0
        
        AP = points - Ai
        ui = (AP@ti)/denom
        ui = np.clip(ui,0,1)
        Qi = Ai + ui[:, None]*ti
        d2 = np.sum((points-Qi)**2, axis=1)
    
        mask = d2 < best_d2
        if np.any(mask):
            best_d2[mask] = d2[mask]
            best_seg[mask] = i
            best_u[mask] = ui[mask]
            best_Q[mask] = Qi[mask]
    r = np.sqrt(best_d2)
    s = S[best_seg] + best_u*L[best_seg]

    idx_right = np.searchsorted(S, s, side="left")
    idx_left  = np.clip(idx_right - 1, 0, len(S)-1)
    idx_right = np.clip(idx_right,     0, len(S)-1)

    vertex_id = np.where(np.abs(S[idx_right] - s) < np.abs(S[idx_left] - s),idx_right, idx_left).astype(np.int32)
    
    # Transform coordinates
    t_v = np.zeros((K+1,3), dtype=float)
    n_v = np.zeros((K+1,3), dtype=float)
    b_v = np.zeros((K+1,3), dtype=float)
    t_v[0], n_v[0], b_v[0] = t[0], n[0], b[0]
    t_v[-1], n_v[-1], b_v[-1] = t[-1], n[-1], b[-1]
    for i in range(1, K):
        tv = t[i-1] + t[i]
        tv /= np.linalg.norm(tv)
        nv = n[i-1]
        nv -= np.dot(nv, tv) * tv
        nv /= np.linalg.norm(nv)
        bv = np.cross(tv, nv)
        t_v[i], n_v[i], b_v[i] = tv, nv, bv
    vid = vertex_id
    tv = t_v[vid]
    nv = n_v[vid]
    bv = b_v[vid]

    w = points - best_Q 
    w_perp = w - (np.sum(w*tv, axis=1)[:, None])*tv
    x = np.sum(w_perp*nv, axis=1)
    y = np.sum(w_perp*bv, axis=1)
    theta = np.arctan2(y,x)
    
    return s, theta, r


def construct_tunnel(universe: mda.Universe, endpoints: list, type: str ="resid"):
    
    axis = np.asarray(get_axis(universe, endpoints, type))
    tunnel_atoms = get_atoms(universe, axis)
    cyl_coor = get_projection(tunnel_atoms.positions, axis)

    return tunnel_atoms, cyl_coor
