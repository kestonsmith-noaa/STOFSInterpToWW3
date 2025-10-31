import numpy as np
import time
from datetime import datetime

def IsInElement(x,y,xp,yp):
    IsIn=False
    c1 = (x[1] - x[0]) * (yp - y[0]) - (y[1] - y[0]) * (xp - x[0])
    c2 = (x[2] - x[1]) * (yp - y[1]) - (y[2] - y[1]) * (xp - x[1])
    c3 = (x[0] - x[2]) * (yp - y[2]) - (y[0] - y[2]) * (xp - x[2])
    if ( ( c1 > 0 and c2 > 0 and c3 > 0) or ( c1 < 0 and c2 < 0 and c3 < 0) ):
        IsIn=True
    if (c1*c2*c3 == 0): # include points on triangle
        IsIn=True
    return IsIn

def FindElement(x,y,e,xi,yi):
    nd=len(e.shape)
    e=np.squeeze(e)
    j=-9999
    ne=e.shape[0]
    #print(nd)
    #print(ne)
    #print(e)
    if nd > 1: 
        xc = np.squeeze(np.mean(x[e-1], axis=1))
        yc = np.squeeze(np.mean(y[e-1], axis=1))
        DistanceToElements = np.abs((xi + 1j*yi) - (xc + 1j*yc))
        n=0
        while (j < 0 and n < ne)  :
            n=n+1
            j = np.argmin( DistanceToElements )
            xl = np.squeeze(x[e[j,:]-1])
            yl = np.squeeze(y[e[j,:]-1])
            IsIn=IsInElement(xl,yl,xi,yi)
            if IsIn :
                return j
            else:
                DistanceToElements[j]=float('inf')
                j=-9999
    else:
        xl = np.squeeze(x[e-1])
        yl = np.squeeze(x[e-1])
        IsIn=IsInElement(xl,yl,xi,yi)
        if IsIn:
            j=0
            return j
        else:
            j=-9999
            return j
        
    return j

def compute_mesh_to_mesh_interp_weights(x, y, e, xi, yi, flout, MaxDist):
    """
    Compute mesh-to-mesh interpolation weights using barycentric coordinates.
    
    Parameters:
    -----------
    x : array-like
        X coordinates of source mesh nodes
    y : array-like
        Y coordinates of source mesh nodes
    e : array-like
        Element connectivity matrix (n_elements x 3 for triangular elements)
    xi : array-like
        X coordinates of target points
    yi : array-like
        Y coordinates of target points
    flout : str
        Output filename for saving results
        
    Returns:
    --------
    weights : ndarray
        Interpolation weights (n_points x 3)
    nodes : ndarray
        Node indices for each point (n_points x 3)
    elenum : ndarray
        Element number for each point (n_points,)
    Dist2EleCenter : ndarray
        Distance to element center for each element
    """
    
    # Convert to numpy arrays
#    MaxDist=0.25 # maximum distance to use in interpolation
    
    x = np.asarray(x)
    y = np.asarray(y)
    e = np.asarray(e)
    xi = np.asarray(xi)
    yi = np.asarray(yi)
    
    # Compute element centers (mean of vertices)
    # MATLAB: mean(x(e'))' with 1-based indexing
    # Python: use 0-based indexing
    xc = np.mean(x[e-1], axis=1)
    yc = np.mean(y[e-1], axis=1)

    # Initialize output arrays
    n_points = len(xi)
    n_elements = len(e)
    weights = np.zeros((n_points, 3))
    nodes = np.zeros((n_points, 3), dtype=int)
    elenum = np.zeros(n_points, dtype=int)
    Dist2EleCenter = np.zeros(n_elements)
    
    t0 = time.time()
    UseNearestEle = False
    
    with open(flout+".csv", 'w') as f:    
    
        for k in range(n_points):
            #print("on point number: " + str(k))
            # Find nearest element center using complex number distance
            # MATLAB: abs(xi(k)+i*yi(k)-xc-i*yc)
            jg=np.where( np.all([( np.abs(xi[k]-xc) < MaxDist ), (np.abs(yi[k]-yc) < MaxDist )], axis=0) )
            jg=np.squeeze(jg)
#            if jg.shape[0]==0 :
     #       print(jg)
     #       print(e[jg,:])
            if jg is None: 
                print("No elements within Maxdist, using closest element")
                distances = np.abs((xi[k] + 1j*yi[k]) - (xc + 1j*yc))
                j = np.argmin(distances)
            else:
                if UseNearestEle :
                    distances = np.abs((xi[k] + 1j*yi[k]) - (xc[jg] + 1j*yc[jg]))
                    j = np.argmin(distances)
                    j=jg[j]
                else: 
                    j=FindElement(x,y,e[jg,:],xi[k],yi[k])
                    if j < 0:
                        print("couldn't find element for point : "+ str(k))
                        distances = np.abs((xi[k] + 1j*yi[k]) - (xc + 1j*yc))
                        j = np.argmin(distances)
                        print("using closest element at distance : "+ str(distances[j]))
                        print(" ")
                    else:
                        j=jg[j]

            elenum[k] = j
            Dist2EleCenter[k] = np.abs((xi[k] + 1j*yi[k]) - (xc[j] + 1j*yc[j]))
            #print("distance: " + str(distances[j]) + ", elem: "+str(j) )
            
            x0 = xi[k]
            y0 = yi[k]
            
            # Get vertices of the element
            xl = x[e[j, :]-1]
            yl = y[e[j, :]-1]
            
            # Compute barycentric coordinates using shoelace formula for areas
            # Area of triangle formed by vertices 1, 2, and point (a3)
            xt = np.array([xl[0], xl[1], x0, xl[0]])
            yt = np.array([yl[0], yl[1], y0, yl[0]])
#           a3 = -np.dot(xt[1:3] - xt[0:2], yt[0:2] + yt[1:3]) / 2
            a3 = -np.dot(xt[1:4] - xt[0:3], yt[0:3] + yt[1:4]) / 2
            
            # Area of triangle formed by vertices 3, 1, and point (a2)
            xt = np.array([xl[2], xl[0], x0, xl[2]])
            yt = np.array([yl[2], yl[0], y0, yl[2]])
#            a2 = -np.dot(xt[1:3] - xt[0:2], yt[0:2] + yt[1:3]) / 2
            a2 = -np.dot(xt[1:4] - xt[0:3], yt[0:3] + yt[1:4]) / 2
            
            # Area of triangle formed by point, vertices 2, 3 (a1)
            xt = np.array([x0, xl[1], xl[2], x0])
            yt = np.array([y0, yl[1], yl[2], y0])
#            a1 = -np.dot(xt[1:3] - xt[0:2], yt[0:2] + yt[1:3]) / 2
            a1 = -np.dot(xt[1:4] - xt[0:3], yt[0:3] + yt[1:4]) / 2
            
            # Normalize to get barycentric weights
            total_area = a1 + a2 + a3
            weights[k, :] = [a1, a2, a3] / total_area
            if (Dist2EleCenter[k]>MaxDist):  # Extrapolate values to 0 if outside
                weights[k, :] = [0., 0., 0.] # maximum interpolation distance

            nodes[k, :] = e[j, :]
            lineN = ', '.join(str(val) for val in nodes[k, :])
            lineW = ', '.join(str(val) for val in weights[k, :])
            f.write(lineN + ", " + lineW + '\n')
            # Progress reporting every 100 iterations
            if (k + 1) % 10 == 0:
                t1 = time.time()
                time_per_iter = (t1 - t0) / (k + 1)
                time_remaining = (n_points - k - 1) * time_per_iter / 60
                print(f"Progress: {k+1}/{n_points}")
                print(f"Time remaining: {time_remaining:.2f} minutes")

    # Save results to file (using NumPy's format) 
# These files are O(100) times the size of the ascii files with same data 
#        np.savez(flout, 
#             weights=weights, 
#             nodes=nodes, 
#             elenum=elenum, 
#             Dist2EleCenter=Dist2EleCenter)
    
    print(f"Results saved to {flout}.npz")
    
    return weights, nodes, elenum, Dist2EleCenter

def InterpolateField2Nodes(nodes,weights, f):
    fi=np.zeros(weights.shape[0])
    for k in range(weights.shape[0]):
        fl=f[nodes[k,:]]
        wl=weights[k,:]
        fi[k]=np.dot(wl,fl)
        if (not np.abs(fi[k]) > 0.):
           fi[k]=0.
        if ( fi[k]>np.max(fl)  ):
           fi[k]=np.max(fl)
        if ( fi[k]<np.min(fl)  ):
           fi[k]=np.min(fl)


    return fi


    """

# Example usage:
if __name__ == "__main__":
    # Example with a simple triangular mesh
    # Source mesh nodes
    x = np.array([0, 1, 0.5, 1.5])
    y = np.array([0, 0, 1, 1])
    
    # Element connectivity (0-based indexing)
    e = np.array([[0, 1, 2],
                  [1, 3, 2]])
    
    # Target points
    xi = np.array([0.5, 0.8])
    yi = np.array([0.3, 0.5])
    
    # Compute weights
    weights, nodes, elenum, dist = compute_mesh_to_mesh_interp_weights(
        x, y, e, xi, yi, 'output_data'
    )
    
    print("\nResults:")
    print(f"Weights:\n{weights}")
    print(f"Nodes:\n{nodes}")
    print(f"Element numbers: {elenum}")
"""
