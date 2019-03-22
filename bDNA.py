import numpy as np

def angle(v1, v2):
    """ calculate angle 
        input: vector1 (3D numpy array), vector2 (3D numpy array)
        out: angle in radians from interval [0,pi] 
    """
    angle = np.arccos(np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2)))
    return angle

def plane_unitnormal_array(points_1, points_2, points_3):
    """ calculate plane unit vector
        in: 3 input points 3D-array
        returns: numpy 3D vector
    """
    nhats = []
    for p1, p2 ,p3 in zip(points_1, points_2, points_3):
        nhat = np.cross((p2-p1), (p3-p1))
        nhat = nhat/np.sqrt( np.dot( nhat,nhat) )
        nhats.append(nhat)
    return nhats


def distance_array(group1, group2):
    
    distances = []
    for i, g1 in enumerate(group1):
        for g2 in group2[i+1:]:
            distances.append(np.linalg.norm(g1-g2))
    
    return distances


def distance_array_normal(group1, nhat, group2):
    
    distances_normal = []
    for i, g1 in enumerate(group1):
        for g2 in group2[i+1:]:
            regular = np.linalg.norm(g1-g2)
            if (regular < 10.) :
                normal = np.abs(np.dot(nhat[i], (g2 - g1)))
                distances_normal.append(normal)
            else:
                distances_normal.append(regular)

    d = np.asarray(distances_normal, dtype=np.float64)
    return d

    # Define model function to be used to fit to the data above:
def gauss(x, A, mu, sig):  
    return A*np.exp(-(x-mu)**2/(2.*sig**2))

def func(x, *p):
    y = np.zeros_like(x)
    for i in range(0, len(p), 3):
        A = p[i]
        mu = p[i+1]
        sig = p[i+2]
        y = y + gauss(x, A, mu, sig)
    return y


