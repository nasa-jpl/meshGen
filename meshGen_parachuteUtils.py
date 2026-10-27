# Definitions
EPS = 1.0e-14



def unit_vector (a):
    """ Convert vector a to unit vector
    """
    import numpy as np
    tol = 1.0e-16
    norm = np.linalg.norm (a)
    if norm <= tol:
        return a
    else:
        return a / norm 


   
    
def createRotationMatrix (a, b):
    """ 
    Create a rotation matrix that rotates vector a onto vector b
    - Based on method described here: https://math.stackexchange.com/a/476311
    """
    import numpy as np
    
    # Convert input vectors to numpy arrays
    a = np.array (a)
    b = np.array (b)
    
    # Normalize vectors
    aNorm = a / np.linalg.norm (a)
    bNorm = b / np.linalg.norm (b)
  
    # Check if vectors are parallel (if so, cross product will be zero)
    if np.linalg.norm (np.cross (aNorm, bNorm)) <= EPS:
        if np.dot (aNorm, bNorm) >= 0.0:
            return np.identity(3)       # Return identity -> parallel and pointing in same direction
        else:
            return -np.identity(3)      # Return -identity -> parallel and pointing in opposite direction
            
    # If not parallel, calculate rotation matrix
    else:
        v = np.cross (a, b)                         # Get vector perpendicular to a and b
        s = np.linalg.norm (v)                      # Sin of angle between vectors
        c = np.dot (a, b)                           # Cos of angle between vectors
        v_x = np.array([[ 0.0,  -v[2],  v[1]],      # Skew-symmetric cross product matrix of v   
                        [ v[2],  0.0,  -v[0]],
                        [-v[1],  v[0],  0.0]])
        R = np.identity(3) + v_x + np.matmul (v_x, v_x) * (1.0 - c * c) / (s * s)
    return R    
  



def calcTriangleNormal (node1, node2, node3):
    """ Calculate normal vector for triangle defined by 3 nodes
    """
    import numpy as np
    
    # Calc vectors and normal
    v1 = np.array(node2) - np.array(node1)
    v2 = np.array(node3) - np.array(node1)
    normal = np.cross(v1, v2)
    
    # Normalize vector and return result
    norm = np.linalg.norm(normal)
    if norm <= EPS:
        raise ValueError("Error in calcTriangleNormal: two sides are collinear so normal vector is undefined")
    return normal / norm
   



def calcQuadNormal (node1, node2, node3, node4):
    """ Calculate normal vector for quadrangle defined by 4 nodes
    """
    import numpy as np
    
    # Calc vectors and normal
    v1 = np.array(node3) - np.array(node1)
    v2 = np.array(node4) - np.array(node2)
    normal = np.cross(v1, v2)
    
    # Normalize vector and return result
    norm = np.linalg.norm(normal)
    if norm <= EPS:
        raise ValueError("Error in calcQuadrangleNormal: two sides are collinear so normal vector is undefined")
    return normal / norm   
    
    
    
def calcVectorRotationAngles (v1, v2, type, flagDegrees=False):
    """
    Calculates parameters to rotate vector v1 to v2
    Inputs:
        v1, v2:         vectors 1 and 2 for which to calculate rotation angles
        flagDegrees:    return results in degrees or rad (default)
        type:           type of result: rotVec (axis-angle representation) or euler (Euler angles)
    """
    import numpy as np
    from scipy.spatial.transform import Rotation as R

    # Normalize vectors
    v1 = v1 / np.linalg.norm(v1)
    v2 = v2 / np.linalg.norm(v2)

    # Get rotation aligning v1 to v2
    rotation, _ = R.align_vectors([v2], [v1])

    # Extract Euler angles (XYZ convention)
    if type == "rotVec":
        return rotation.as_rotvec(degrees=flagDegrees)
    elif type == "euler":
        return rotation.as_euler('xyz', degrees=flagDegrees)
    else:
        raise Exception ("Error in calcVectorRotationAngles: provide correct type input for result format")