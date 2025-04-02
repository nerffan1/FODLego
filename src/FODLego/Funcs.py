from  typing import List
import numpy as np
from scipy.spatial.transform import Rotation as R

def AddNormals(vectors: list[np.array]) -> np.ndarray:
    """This function normalizes all the given vectors, adds their normal,
    and normalizes the resultant vector.
    
    This function is motivated by the observation that FFODs are closer to the
    bisecting angle in double FFODs, rather than a weighted """
    free_dir_norm = [1/np.linalg.norm(x) for x in vectors]
    free_dir = vectors * np.reshape(free_dir_norm,(len(vectors),1))
    free_dir = free_dir.sum(0)
    free_dir /= np.linalg.norm(free_dir)
    return free_dir
            
def RotateVec(V: np.ndarray, axis: np.ndarray):
    rot = R.from_rotvec(axis)
    return np.matmul(rot.as_matrix(),V)

def RotatePoints(n:int,fod0:np.ndarray,axis:np.ndarray) -> List[np.ndarray]:
    """
    This function creates n points in a circle, starting from your fod0 (i.e. the first fod in the circle).

    n: Number of points to equally distribute on a circle
    ogp: The Original Point, the first element in your circle
    axis: The axis of rotation
    """
    assert len(axis) == 3, "The array must have 3 dimensions"
    fodsRotated = [fod0]
    step = (2*np.pi)/n
    for i in range(1,n):
        rot = R.from_rotvec((step*i)*axis)
        fod = np.matmul(rot.as_matrix(),fod0)
        fodsRotated.append(fod)
    return fodsRotated

def RotateNormals(n: int, firstdir: np.ndarray, axis: np.ndarray) -> List[np.ndarray]:
    # Assertions
    assert len(axis) == 3, "The array must have 3 dimensions"
    assert len(firstdir) == 3, "The array must have 3 dimensions"

    # Variables
    fodsRotated = [firstdir]
    step = (2*np.pi)/n

    # Rotate input vector and create a new coordinate from that
    for i in range(1,n):
        rot = R.from_rotvec((step*i)*axis)
        fod = np.matmul(rot.as_matrix(),firstdir)
        fodsRotated.append(fod)
    return fodsRotated

def RandomPerpDir(ref: np.ndarray) -> np.ndarray:
    """
    This function returns a random perpendicular direction with respect to your reference direction.

    Ref: Reference direction
    """
    if ref[0] == 0: return np.array([1.0,0.0,0.0])
    elif ref[1] == 0: return np.array([0.0,1.0,0.0])
    elif ref[2] == 0: return np.array([0.0,0.0,1.0])
    else:
        b_z = -(10*ref[0] + 2*ref[1])/ref[2]
        randperp = np.array([10,2,b_z])
        return normalize(randperp)     
    
def AngleBetween(A: np.ndarray, B: np.ndarray) -> float:
    """
    Returns the angle between vector A and vector B using simple dot product properties. 
    """
    ab = np.dot(A,B)
    ab_abs = np.linalg.norm(A)*np.linalg.norm(B)
    return np.arccos(ab/(ab_abs))


####LAMBDA FUNCTIONS######
normalize = lambda v: v/np.linalg.norm(v)
tofrom = lambda v, w: v - w
dist = lambda a,b: np.linalg.norm(a-b)
exists = lambda typ,l: any(isinstance(fod, typ) for fod in l)
