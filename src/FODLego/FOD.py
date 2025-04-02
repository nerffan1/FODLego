import numpy as np
from typing import List
from numpy.linalg import norm

class FOD:
    def __init__(self, pos: np.ndarray = np.zeros(3), ch=True) -> None:
        self.mPos = pos
        self.mSiblings = ()
        self.mAssocFOD: FOD = []
        self.mChannel = ch

    def AddSibling(self, *siblings):
        """
        Add associated FODs, siblings, to the FOD. For example, the 3 FODs in a TBFOD would be the siblings, where each points to the other 2.
        """
        self.mSiblings = siblings

    #################### Operator Overloading #################### 
    def __str__(self) -> str:
        return str(self.mPos)

    def __mul__(self, factor: float):
        self.mPos *= factor
        return self

    def __add__(self, shift: float):
        """
        Overloading addition operation to add shift the mPos value. Mutates instance.
        """
        return (self.mPos + shift)

    def __sub__(self, shift: float):
        """
        Overloading addition operation to add a negative shift the mPos value. Mutates instance.
        """
        return (self.mPos - shift)

class CFOD(FOD):
    """
    This class is essentially the same as <FOD>. For categorization purposes, it was better to
    create this class instead of creating a Boolean member inside the FOD class that made it a
    CFOD.
    """
    def __init__(self, atom, pos: np.ndarray = np.zeros(3), ch=True) -> None:
        super().__init__(pos,ch)
        self.mAtom = atom
        self.mR = np.linalg.norm(self.mPos - atom.mPos)
