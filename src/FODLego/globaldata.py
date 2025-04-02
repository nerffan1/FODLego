#Description: The GlobalData class is dedicated to holding any static information that will be used
#  in the heuristic method of FODs. Some of the data it includes are:
#  - Number of electrons in closed shells
#  - Access to elements information found in the elements2.0 file
#RoadMap: Will eventually require access to average FOD-Tetrahedra radii (along those lines...), so that 
#  we can have a good guess at their positions in space. If there must be hybridization of FODs and Tetrahedra
#  combine, then we can also combine the shared FODs (i.e. the bonding FODs), and use the monoatomic data for reference.
#  Currently more time to look into heuristics is necessary.
#Author: Angel-Emilio Villegas S.
from numpy import where, genfromtxt, sqrt
from os import path
from sys import argv
from numpy import array


class GlobalData:

    #Public Functions
    @staticmethod
    def GetFullElecCount(group: int, period: int):
        """
        Receive Group number and return the amount of electrons that the atoms contains
        when its valence is fully completed. E.g., Z = 6, then we will return 10 since 10
        closes the subshell in that row of the periodic table.
        TODO: Implement the metals     
        TODO: Implement 1s shell logic   
        """
        if group <= 2:
            if period == 1: return 2
            elif period == 2: return 4
            elif period == 3: return 12
            elif period == 4: return 20
            elif period == 5: return 38
            elif period == 6: return 56
            else: -1 
        elif group > 2 and group <= 18: 
            if period == 1: return 2
            elif period == 2: return 10
            elif period == 3 : return 18
            elif period == 4: return 36    
            elif period == 5: return 54
            elif period == 6: return 86

    @staticmethod
    def GetRow(Z: int) -> int:
        row4 = array([22,40,72])
        row13 = array([5,13,31,49,81])
        if Z in [1,3,11,19,37,55]:
            return 1
        elif Z in [4,12,20,38,56]:
            return 2
        elif Z in [21,39]:
            return 3
        elif Z in row4:
            return 4
        elif Z in row4 + 1:
            return 5
        elif Z in row4 + 2:
            return 6
        elif Z in row4 + 3:
            return 7
        elif Z in row4 + 4:
            return 8
        elif Z in row4 + 5:
            return 9
        elif Z in row4 + 5:
            return 10
        elif Z in row4 + 5:
            return 11
        elif Z in row4 + 6:
            return 12
        elif Z in row13:
            return 13
        elif Z in row13 + 1:
            return 14
        elif Z in row13 + 2:
            return 15
        elif Z in row13 + 3:
            return 16
        elif Z in row13 + 4:
            return 17
        elif Z in row13 + 5 or Z==2:
            return 18
        else:
            return -1

    @staticmethod
    def GetPeriod(Z: int) -> int:
       p1 = array([1,2])
       p2 = array([3,10])
       p4 = array([19,36])

       if Z >= p1[0] and Z <= p1[1]: return 1
       elif Z >= p2[0] and Z <= p2[1]: return 2
       elif Z >= p2[0]+8 and Z <= p2[1]+8: return 3
       elif Z >= p4[0] and Z <= p4[1]: return 4
       elif Z >= p4[0]+18 and Z <= p4[1]+18: return 5
       else: return -1

    ############Class Variables############
    mElementInfo = []
    mElemNames = []
    mClosedGroups = [2,12,18]
    mShellCount = [0,2,8,8,18,18]

    # The following ladder is based of various monoatomic calculations.
    #Think of a scheme that places the beginning of a 
    mGeo_Ladder = { 2: ['pt'], #1s
                        4: ['pt','pt'],
                        7: ['pt','tetra', 'tetra'],
                        8: ['pt','tetra','pt'],
                        10: ['pt','tetra'], #1s2s2p 
                        18: ['pt','tetra', 'tetra'], #1s2s2p3s3p
                        20: ['pt', 'tetra', 'triaug', 'pt'], 
                        30: ['pt', 'tetra', 'triaug', 'pt'],
                        36: ['pt', 'tetra', 'triaug', 'tetra'], #1s2s2p3s3p3d4s4p
                        54: ['pt', 'tetra', 'triaug', 'triaug', 'tetra'] }
    mShellShapes = {1: ['point'], 4: ['a_tetra'], 9: ['triaugmented']}
    #This ladder is based of the 
    mElecConfLadder= [2,2,6,2,6,2,10,6,2,10,6,2,10,6]
    #Geometries for known shell structures
    mTriPlane = [[0,sqrt(3)/2,0],
                 [-sqrt(3)/2,-sqrt(3)/2,0],
                 [-sqrt(3)/2,-sqrt(3)/2,0]]
    
    #Molecules
    #Radii of Tetrahedra obtained by closing the shells of 
    #several atoms to close the sp3 shell. E.g., for Z=5, 5 
    # electrons were added. 
    mRadii = {
        10: {
            5: 1.0724622240214408,
            6: 0.9326307459817831,
            7: 0.8245934058016624,
            8: 0.6807571999707002,
            9: 0.6127101178811892, 
            10: 0.6127101178811892,
            11: 0.5005890731191966,
            13: 0.3336349313415564,
            14: 0.29211200538406207,
            15: 0.2619356755640106,
            16: 0.24136512980895986,
            17: 0.2037282349524298,
            18: 0.2037282349524298,
            31: 0.09504603503116242,
            32: 0.09131627498496026,
            33: 0.08749593856607064,
            34: 0.0782336953726361,
            35: 0.08107096633476632,
            36: 0.0782336953726361,
            51: 0.05127235247408605, #TODO: Redo this calculation!
            52: 0.05011655110416162, 
            53: 0.04900681619768362,
            54: 0.04795854985412319
        },
        18: {
            13: 2.9600771375101322,
            14: 1.571379042838154,
            15: 1.179809249943448,
            16: 0.956759529738896,
            17: 0.7103844384774737, 
            18: 0.7103844384774737
        },
        36: { 

        },
        54: {
            51: 3.5394613534573764,
            52: 1.5838395609028515,
            53: 1.5836443998406289,
            54: 1.4501949651948935

        }
    }
    #Average Edge length of FOD-FOD distances in the outmost shell, of this amount of 
    mVert = {
        10: {
            5: 1.751301962745536,
            6: 1.5229774321557852,
            7: 1.3465512547238012,
            8: 1.1116706527406452,
            9: 1.0005509960170376,
            10: 1.0005509960170376,
            11: 0.8174598158123096
        },  
        18: {
            13: 4.173993393062307,
            14: 2.5660512319934154,
            15: 1.9266204374514642,
            16: 1.5623817696036548,
            17: 1.1600529303222396,
            18: 1.1600529303222396
        }
    }
    AU2ANG = 0.529177249
    ANG2AU = 1.8897259886
    mAtoms = []
    mFODs = []
    mBFODs = []
