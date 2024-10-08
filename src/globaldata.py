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

class GlobalData:
    
    def __init__(self):
        GlobalData.mElementInfo = self.LoadElements()
        GlobalData.mElemNames = self.LoadNames()
    
    #Public Functions
    def LoadElements(self):
        """
        Loads information from CSV table. This function can probably be deprecated by
        replacing the information offered by this table using RDKit.
        """
        # Get the path to src directory to access elements2.0 file
        pth = path.dirname(argv[0])
        with open(f'{pth}/elements2.0', mode ='r') as file:
            return genfromtxt(file, delimiter=',', encoding=None, dtype=None)

    def LoadNames(self):
        """
        Names are loaded in their abbreviated form, since XYZ files receive them this way
        """
        return GlobalData.mElementInfo[1:,2]

    def GetElementAtt(name: str, attr: str):
        """
        Gets the attribute found in the elements2.0 file

        Parameters:
        name (str): The name of the atom, in abbreviated form.
        attr (str): The attribute you want of chosen atom, found in the first row of elements2.0

        """
        if attr == "AtomicNumber":
            attrib_i = where(GlobalData.mElemNames == name)  
            # We must offset index by +1 because the array starts at zero instead of 1
            return attrib_i[0].item() + 1      
        else:
            element_i = where(GlobalData.mElemNames == name)[0]
            attrib_i = where(GlobalData.mElementInfo[0] == attr)[0]
            #We offset by 1 because the Info list has attribute names on row 1
            return GlobalData.mElementInfo[element_i + 1,attrib_i]
        
    def GetZAtt(Z: int, attr: str):
        attrib_i = where(GlobalData.mElementInfo[0] == attr)
        #We offset by 1 because the Info list has attribute names on row 1
        return GlobalData.mElementInfo[Z,attrib_i[0].item()]

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

    #Degugging tests
    def _debug_samplenames():
        for att in ["AtomicNumber", "Group", "Period", "Element","Metal", "NumberofShells","NumberofValence"]:
            print(f'{att}: {GlobalData.GetElementAtt("Ga", att)}')
            print(GlobalData.GetZAtt(31, att))

    ############Class Variables############    
    mElementInfo = []
    mElemNames = []
    mClosedGroups = [2,12,18]
    mShellCount = [0,2,8,8,18,18]

    # The following ladder is based of various monoatomic calculations.
    #Think of a scheme that places the beginning of a 
    mGeo_Ladder = { 2: ['point'], 
                        4: ['point','point'],
                        10: ['point','tetra'], 
                        18: ['point', 'tetra', 'tetra'],
                        20: ['point', 'tetra', 'triaug_val', 'point'], 
                        30: ['point', 'tetra', 'triaug', 'point'],
                        36: ['point', 'tetra', 'triaug', 'tetra'],
                        54: ['point', 'tetra', 'triaug', 'triaug', 'tetra'] }
    mShellShapes = {1: ['point'], 4: ['tetra'], 9: ['triaugmented']}
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
