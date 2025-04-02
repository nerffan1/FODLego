# FODLego
This program is meant to predict the FOD Positions of molecules for use in 
FLOSIC calculations.

You can see the program in action in the following [section](https://nerffan1.github.io/FLOSICDOC/Reference/Generation.html) of the FLOSIC documentation.

## Implementation Description:
Currently the bond orders are being predicted using RDKit. The Bonding FODs 
(BFODs) and Free FODs (FFODs) are selected based of monoatomic calculations and based of the placement of BFODs. A few rules are based of empirical 
rules.

## Prerequisites:
In order to run, please install the following packages in python 

```
$pip install rdkit scipy
```

# Contact:
Please contact me at my university email ville2a@cmich.edu
