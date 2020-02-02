# workflows
This is a playground for me aiming at setting up my own workflow
to automatize a gimic integral calculation over a bond. 

You need in the same directory the following files:

coord.xyz
get_int.py
func.py
XDENS
MOL

Usage of the script is: 
python get_int.py --atoms X,Y --fix Z 

XY are the atom indices of the bond you are interested in
Z  is the atom index of the fixpoint which should be on the 
left side when the bond XY is oriented to the right.
