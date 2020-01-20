#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 17 09:27:01 2020

@author: bellaiche
"""

import sys
import PDBOperation as pdbo

if __name__ == "__main__":

    pdb1 = pdbo.PDBOperation(sys.argv[1])
    #pdb1 = pdb1.write_new_pdb()
    # pdb1.print_coordinates()
    for atom in pdb1.atoms:
    	print(str(atom) + "    " + pdb1.atoms[atom])
