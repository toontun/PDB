#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 17 09:38:20 2020

@author: bellaiche
"""

class PDBOperation:
    
    def __init__(self, pdbpath):
        self.text=""
        self.pdb = pdbpath
        self.atoms = {}
        with open(pdbpath, "r") as PDB:
            for line in PDB:
                list_line=line.split()
                id = list_line[0]
                if "ANISOU" not in id:
                    self.text = self.text + line
                if id[0] == "ATOM":
                    self.atoms[int(list_line[5])] = {list[6:8]}
    
    def write_new_pdb(self):
        with open(self.pdb, "w") as PDB:
            PDB.write(self.text)
    
    def rename_chains(self):
        #alphabet ["A", "B", "C", "D", ]
        pass

    def print_coordinates(self):
        for atom in self.atoms:
            print(atom)
