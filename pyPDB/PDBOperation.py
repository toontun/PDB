#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 17 09:38:20 2020

@author: bellaiche
"""
import numpy as np
from scipy.spatial import distance

class PDBOperation:
    
    def __init__(self, pdbpath):
        ##EXCEPTION A GERER AVEC LINE.SPLIT, en effet des fois on n'a pas d'espace
        ##entre temperature factor et occupency
        self.text=""
        self.pdb = pdbpath
        self.coordinates = {}
        self.atom_type = {}
        self.res_type = {}
        self.chains = {}
        self.res_number = {}
        self.occupency = {}
        self.temperature_factors = {}
        self.element_symbol = {}
        with open(pdbpath, "r") as PDB:
            for line in PDB:
                list_line=line.split()
                id = list_line[0]
                if "ANISOU" not in id:
                    self.text = self.text + line
                if id == "ATOM":
                    atom_number = int(list_line[1])
                    self.coordinates[atom_number] =\
                    np.array([float(list_line[6]),
                        float(list_line[7]),
                        float(list_line[8])])
                    self.atom_type[atom_number] = list_line[2]
                    self.res_type[atom_number] = list_line[3]
                    self.chains[atom_number] = list_line[4]
                    self.res_number[atom_number] = list_line[5]
                    self.occupency[atom_number] = list_line[9]
                    self.temperature_factors[atom_number] = list_line[10]
                    self.element_symbol[atom_number] = list_line[11]

        self.matroordinates = np.zeros((len(self.coordinates),3))
        i = 0
        for k,v in self.coordinates.items():
            self.matroordinates[i] = v
            i+=1

    def write_new_pdb(self):
        with open(self.pdb, "w") as f:
            f.write(self.text)

    def print_coordinates(self):
        for atom in self.coordinates:
            print("{} : {}".format(atom, self.coordinates[atom]))

    def translate(self, xtranslate=0, ytranslate=0, ztranslate = 0):
        translation_array = np.array([xtranslate, ytranslate, ztranslate])
        i = 0
        for k, v in self.coordinates.items():
            self.coordinates[k] += translation_array
            self.matroordinates[i] += translation_array
            i+=1

    def get_min_distance_between_PDB(self, pdb):
        distance_matrix = distance.cdist(self.matroordinates, pdb.matroordinates)
        return np.amin(distance_matrix)



