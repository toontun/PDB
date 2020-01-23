#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 17 09:38:20 2020

@author: bellaiche

Module to mange PDB files: read it, translate coordinates, get distance from
another PDB and write a new pdb file.
"""
import numpy as np
from scipy.spatial import distance

class PDBOperation:
    
    def __init__(self, pdbpath):
        """
        Init an PDBOperation object by parsing a PDB file. Just give the path
        to the PDB file of your choice. Coordinates are also stored in a numpy
        matrix for optimization while distance is calculated. Other data are stored
        in dictionaries to keep the right information on each atoms for further
        manipulations.

        ARGUMENTS:

        pdbpath: the path to the pdb.

        """
        self.text=""
        self.pdb = pdbpath
        self.id = {}
        self.coordinates = {}
        self.atom_type = {}
        self.res_type = {}
        self.chains = {}
        self.res_number = {}
        self.occupency = {}
        self.temperature_factors = {}
        self.element_symbol = {}
        self.atom_charge = {}
        with open(pdbpath, "r") as PDB:
            for line in PDB:
                id = line[0:6].strip()               
                if "ANISOU" not in id:
                    self.text = self.text + line
                if id=="ATOM" or id=="HETATM":
                    atom_number = int(line[6:11].strip())
                    self.id[atom_number]= id
                    self.coordinates[atom_number] =\
                    np.array([float(line[30:38].strip()),
                        float(line[38:46].strip()),
                        float(line[46:54].strip())])
                    self.atom_type[atom_number] = line[12:16].strip()
                    self.res_type[atom_number] = line[17:20].strip()
                    self.chains[atom_number] = line[21:22].strip()
                    self.res_number[atom_number] = int(line[22:26].strip())
                    self.occupency[atom_number] = float(line[54:60].strip())
                    self.temperature_factors[atom_number] = float(line[60:66].strip())
                    self.element_symbol[atom_number] = line[76:78].strip()
                    self.atom_charge[atom_number] = line[78:80].strip()

        self.matroordinates = np.zeros((len(self.coordinates),3))
        # created for optimize distance calcul.
        i = 0
        for k,v in self.coordinates.items():
            self.matroordinates[i] = v
            i+=1

    def write_new_pdb(self, pdbpath = None):
        """
        Thanks to Pierre Poulain for the formated string line in PDB files:
        http://cupnet.net/pdb-format/.
        This function permits to write a new pdb file at the location of your choice.
        If no path is given, the old PDB will be overwrite.

        ARGUMENTS:

        pdbpath: the path to the new pdb.

        """
        if pdbpath is None:
        	pdbpath = self.pdb
        	
        f = open(pdbpath, "w")
        for k, v in self.coordinates.items():
            f.write("{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}\n"\
                .format(self.id[k], k, self.atom_type[k], " ",\
                self.res_type[k], self.chains[k], self.res_number[k], " ",\
                v[0], v[1], v[2], self.occupency[k],\
                self.temperature_factors[k], self.element_symbol[k], self.atom_charge[k]))
        f.close()



                

    def print_coordinates(self):
        """
        Simple function to print coordinates, no arguments.
        """
        for atom in self.coordinates:
            print("{} : {}".format(atom, self.coordinates[atom]))

    def translate(self, xtranslate=0, ytranslate=0, ztranslate = 0):
        """
        Function to translate x, y or z coordinates.

        ARGUMENTS:

        xtranslate: value of the translation on x axis.
        ytranslate: value of the translation on y axis.
        ztranslate: value of the translation on z axis.

        If no argument is given, default is 0.

        """
        translation_array = np.array([xtranslate, ytranslate, ztranslate])
        i = 0
        for k, v in self.coordinates.items():
            self.coordinates[k] += translation_array
            self.matroordinates[i] += translation_array
            i+=1

    def rotateX(self, angle):
        sin_t = np.sin(angle)
        cos_t = np.cos(angle)
        for point in self.matroordinates:
            y = point[1]
            z = point[2]
            point[1] = (y * cos_t - z * sin_t)*2
            point[2] = (z * cos_t + y * sin_t)*2

    def rotateY(self, angle):
        sin_t = np.sin(angle)
        cos_t = np.cos(angle)
        for point in self.matroordinates:
            x = point[0]
            z = point[2]
            point[0] = (x * cos_t + z * sin_t)*2
            point[2] = (z * cos_t - x * sin_t)*2

    def rotateZ(self, angle):
        sin_t = np.sin(angle)
        cos_t = np.cos(angle)
        for point in self.matroordinates:
            x = point[0]
            y = point[1]
            point[0] = (x * cos_t - y * sin_t)*2
            point[1] = (x * sin_t + y * sin_t)*2



    def get_min_distance_between_PDB(self, pdb):
        """
        Function which returns the minimum distance between two set of
        coordinates in PDBOperation object.

        ARGUMENTS:

        pdb: an other PDBOperation object. 
        """
        distance_matrix = distance.cdist(self.matroordinates, pdb.matroordinates)
        return np.amin(distance_matrix)