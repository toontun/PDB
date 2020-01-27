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

    """
    Class which contains PDB data.
    """

    def __init__(self, pdbpath):
        """
        Init an PDBOperation object by parsing a PDB file. Just give the path
        to the PDB file of your choice. Coordinates are also stored in a numpy
        matrix for optimization while distance is calculated. Other data are
        stored in dictionaries to keep the right information on each atoms for
        further manipulations.

        ARGUMENTS:

        pdbpath: the path to the pdb.

        """

        # pylint: increase [design] max-attributes to 12.
        # 12 is reasonable in this case.

        self.pdb_file = pdbpath
        self.line_start = {}
        self.coordinates = {}
        self.atom_type = {}
        self.res_type = {}
        self.chains = {}
        self.res_number = {}
        self.occupency = {}
        self.temperature_factors = {}
        self.element_symbol = {}
        self.atom_charge = {}
        with open(pdbpath, "r") as opened_pdb:
            for line in opened_pdb:
                line_start = line[0:6].strip()
                if line_start in ("ATOM", "HETATM"):
                    atom_number = int(line[6:11].strip())
                    self.line_start[atom_number] = line_start
                    self.coordinates[atom_number] =\
                        np.array([float(line[30:38].strip()),
                                  float(line[38:46].strip()),
                                  float(line[46:54].strip())])
                    self.atom_type[atom_number] = line[12:16].strip()
                    self.res_type[atom_number] = line[17:20].strip()
                    self.chains[atom_number] = line[21:22].strip()
                    self.res_number[atom_number] = int(line[22:26].strip())
                    self.occupency[atom_number] = float(line[54:60].strip())
                    self.temperature_factors[atom_number] =\
                        float(line[60:66].strip())
                    self.element_symbol[atom_number] = line[76:78].strip()
                    self.atom_charge[atom_number] = line[78:80].strip()

        self.matroordinates = np.zeros((len(self.coordinates), 3))
        # created to optimize distance calcul.
        i = 0
        for key in self.coordinates:
            self.matroordinates[i] = self.coordinates[key]
            i += 1

    def write_new_pdb(self, pdbpath=None):
        """
        Thanks to Pierre Poulain for the formated string line in PDB files:
        http://cupnet.net/pdb-format/.
        This function permits to write a new pdb file at the location of your
        choice. If no path is given, the old PDB will be overwrite.

        ARGUMENTS:

        pdbpath: the path to the new pdb.

        """
        if pdbpath is None:
            pdbpath = self.pdb_file

        opened_file = open(pdbpath, "w")
        for k, value in self.coordinates.items():
            opened_file.write("{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}\
                {:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}\
                {:6.2f}          {:>2s}{:2s}\n"\
                .format(self.line_start[k], k, self.atom_type[k], " ",\
                        self.res_type[k], self.chains[k], self.res_number[k],
                        " ", value[0], value[1], value[2], self.occupency[k],\
                        self.temperature_factors[k], self.element_symbol[k],\
                        self.atom_charge[k]))
        opened_file.close()

    def print_coordinates(self):
        """
        Simple function to print coordinates, no arguments.
        """
        for atom in self.coordinates:
            print("{} : {}".format(atom, self.coordinates[atom]))

    def translate(self, xtranslate=0, ytranslate=0, ztranslate=0):
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
        for key in self.coordinates:
            self.coordinates[key] += translation_array
            self.matroordinates[i] += translation_array
            i += 1

    def get_distpdb(self, pdb):
        """
        Function which returns the minimum distance between two set of
        coordinates in PDBOperation object.

        ARGUMENTS:

        pdb: an other PDBOperation object.
        """
        distance_matrix = distance.cdist(self.matroordinates,
                                         pdb.matroordinates)
        return np.amin(distance_matrix)
