#Last Updated: 21 April 2021

#Script for making database to run Jeanne scripts
#Create directory for each species, and subject
#Copy fine white, long, lat meshes from each subject in the BV database

#IMPORTANT: Change directory names and subject names according to species
#Command: /hpc/meca/users/auzias/brainvisa_4.5.0_patch_bugfix_2020_01_15/bin/python MakeDatabase.py 


import numpy as np
import os
import sys
from soma import aims

WORKDIR = "/hpc/meca/users/loh.k/interspecies_hiphop"
BVDIR = "/hpc/meca/data/Baboons/BV_Baboons"
CENTER = "Adrien"
SPECIES = "baboon"

SubNameList = [i for i in os.listdir(BVDIR + "/" + CENTER) if '.' not in i]
#print(len(SubNameList))

# Make species output directory

os.system('mkdir ' + (WORKDIR + "/" + SPECIES))

for s in range(len(SubNameList)):

    SubName = (SubNameList[s].split("_")[4]) #baboons
    #SubName = (SubNameList[s]) #macaques

    # Filenames to copy
    LwhiteMesh = (BVDIR + "/" + CENTER + "/" + SubNameList[s] 
            + "/t1mri/default_acquisition/default_analysis/segmentation/mesh/"
            + SubNameList[s] + "_Lwhite_fine.gii") 
    RwhiteMesh = (BVDIR + "/" + CENTER + "/" + SubNameList[s] 
            + "/t1mri/default_acquisition/default_analysis/segmentation/mesh/"
            + SubNameList[s] + "_Rwhite_fine.gii")
    LwhiteLon = (BVDIR + "/" + CENTER + "/" + SubNameList[s] 
            + "/t1mri/default_acquisition/default_analysis/segmentation/mesh/"
            + "surface_analysis/" + SubNameList[s] + "_Lwhite_lon.gii")
    RwhiteLon = (BVDIR + "/" + CENTER + "/" + SubNameList[s] 
            + "/t1mri/default_acquisition/default_analysis/segmentation/mesh/"
            + "surface_analysis/" + SubNameList[s] + "_Rwhite_lon.gii")
    LwhiteLat = (BVDIR + "/" + CENTER + "/" + SubNameList[s] 
            + "/t1mri/default_acquisition/default_analysis/segmentation/mesh/"
            + "surface_analysis/" + SubNameList[s] + "_Lwhite_lat.gii")
    RwhiteLat = (BVDIR + "/" + CENTER + "/" + SubNameList[s] 
            + "/t1mri/default_acquisition/default_analysis/segmentation/mesh/"
            + "surface_analysis/" + SubNameList[s] + "_Rwhite_lat.gii")

    SubjOutDir = WORKDIR + "/" + SPECIES + "/" + SubName + "_" + SPECIES

    os.system("mkdir " + SubjOutDir)
    os.system("cp " + LwhiteMesh + " " + SubjOutDir + "/" + SubName + "_Lwhite_fine.gii")
    os.system("cp " + RwhiteMesh + " " + SubjOutDir + "/" + SubName + "_Rwhite_fine.gii")
    os.system("cp " + LwhiteLon + " " + SubjOutDir + "/" + SubName + "_Lwhite_lon.gii")
    os.system("cp " + RwhiteLon + " " + SubjOutDir + "/" + SubName + "_Rwhite_lon.gii")
    os.system("cp " + LwhiteLat + " " + SubjOutDir + "/" + SubName + "_Lwhite_lat.gii")
    os.system("cp " + RwhiteLat + " " + SubjOutDir + "/" + SubName + "_Rwhite_lat.gii")

    print("Files copied for: " + SubNameList[s])


