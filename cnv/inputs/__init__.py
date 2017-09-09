import os

def get_true_sight_inherited_disease_bed():
    return os.path.join(os.path.realpath(os.path.dirname(__file__)), 'TruSight_Inherited_Disease_Manifest_A.bed')

def get_true_sight_one_bed():
    return os.path.join(os.path.realpath(os.path.dirname(__file__)), 'TruSight-One-BED-May-2014.txt')

def get_dmd_exons():
    return os.path.join(os.path.realpath(os.path.dirname(__file__)), 'primary_DMD_exons.bed')