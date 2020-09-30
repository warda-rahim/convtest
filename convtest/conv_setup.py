#! /usr/bin/env python3

import click
import os
import shutil
import fileinput
import subprocess
import re
import numpy as np
from natsort import natsorted
import pandas as pd
from pandas import DataFrame
import csv
from pathlib import Path


def make_cutoff_folders(path_cutoff, einputs):
    """ Make cutoff folders 

    Args:
        path_cutoff (:obj:'str'): Path to folder inside which subfolders for cutoff convergence will be created.
        einputs (:obj:'list'): List of cutoff/encut values for convergence test.

    """
    for i in einputs:
        path = os.path.join(path_cutoff, i)
        if not os.path.exists(path):
           os.makedirs(path)



def make_kmesh_folders(path_kmesh, kinputs):
    """ Make kmesh folders

    Args:
        path_kmesh (:obj:'str'): Path to folder insider which subfolders for kmesh convergence will be created.
        kinputs (:obj:'list'):  List of kmesh values for convergence test.

    """

    sum_k = sum(len(i.replace(" ", "")) for i in kinputs)
    if sum_k % 3 == 0:
        pass
    else:
        raise Exception('Invalid number of kmeshes\n'
                           'Each element in the array should contain 3 integers')
    
    for i in kinputs:
        path = os.path.join(path_kmesh, i.replace(" ", "x"))
        if not os.path.exists(path):
           os.makedirs(path)
    


def copy_files(path_cutoff, einputs, path_kmesh, kinputs):
    """ Copy input files inside each folder for cutoff and kmesh convergence test.

    Args:
        path_cutoff (:obj:'str'): Path to folder inside which subfolders for cutoff convergence are present.
        einputs (:obj:'list'): List of cutoff/encut values for convergence test.
        path_kmesh (:obj:'str'): Path to folder insider which subfolders for kmesh convergence are present.
        kinputs (:obj:'list'):  List of kmesh values for convergence test.
    
    """
   
    destinations  = [os.path.join(path_cutoff, i) for i in einputs]
    destinations.extend(os.path.join(path_kmesh, i.replace(" ","x")) for i in kinputs)

    files = ['POSCAR', 'POTCAR', 'KPOINTS', 'INCAR', 'job']

    for dest in destinations:
        for f in files:
            shutil.copy(f, dest)
    


def update_values(filepath, save_path, value):
    """ Create a new file with updated values.
    
    Args:
        filepath (:obj:'str'): Path to the input file.
        save_path (:obj:'str'): Path to the output/updated file.
        value (:obj:'str'): New/updated cutoffs and kmeshes.

        """

    with open(filepath, 'r') as fin:               # Read the input file
        lines = fin.readlines()

    if filepath.name == 'INCAR':
        with open(save_path, 'w') as fout:          # Write the output INCAR file with updated value
            for line in lines:
                if line.strip().lower().startswith('encut'):
                    line = 'encut =' + value + '\n'
                fout.write(line)

    elif filepath.name == 'KPOINTS':               # Write the output KPOINTS file with updated value
        with open(save_path, 'w') as fout:
            lines[-1] = ' '.join([i for i in value.split('x')])
            for line in lines:
                fout.write(line)




