#! /usr/bin/env python3

from pymatgen.io.vasp import Vasprun, Outcar
import os
import numpy as np
from natsort import natsorted
import pandas as pd
from pandas import DataFrame
import csv


def extract_cutoff_energies(path_cutoff):
    """ Extracts energy values from each cutoff/encut convergence folder.

    Args:
        path_cutoff: (:obj:'str'): Path to folder containing folders for cutoff convergence
     
    Returns:
        :obj:'list': List of dictionaries containing encuts/cutoffs and corresponding energy values.

    """
    
    data_dicts_encut = []

    for root, dirs, files in os.walk(path_cutoff):	
        for name in files:
            if name == 'vasprun.xml':  
                path = os.path.join(root,name)
                vasprun = Vasprun(path)
                energy = vasprun.final_energy
                encut = os.path.basename(os.path.dirname(path))
                data_dicts_encut.append({'Encut': encut, 'Energy': energy})
                
    data_dicts_encut = sorted(data_dicts_encut, key=lambda k: k['Encut'])

    return data_dicts_encut


def extract_kmesh_energies(path_kmesh):
    """ Extracts energy values from each kmesh convergence folder.
    
    Args:
        path_kmesh: (:obj:'str'): Path to folder containing folders for kmesh convergence.
    
    Returns:
        :obj:'list': List of dictionaries containing kmeshes and corresponding energy values.

    """

    data_dicts_kpt = []

    for root, dirs, files in os.walk(path_kmesh):
        for name in files:
            if name == 'vasprun.xml':
                path = os.path.join(root,name)
                vasprun = Vasprun(path)
                energy = vasprun.final_energy
                kpt = os.path.basename(os.path.dirname(path))
                data_dicts_kpt.append({'Mesh': kpt, 'Energy': energy})
    
    data_dicts_kpt = natsorted(data_dicts_kpt, key=lambda k: k['Mesh'])
    
    return data_dicts_kpt



def create_dataframes(cutoff_dicts, kmesh_dicts, path_kmesh):
    """ Creates dataframes from lists containing dictionaries with cutoffs, kmeshes and energies.

    Args:
        cutoff_dicts: (:obj:'list') List of dictionaries with cutoffs and energies.
        kmesh_dicts: (:obj:'list') List of dictionaries with kmeshes and energies.
        path_kmesh: (:obj:'str') Path to folder containing folders for kmesh convergence.

    Returns:
        obj:'dtype:object': Two dataframes with cutoffs and kmeshes as columns respectively,
                        along with columns with energy values.

    """

    df1 = DataFrame(cutoff_dicts)
    df2 = DataFrame(kmesh_dicts)
    df2 = df2[['Mesh', 'Energy']]
    
    for root, dirs, files in os.walk(path_kmesh):       # finds number of atoms in the structure using vasprun.xml file
        for name in files:
            if name == 'vasprun.xml':
               path = os.path.join(root,name)
               vasprun = Vasprun(path)
               num_atoms = vasprun.structures[-1].composition.num_atoms
    
    df1['E per atom/meV'] = df1['Energy'] / num_atoms * 1000
    df2['E per atom/meV'] = df2['Energy'] / num_atoms * 1000
    
    df1['E diff per atom/meV'] = df1['E per atom/meV'].diff(periods=1)
    df2['E diff per atom/meV'] = df2['E per atom/meV'].diff(periods=1)

    return df1, df2



def print_converged_values(df_cutoff, df_kmesh):

    """ Prints the converged value of cutoff/encut and kmesh.

    Args:
        df_cutoff: (:obj:'dataframe') Dataframe with cutoffs and energy values.
        df_kmesh: (:obj:'dataframe') Dataframe with kmeshes and energy values.

    Returns:
        obj:'str': Converged cutoff and kmesh values.

    """

    index1 = df_cutoff[df_cutoff['E diff per atom/meV'].abs() <= 3].index[0]
    index2 = df_kmesh[df_kmesh['E diff per atom/meV'].abs() <= 1].index[0]
    
    converg_cutoff = df_cutoff.iloc[index1]['Encut']
    converg_kmesh = df_kmesh.iloc[index2]['Mesh']
    
    print("Converged encut is: {}".format(converg_cutoff) + '\n' + "Converged kmesh is: {}".format(converg_kmesh))
    return converg_cutoff, converg_kmesh



def plot_conv(axis, cutoff_energy_file, kmesh_energy_file, converg_encut, converg_kmesh):
    """ Generates convergence plots.

    Args:
        axis: (:obj:'matplotlib.pyplot'): Axis to plot on.
        cutoff_energy_file: (:obj:'file') csv file with cutoffs and energy values.
        kmesh_energy_file: (:obj:'file') csv file with kmeshes and energy values.

     Returns:
         :obj:'matplotlib.pyplot': Axis with convergence plots.

    """

    with open(cutoff_energy_file) as csvfile:
        reader = csv.reader(csvfile)
    
        encut, energy1 = [], []
    
        next(reader)
        
        for row in reader:
            encut.append(row[0])
            energy1.append(round(float(row[2]),1))
    
    with open(kmesh_energy_file) as csvfile:
        reader = csv.reader(csvfile)
    
        kmesh, energy2 = [], []
    
        next(reader)
    
        for row in reader:
            kmesh.append(row[0]) 
            energy2.append(round(float(row[2]),1))
    
    index_encut = encut.index(converg_encut)
    index_kpt = kmesh.index(converg_kmesh)

    colours_encut = ["#CC3366" for i in range(len(encut))]
    colours_encut[index_encut] = "#9FC131"
    
    colours_kpt = ["#CC3366" for i in range(len(kmesh))]
    colours_kpt[index_kpt] = "#9FC131"
    
    axis[0].set_xlabel('Plane wave cutoff (eV)', fontsize=50)
    axis[0].set_ylabel(r'$\mathregular{Energy\ (meV\ atom^{-1}}$)', fontsize=50)
    axis[1].set_xlabel('k-point mesh', fontsize=50)
    axis[1].set_ylabel(r'$\mathregular{Energy\ (meV\ atom^{-1}}$)', fontsize=50)

    axis[0].set_xticks(np.arange(len(encut)))
    axis[0].set_xticklabels(encut, rotation=45, ha="right")
    axis[1].set_xticks(np.arange(len(kmesh)))
    axis[1].set_xticklabels(kmesh, rotation=45, ha="right")
 
    axis[0].plot(encut, energy1, linewidth=2.5, color="#000000", marker='h', markersize=20, markerfacecolor="#CC3366", markeredgecolor="#000000")
    for i in range(len(encut)):
        axis[0].plot(encut[i], energy1[i], color="#000000", marker='h', markersize=20, markerfacecolor=colours_encut[i], markeredgecolor="#000000")  
    
    axis[1].plot(kmesh, energy2, linewidth=2.5, color="#000000", marker='h', markersize=20, markerfacecolor="#CC3366", markeredgecolor="#000000")
    for i in range(len(kmesh)):
        axis[1].plot(kmesh[i], energy2[i], color="#000000", marker='h', markersize=20, markerfacecolor=colours_kpt[i], markeredgecolor="#000000")
   
    return axis




