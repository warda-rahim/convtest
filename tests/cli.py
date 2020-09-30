#! /usr/bin/env python3

import click
import os
import numpy as np
from pathlib import Path

from convtest.write_files import prepare_incar, write_kpoints
from convtest.conv_setup import make_cutoff_folders, make_kmesh_folders, copy_files, update_values
from convtest.conv_plot import extract_cutoff_energies, extract_kmesh_energies, create_dataframes, print_converged_values, plot_conv


@click.group('conv_test')
def conv_test():
    """ Prepare convergence testing and analyse data """
    pass

@conv_test.command('generate_input_files')
@click.argument('encut', type=int)
@click.argument('kmesh', nargs=3)
@click.option('-f', '--functional', default='ps', help='functional for convergence test')
@click.option('--aexx', type=float, default=0.25, help='amount of HF-exchange for hybrid functional (AEXX)')
@click.option('--systype', default='semiconductor', help='type of system metal or semiconductor')
@click.option('--ismear', type=int, default='0', help='type of smearing')
@click.option('--sigma', type=float, default=0.05, help='smearing value in eV')
@click.option('--sym', help='symmetry scheme')
@click.option('--algo', default='Normal', help='electronic minimization algorithm')
@click.option('--prec', default='Accurate', help='precision-mode')
@click.option('--computer', help='cluster where calculations will be submitted')
@click.option('--ncore', type=int,  help='number of core working per band')
@click.option('--kpar', type=int, help='kpoint parallelism')
def generate_input_files(kmesh, encut, functional, aexx, systype, ismear, sigma, sym, algo, prec, computer, ncore, kpar):
    """ Generates input files for convergence tests 
        
        (encut is the cutoff used for all kmesh conv_tests)
        (kmesh is the kpoints mesh used for all cutoff conv_tests)
    """

    prepare_incar(encut, functional, aexx, systype, ismear, sigma, sym, algo, prec, computer, ncore, kpar)
    write_kpoints(kmesh)



path_cutoff = os.path.join(os.getcwd(), 'cutoff_conv')
path_kmesh = os.path.join(os.getcwd(), 'kmesh_conv')


@conv_test.command('generate_folders')
@click.option('--emin', type=int, default=200, help='minimum encut for cutoff conv_test')
@click.option('--emax', type=int, default=700, help='maximum encut for cutoff conv_test')
@click.option('--estep', type=int, default=50, help='step for encut for cutoff conv_test')
@click.option('-k', '--kinputs', default=['2 2 2', '4 4 4', '6 6 6', '8 8 8'], multiple=True, help='Kmesh values for conv_tests')
def generate_folders(emin, emax, estep, kinputs):
    """ Creates folders for cutoff and kmesh convergence tests """

    # Creates folders for cutoff convergence tests
    if not os.path.isdir(path_cutoff):      
        os.makedirs(path_cutoff)
    einputs = list(map(str, np.arange(emin, emax+estep, estep)))
    make_cutoff_folders(path_cutoff, einputs)

   
    # Creates folders for cutoff convergence tests
    if not os.path.isdir(path_kmesh):
        os.makedirs(path_kmesh)
    make_kmesh_folders(path_kmesh, kinputs)

    # Copies files to each folder
    copy_files(path_cutoff, einputs, path_kmesh, kinputs)
    

    p_cutoff = Path(path_cutoff)
    # Lists subdirectories in p_cutoff
    dir_cutoff = [x for x in p_cutoff.iterdir() if x.is_dir()]

    for subfolder in dir_cutoff:
        incar = subfolder / 'INCAR'
        # Update the INCAR file if it exists in the folder
        if not incar.is_file():
            continue
        else:
            update_values(incar, subfolder / 'INCAR', subfolder.name)


    p_kmesh = Path(path_kmesh)
    # Lists subdirectories in p_kmesh
    dir_kmesh  = [x for x in p_kmesh.iterdir() if x.is_dir()]

    for subfolder in dir_kmesh:
        kpoints = subfolder / 'KPOINTS'
        # Update the KPOINTS file if it exists in the folder
        if not kpoints.is_file():
            continue
        else:
            update_values(kpoints, subfolder/ 'KPOINTS', subfolder.name)



@conv_test.command('analyse') 
def analyse_conv():
    """ Generates plots and csvfiles for output convergence test data """

    # Extract energies
    cutoff_energy_dicts = extract_cutoff_energies(path_cutoff)
    kmesh_energy_dicts = extract_kmesh_energies(path_kmesh)

    # Create dataframes
    df_cutoff, df_kmesh = create_dataframes(cutoff_energy_dicts, kmesh_energy_dicts, path_kmesh)

    # Print converged values
    cutoff_conv, kmesh_conv = print_converged_values(df_cutoff, df_kmesh)

    # Export data to a csv file
    cutoff_energy_file = 'cutoff-energy.csv'
    kmesh_energy_file = 'kmesh-energy.csv'

    df_cutoff.to_csv(cutoff_energy_file, index=False)
    df_kmesh.to_csv(kmesh_energy_file, index=False)

    # Plot convergence data
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec
    from cycler import cycler
    import matplotlib.ticker as ticker
    
    #plt.style.use("pretty2")
    
    mpl.rcParams['axes.linewidth'] = 2
    
    fig = plt.figure(figsize=(30, 14))
    
    gs = GridSpec(1, 2, wspace=0.5)
    ax = ['', '']
    
    ax[0] = fig.add_subplot(gs[:,:-1])
    ax[1] = fig.add_subplot(gs[:,-1])
    
    fig.align_xlabels()
    
    for i in range(2):
        ax[i].tick_params(axis='both', which='major', direction='in', length=20, labelsize=42, pad=15)
        ax[i].tick_params(axis='both', which='minor', direction='in', length=8, pad=15)
        ax[i].yaxis.set_major_locator(ticker.MaxNLocator(6))
        ax[i].yaxis.set_minor_locator(ticker.AutoMinorLocator(2))
    
    ax = plot_conv(ax, cutoff_energy_file, kmesh_energy_file, cutoff_conv, kmesh_conv)
    
    plt.subplots_adjust(left=0.12, right=0.98, bottom=0.25, top=0.95)
    plt.savefig('encut-kpt-converg-test.pdf')
    plt.savefig('encut-kpt-converg-test.png')

if __name__ == '__main__':
    conv_test()
