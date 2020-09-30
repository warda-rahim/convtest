#! /usr/bin/env python3

import click
import os
import warnings
import numpy as np
from natsort import natsorted
import pandas as pd
from pandas import DataFrame
import csv


def prepare_incar(encut, functional, aexx, systype, ismear, sigma, sym, algo, prec, computer, ncore, kpar):
    """ Write INCAR file 

    (ENCUT is the cutoff value)
    """

    defncore = {'archer': 12,
                'michael': 10,
                'thomas': 12,
                'kathleen': 10,
                'grace': 8,
                'default': 8}
    
    if computer is None:
        try:
            computer = os.environ['computer'].lower()
        except Exception:
            pass
    else: 
        computer = computer.lower()
    
    if ncore is None:
        try:
            ncore = defncore[computer]
        except Exception:
            ncore = defncore['default']
    else:
        ncore = ncore
    
    
    gga = ['lda', 'pbe', 'pbesol', 'ps', 'scan', 'tb', 'mbj']
    hf = ['hf', 'hse', 'pbe0']
    
    if functional.lower() in gga:
        isym = '2'
        algo = 'Normal'
    elif functional.lower() in hf:
        isym = '3'
        algo = 'All'
  

    if sym is not None:
       isym = sym
    if algo is not None:
       algo = algo
   
    
    if systype == 'metal':
        ismear = '1'
        sigma = '0.2'
    else:
        ismear = '0'
        sigma = '0.05'
      
    if ismear is not None:
        ismear = ismear
    if sigma is not None:
        sigma = sigma
      
     
    INCAR = ['ENCUT = {}     ! Plane wave cutoff'.format(encut),
             'EDIFF = 1E-5   ! SCF convergence', 
             'ISTART = 1     ! Read existing wavefunction if exists',
             'NSW = 0        ! No relaxation (single-point calculations for convergence testing',
             'IBRION = -1    ! Determines how ions are moved (no ion movements with -1)',
             'ISMEAR = {}    ! Type of smearing'.format(ismear),
             'SIGMA = {}     ! Smearing value in eV'.format(sigma),
             'ISYM = {:<5}      ! Symmetry: 0-none; 2-GGA; 3-hybrid'.format(isym),
             'PREC =  {}     ! Precision level'.format(prec),
             'ALGO = {}      ! Electronic minimisation algorithm'.format(algo),
             'LASPH = .TRUE. ! Non-spherical elements',
             'KPAR = {}      ! Number of kpoints that are to be parallelised'.format(kpar),
             'NCORE = {}     ! Number of cores working on a single orbital'.format(ncore)]
      
    functionals = {
    'hf':      ['\nHartree-Fock',
                '  LHFCALC = .TRUE.       ! Activate HF',
                '  AEXX = 1               ! 100 % HF',
                '  ALDAC = 0              ! No LDA correllation',
                '  AGGAC = 0              ! No correllation correction\n'],
      
    'hse':     ['\nHeyd-Scuseria-Ernzerhof',
                '  GGA = PE               ! PBE-based',
                '  LHFCALC = .TRUE.       ! Activate HF',
                '  PRECFOCK = FAST        ! Controls the FFT grids used in exact exchange routines',
                '  AEXX = {:<8.2f}        ! Proportion HF exchange'.format(aexx),
                '  HFSCREEN = 0.207       ! Screened exchange\n'],
      
    'lda':     [],
      
    'pbe':     ['\nPerdew-Burke-Ernzerhof',
                '  GGA = PE               ! PBE-based\n'],
      
    'pbe0':    ['\nPerdew-Burke-Ernzerhof 0',
                '  GGA = PE               ! PBE based',
                '  LHFCALC = .TRUE.       ! Activate HF',
                '  PRECKFOCK = FAST       ! HF FFT grid',
                '  AEXX = {:<8.2f}        ! Proportion HF exchange\n'.format(aexx)],
      
    'pbesol':  ['\nPerdew-Burke-Ernzerhof for solids',
    '  GGA = PS               ! PBEsol\n'],
      
    'scan':    ['\nStrongly constrained and appropriately normed semilocal density functional',
                '  METAGGA = SCAN         ! SCAN',
                '  LMIXTAU = .TRUE.       ! Kinetic energy mixing\n'],
      
    'tb':      ['\nTran-Blaha (modified Becke-Johnson)',
                '  METAGGA = MBJ          ! TB/mBJ\n']}
      
    functionals['ps'] = functionals['pbesol']
    functionals['mbj'] = functionals['tb']

    if functional.lower() == 'lda':
        warnings.warn('Remember to use LDA POTCARS')

   
    with open('INCAR', 'w') as f:
        f.write('\n'.join(INCAR))
        if functional:
            f.write('\n'.join(functionals[functional.lower()])) 



def write_kpoints(kmesh):
    """ Write KPOINTS file 

    (Kmesh is the kpoints mesh)
    """
    
    with open('KPOINTS', 'w') as f:
        f.write("Automatic mesh \n0 \nGamma \n{} {} {}".format(kmesh[0], kmesh[1], kmesh[2]))

