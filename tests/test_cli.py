#! /usr/bin/env python3

from click.testing import CliRunner
from cli import conv_test
from tempfile import NamedTemporaryFile

# Unit Test

def test_generate_files_in_cli():
    expected = '\n'.join(['ENCUT = 500     ! Plane wave cutoff',
                          'EDIFF = 1E-5   ! SCF convergence',
                          'ISTART = 1     ! Read existing wavefunction if exists',
                          'NSW = 0        ! No relaxation (single-point calculations for convergence testing',
                          'IBRION = -1    ! Determines how ions are moved (no ion movements with -1)',
                          'ISMEAR = 0    ! Type of smearing',
                          'SIGMA = 0.05     ! Smearing value in eV',
                          'ISYM = 2          ! Symmetry: 0-none; 2-GGA; 3-hybrid',
                          'PREC =  Accurate     ! Precision level',
                          'ALGO = Normal      ! Electronic minimisation algorithm',
                          'LASPH = .TRUE. ! Non-spherical elements',
                          'KPAR = 2      ! Number of kpoints that are to be parallelise',
                          'NCORE = 8     ! Number of cores working on a single orbital',
                          'Perdew-Burke-Ernzerhof for solids',
                          'GGA = PS               ! PBEsol'])
                          
    runner = CliRunner()

    with NamedTemporaryFile(mode='w+') as tf:
        result = runner.invoke(conv_test, ['generate_files', '500', '2 2 2', tf.name])
        assert result.exit_code == 0
        tf.seek(0)
        assert tf.read() == expected

test_generate_files_in_cli()
