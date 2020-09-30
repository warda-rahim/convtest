from click.testing import CliRunner
from cli import conv_test

# Unit Test

def test_generate_files_in_cli():
    runner = CliRunner()
    with runner.isolated_filesystems():
        with open('INCAR_test.txt', 'w') as f:
            f.write('ENCUT = 500     ! Plane wave cutoff\n 
                    EDIFF = 1E-5   ! SCF convergence\n 
                    ISTART = 1     ! Read existing wavefunction if exists\n
                    NSW = 0        ! No relaxation (single-point calculations for convergence testing\n 
                    IBRION = -1    ! Determines how ions are moved (no ion movements with -1)\n
                    ISMEAR = 0    ! Type of smearing\n
                    SIGMA = 0.05     ! Smearing value in eV\n
                    ISYM = 2          ! Symmetry: 0-none; 2-GGA; 3-hybrid\n
                    PREC =  Accurate     ! Precision level\n
                    ALGO = Normal      ! Electronic minimisation algorithm\n
                    LASPH = .TRUE. ! Non-spherical elements\n
                    KPAR = 2      ! Number of kpoints that are to be parallelise\n
                    NCORE = 8     ! Number of cores working on a single orbital\n
                    Perdew-Burke-Ernzerhof for solids\n
                      GGA = PS               ! PBEsol\n')


    result = runner.invoke(conv_test, ['500', '2 2 2', 'generate_files'])
    assert result.exit_code == 0
    self.assertTrue(filecmp.cmp('INCAR', 'INCAR_test.txt))

