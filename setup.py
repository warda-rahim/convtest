from setuptools import setup, find_packages

setup(
    name='ConvTest',
    version='1.0.0',
    description='Cutoff and kmesh convergence testing for ab initio '
                 'solid-state calculations using VASP',
    long_description=open('README.md').read(),
    author='Warda Rahim',
    author_email='wardarahim25@gmail.com',
    url="https://github.com/warda-rahim/convergence-test",
    license='MIT',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.5',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Topic :: Scientific/Engineering :: Physics'
        ],
    keywords='chemistry dft vasp pymatgen',
    packages=find_packages(),
    install_requires=['numpy', 'matplotlib'],
    extras_require={':python_version=="3.5"': [
                    'pymatgen >=2016.12.30, <=2019.6.20']
                    },
    entry_points='''
        [console_scripts]
        conv_test=convtest.cli:conv_test
        ''',
)
