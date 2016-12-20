from __future__ import division, absolute_import, print_function

import os
import sys
from setuptools import setup
from numpy.distutils.core import setup
from numpy.distutils.misc_util import Configuration
from setuptools.command.test import test as TestCommand


class PyTest(TestCommand):
    def finalize_options(self):
        TestCommand.finalize_options(self)
        self.verbose = True

    def run_tests(self):
        import pytest
        errno = pytest.main(self.test_args)
        sys.exit(errno)


def extract_version(module='mocsy'):
    version = None
    fdir = os.path.dirname(__file__)
    fnme = os.path.join(fdir, module, '__init__.py')
    with open(fnme) as fd:
        for line in fd:
            if (line.startswith('__version__')):
                _, version = line.split('=')
                # Remove quotation characters.
                version = version.strip()[1:-1]
                break
    return version

rootpath = os.path.abspath(os.path.dirname(__file__))


def read(*parts):
    return open(os.path.join(rootpath, *parts), 'r').read()


long_description = '{}\n{}'.format(read('README.rst'), read('ChangeLog'))
LICENSE = read('LICENSE')

with open('requirements.txt') as f:
    require = f.readlines()
install_requires = [r.strip() for r in require]


sources = [
    'src/singledouble.f90',
    'src/p80.f90',
    'src/sw_adtg.f90',
    'src/sw_ptmp.f90',
    'src/sw_temp.f90',
    'src/constants.f90',
    'src/p2fCO2.f90',
    'src/phsolvers.f90',
    'src/rho.f90',
    'src/rhoinsitu.f90',
    'src/tis.f90',
    'src/tpot.f90',
    'src/varsolver.f90',
    'src/vars.f90',
    'src/depth2press.f90',
    'src/f2pCO2.f90',
    'src/buffesm.f90',
    'src/DNAD.f90',
    'src/derivauto.f90',
    'src/derivnum.f90',
    'src/errors.f90',
    'src/gasx.f90'
    ]

config = Configuration('')

config.add_extension('_mocsy',
                     sources=sources)

setup(name='mocsy',
      version=extract_version(),
      license=LICENSE,
      long_description=long_description,
      classifiers=['Development Status :: 5 - Production/Stable',
                   'Environment :: Console',
                   'Intended Audience :: Science/Research',
                   'Intended Audience :: Developers',
                   'Intended Audience :: Education',
                   'License :: OSI Approved :: MIT License',
                   'Operating System :: OS Independent',
                   'Programming Language :: Python',
                   'Topic :: Scientific/Engineering',
                   'Topic :: Education',
                   ],
      description='Routines to model ocean carbonate system thermodynamics',
      url='https://github.com/jamesorr/mocsy.git',
      platforms='any',
      keywords=['carbonate', 'ocean acidification', 'CO2SYS'],
      install_requires=install_requires,
      packages=['mocsy'],
      ext_package='mocsy',
      tests_require=['pytest'],
      cmdclass=dict(test=PyTest),
      **config.todict()
      )
