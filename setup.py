import os
from setuptools import setup
from numpy.distutils.core import Extension, setup


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


long_description = '{}\n{}'.format(read('README.rst'), read('CHANGES.txt'))
LICENSE = read('LICENSE.txt')

with open('requirements.txt') as f:
    require = f.readlines()
install_requires = [r.strip() for r in require]


setup(name="mocsy",
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
      ext_modules=[Extension(name='_mocsy',
                             sources=['mocsy/_mocsy.f90'])],
      maintainer='Filipe Fernandes',
      maintainer_email='ocefpaf@gmail.com',
      )
