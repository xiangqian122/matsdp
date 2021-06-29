#!/usr/bin/env python
from setuptools import setup, find_packages

def package_files(directory):
    import os
    paths = []
    for (path, directories, filenames) in os.walk(directory):
        for filename in filenames:
            paths.append(os.path.join('..', path, filename))
    return paths

extra_files = package_files('./matsdp/tests')


setup(name =  'matsdp',
      version = '0.2.4',
      description = 'The materials simulation and data processing toolkit',
      long_description = open('README.rst.txt').read(),
      author = 'Dianwu Wang',
      author_email  = 'dianwuwang@163.com',
      maintainer = 'Dianwu Wang',
      maintainer_email = 'dianwuwang@163.com',
      license = 'BSD',
      url = 'https://github.com/dianwdw/matsdp',
      packages = ['matsdp',
                  'matsdp.vasp',
                  'matsdp.wannier',
                  'matsdp.apt',
                  'matsdp.dvm',
		  'matsdp.pms',
                  ],
      package_dir = {'matsdp': './matsdp',
                     'matsdp.vasp' : './matsdp/vasp',
                     'matsdp.wannier': './matsdp/wannier',
                     'matsdp.apt' : './matsdp/apt',
                     'matsdp.dvm' : './matsdp/dvm',
		     'matsdp.pms' : './matsdp/pms',
                     },
      package_data = {'':extra_files,
                      },
      platforms = ['all'],
      install_requires=[
          'numpy',
          'scipy',
          'matplotlib',
          'scikit-learn'
          ],
      python_requires = '>=3.5',
      classifiers = [
          'Development Status :: 5 - Production/Stable',
          'Operating System :: OS Independent',
          'Intended Audience :: End Users/Desktop',
          'Intended Audience :: Science/Research',
          'License :: OSI Approved :: BSD License',
          'Programming Language :: Python',
          'Programming Language :: Python :: 3.7',
          ],
      keywords = [
          'VASP',
          'DOS',
	  'bands',
          'APT',
          'DVM',
          'data processing',
          ],          
      )
