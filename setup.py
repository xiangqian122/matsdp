#!/usr/bin/env python
from setuptools import setup, find_packages

setup(name =  'matsdp',
      version = '0.1.3',
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
                  'matsdp.apt'
                  ],
      package_dir = {'matsdp': './matsdp',
                     'matsdp.vasp' : './matsdp/vasp',
                     'matsdp.apt' : './matsdp/apt'
                     },
      platforms = ['all'],
      install_requires=[
          'numpy',
          'scipy',
          'matplotlib',
          'scikit-learn'
          ],
      classifiers = [
          'Development Status :: 5 - Production/Stable',
          'Operating System :: OS Independent',
          'Intended Audience :: End Users/Desktop',
          'Intended Audience :: Science/Research',
          'License :: OSI Approved :: BSD License',
          'Programming Language :: Python'
          ],
      keywords = [
          'VASP',
          'DOS',
          'APT',
          'data processing'
          ],          
      )
