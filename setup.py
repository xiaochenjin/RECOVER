from Cython.Build import cythonize
from setuptools import setup
#from distutils.core import setup
from distutils.extension import Extension
import os

defs = os.path.join('src')

long_description = open('./README.md', 'r').read()

#compile cython codes
ext_modules = [
	Extension("find_dK", ["src/find_dK.pyx"]),
    Extension("PairCorrelation_AA", ["src/PairCorrelation_AA.pyx"]),
	Extension("PairCorrelation_AB", ["src/PairCorrelation_AB.pyx"]),
	Extension("KNN_info_AA", ["src/KNN_info_AA.pyx"]),
	Extension("KNN_info_AB", ["src/KNN_info_AB.pyx"]),
	Extension("KNN_track_AA", ["src/KNN_track_AA.pyx"]),
	Extension("KNN_track_AA", ["src/KNN_track_AA.pyx"])
]

setup(name='RECOVER',
      version='1.0',
      description='Retrieve SRO parameters from perturbed APT data',
      long_description=long_description,
      long_description_content_type='text/markdown',
      author='Xiaochen Jin',
      author_email='xcjin@gwmail.gwu.edu',
      platforms='Unix',
      url='https://github.com/xiaochenjin/RECOVER',
      packages=['RECOVER'],
      package_dir={'RECOVER':'src'},
      install_requires=['numpy','scipy>=1.15.2','Cython'],
      extra_requires=['matplotlib'],
      python_requires='>=3.10',
	  ext_modules = ext_modules
)
