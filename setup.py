from setuptools import setup
import sys

req = ['biopython==1.70', 'Cython==0.24', 'pysam==0.13', 'rpy2>=2.7.0,<=2.8.6']

if int(sys.version.split()[0][0]) > 2:
    req = ['biopython>=1.70', 'Cython>=0.24', "pysam@https://github.com/pysam-developers/pysam/tarball/master#egg=pysam", 'rpy2>=2.7.0']

setup(name='Xome_Blender',
      version='0.1',
      description='Generates synthetic cancer genomes with different contamination level and intra-tumor heterogeneity and devoid of any synthetic element',
      url='https://github.com/rsemeraro/XomeBlender',
      author='Roberto Semeraro',
      author_email='robe.semeraro@gmail.com',
      license='LICENSE.txt',
      package_data={'xomeblender': ['*py']},
      packages=['xomeblender'],
      include_package_data=False,
      install_requires=req
)
