from setuptools import setup

setup(name='Xome_Blender',
      version='0.1',
      description='Generates synthetic cancer genomes with different contamination level and intra-tumor heterogeneity and devoid of any synthetic element',
      url='https://github.com/rsemeraro/XomeBlender',
      author='Roberto Semeraro',
      author_email='robe.semeraro@gmail.com',
      license='LICENSE.txt',
      scripts=['xome_blender', 'inxalizer'],
      packages=['xomeblender', 'xomeblender.scripts'],
      install_requires=['biopython==1.70', 'Cython==0.24', 'pysam==0.13', 'rpy2==2.7.0']
)
