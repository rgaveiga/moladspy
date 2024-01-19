import os.path
from setuptools import setup

def read(fname):
    return open(os.path.join(os.path.dirname(__file__),fname)).read()

setup(name='MolAdsPy',
      version='1.0.0',
      description="Constructing complex atomic scale structures with Python.",
      long_description_content_type='text/markdown',
      long_description=read('Desc.MD'),
      install_requires=['multipledispatch'],
      author='Roberto Gomes de Aguiar Veiga',
      url="https://github.com/rgaveiga/moladspy",
      packages=['MolAdsPy'])


