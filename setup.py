from setuptools import setup, find_packages
import os
from cnv import __version__

VERSION = __version__


# Hack to avoid specifying requirements in two places, see:
# http://stackoverflow.com/questions/14399534/how-can-i-reference-requirements-txt-for-the-install-requires-kwarg-in-setuptool
requirements_file = os.path.join(os.path.dirname(__file__), "requirements.txt")
with open(requirements_file) as f:
    all_reqs = f.read().splitlines()

setup(name='cnv',
      version=VERSION,
      url='https://github.com/genepeeks/dmd',
      description='Library and command line scripts to find copy number variation.',
      packages=find_packages(),
      install_requires=required,
      include_package_data=True,
      package_data={'': ['inputs/*.txt', 'inputs/*.bed']},
      dependency_links=[],
      entry_points={
          'console_scripts': ['cnv=cnv.cli:main']
      },
      zip_safe=False)
