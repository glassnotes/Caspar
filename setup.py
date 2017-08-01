from setuptools import setup

setup(name='caspar',
      version='1.0.0',
      description='Factorization of SU(n) transformations in Python.',
      url='https://github.com/glassnotes/Caspar',
      author='Olivia Di Matteo',
      author_email='odimatte@uwaterloo.ca',
      license='BSD',
      package_dir={'': 'src'},
      packages=['caspar'],
      zip_safe=False)
