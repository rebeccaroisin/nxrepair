from distutils.core import setup
setup(
  name = 'nxrepair',
  packages = ['nxrepair'], # this must be the same as the name above
  entry_points = {
              'console_scripts': [
                  'nxrepair = nxrepair.nxrepair:main',                  
              ],              
          },
  version = '0.13',
  description = 'Automated error detection using Nextera Mate Pair data',
  author = 'Rebecca R. Murphy',
  author_email = 'rrm33@cam.ac.uk',
  url='http://pypi.python.org/pypi/nxrepair/',
        license='LICENSE',
    long_description=open('README').read(),
    install_requires=[
        "scipy >= 0.9.0",
        "numpy >= 1.6.1",
        "matplotlib >= 1.1.1rc"
    ],
  )