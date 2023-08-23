from setuptools import setup, find_packages

version='1.6.3'
setup(
    name = 'DFTtoolbox',
    version = version,
    packages = ['DFTtoolbox'],
    description = 'A toolbox to initialize or postpocess several DFT codes',
    scripts = [],
    license='MIT',
    author = 'pipidog',
    author_email = 'pipidog@gmail.com',
    url = 'https://github.com/pipidog/DFTtoolbox',
    download_url = 'https://github.com/pipidog/DFTtoolbox/archive/v'+version+'.tar.gz',
    keywords = ['density-functional-theory','qantum-espresso','elk','abinit'],
    classifiers = ['Topic :: Scientific/Engineering :: Physics'],
    install_requires=['numpy','matplotlib']
)
