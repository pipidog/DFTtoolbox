from setuptools import setup, find_packages

# with open('C:/Users/pipidog/Dropbox/Code/DFTtoolbox/README.md') as file:
    # long_description = file.read()

setup(
    name = 'DFTtoolbox',
    version = '1.0',
    packages = ['DFTtoolbox'],
    description = 'A toolbox to initialize or postpocess several codes',
    scripts = [],
    license='MIT',
    author = 'pipidog',
    author_email = 'pipidog@gmail.com',
    url = 'https://github.com/pipidog/DFTtoolbox',
    download_url = 'https://github.com/pipidog/DFTtoolbox/archive/v1.0.tar.gz',
    keywords = ['density-functional-theory','qantum-espresso','elk','abinit'],
    classifiers = ['Topic :: Scientific/Engineering :: Physics']
)