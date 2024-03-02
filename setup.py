# Copyright 2024, Lorenzo Martelli
# Authors: Lorenzo Martelli, Igor Andriyash
# License: GPL3

from setuptools import setup, find_packages

try:
    import pypandoc
    long_description = pypandoc.convert( './README.md', 'rst')
except (ImportError, RuntimeError):
    long_description = open('./README.md').read()

with open('requirements.txt') as f:
    install_requires = [ line.strip('\n') for line in f.readlines() ]

setup(
    name='FBPIC-EWP',
    version='v1.0.0',
    description=long_description,
    author=['Lorenzo Martelli', 'Igor Andriyash'],
    maintainer='Lorenzo Martelli',
    maintainer_email='lorenzo.martelli295@gmail.com',
    license='GPL3',
    packages=find_packages('.'),
    package_data={"": ['*']},
    tests_require=[],
    cmdclass={},
    install_requires=install_requires,
    include_package_data=True,
    platforms='any',
    url='https://github.com/laumrt/FBPIC-EWP',
    classifiers=[
        'Programming Language :: Python',
        'Natural Language :: English',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'Operating System :: OS Independent',
        'Topic :: Scientific/Engineering :: Physics',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5'],
    zip_safe=False
)
