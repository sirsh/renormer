import os, io

from setuptools import setup, find_packages

setup(
    name='renormer',
    version='0.1',
    author='Sirsh',
    author_email='amartey@gmail.com',
    license='MIT',
    url='git@github.com:mrsirsh/renormer.git',
    keywords='RG, python',
    description='renormer description',
    long_description=('renormer description'),
    packages=['rg'],
    test_suite='nose.collector',
    tests_require=['nose'],
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Developers',
        'License :: MIT',
        'Operating System :: OS Independent',
        'Topic :: Software Development :: Libraries :: Python Modules',
        'Topic :: Communications :: Chat',
        'Topic :: Internet',
        'Programming Language :: Python',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.5',
    ],
)


