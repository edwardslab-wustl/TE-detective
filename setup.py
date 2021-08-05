#!/usr/bin/env python

from setuptools import setup

def readme():
    with open('README.md') as f:
        return f.read()

setup(name='TEdetective',
        version='0.1.0',
        description="s",
        long_description=readme(),
        url='http://github.com/edwardslab-wustl/TE-detective',
        author='Manoj Singh and John Edwards',
        author_email="jredwards@wustl.edu",
        license='GPLv3',
        packages=['TEdetective'],
        classifiers=[
            'Development Status :: 4 - Beta',
            'Intended Audience :: Science/Research',
            'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
            'Operating System :: POSIX',
            'Topic :: Scientific/Engineering :: Bio-Informatics',
            'Programming Language :: Python :: 2',
            'Programming Language :: Python :: 2.6',
            'Programming Language :: Python :: 2.7',
            'Programming Language :: Python :: 3',
            'Programming Language :: Python :: 3.2',
            'Programming Language :: Python :: 3.3',
            'Programming Language :: Python :: 3.4',
        ],
        install_requires=[
            'numpy',
            'pysam',
            'biopython',
        ],
        python_requires='>=2.6, !=3.0.*, !=3.1.*, !=3.2.*, <4',
        entry_points = {
            "console_scripts": [
                "TE_detective=TEdetective.entry:main",
                "TE_detective_eval=TEdetective.eval_entry:main",
                ],
            },
        include_package_data=True,
        zip_safe=False)


#packages=find_packages(exclude=['contrib', 'docs', 'tests*']),
