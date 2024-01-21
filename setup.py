#!/usr/bin/env python

"""
Setup script for the AstroPhotography application.
"""

from setuptools import setup, find_packages

with open('README.md') as readme_file:
    readme = readme_file.read()

with open('RELEASE_NOTES.md') as history_file:
    history = history_file.read()

requirements = [ ]

test_requirements = ['pytest>=3', ]

setup(
    author="Dave Strickland",
    author_email='dave.strickland@gmail.com',
    python_requires='>=3.6',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
    ],
    description="Python workflow f for processing astrophotographic images from amateur telescopes and/or digital cameras",
    entry_points={
        'console_scripts': [
            'dksraw=AstroPhotography.cli:main',
        ],
    },
    install_requires=requirements,
    license="GNU General Public License v3",
    long_description=readme + '\n\n' + history,
    include_package_data=True,
    keywords='AstroPhotography',
    name='AstroPhotography',
    packages=find_packages(include=['AstroPhotography', 'AstroPhotography.*']),
    test_suite='test',
    tests_require=test_requirements,
    url='https://github.com/DaveStrickland/AstroPhotography',
    version='0.5.1',
    zip_safe=False,
)
