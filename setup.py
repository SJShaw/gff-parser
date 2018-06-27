"""Setuptools magic to install gffparser."""
import os
import sys
from setuptools import setup
from setuptools.command.test import test as TestCommand


INSTALL_REQUIRES = [
    'biopython >= 1.71',
    'helperlibs',
]

TESTS_REQUIRE = [
    'pytest',
    'coverage',
    'pylint',
]


class PyTest(TestCommand):
    """Allow running tests via python setup.py test."""

    def finalize_options(self):
        """Test command magic."""
        TestCommand.finalize_options(self)
        self.test_args = []
        self.test_suite = True

    def run_tests(self):
        """Run tests."""
        import pytest
        errcode = pytest.main(self.test_args)
        sys.exit(errcode)


setup(
    name="gffparser",
    version="0.0",
    packages=['gffparser'],
    author='antiSMASH development team',
    author_email='antismash@secondarymetabolites.org',
    description='The antibiotics and Secondary Metabolites Analysis Shell.',
    long_description=open(os.path.join(os.path.dirname(__file__), 'README.md')).read(),
    install_requires=INSTALL_REQUIRES,
    tests_require=TESTS_REQUIRE,
    cmdclass={'test': PyTest},
    url='https://github.com/SJShaw/gff-parser',
    license='GNU Affero General Public License v3 or later (AGPLv3+)',
    classifiers=[
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: GNU Affero General Public License v3 or later (AGPLv3+)',
        'Operating System :: OS Independent',
    ],
    extras_require={
        'testing': TESTS_REQUIRE,
    },
)
