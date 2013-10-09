#!/usr/bin/env python

#-----------------------------------------------------------------------------
#       Copyright (C) 2013 The Mosaic Development Team
#
#  Distributed under the terms of the BSD License.  The full license is in
#  the file LICENSE.txt, distributed as part of this software.
#-----------------------------------------------------------------------------

from distutils.core import setup, Command
from unittest import TextTestRunner
import copy
import os
import sys

package_dir = "lib"
script_dir = "scripts"
test_dir = "tests"
root_dir = os.path.split(os.path.join(os.getcwd(), sys.argv[0]))[0]
commands = {}

##################################################
# Run unit tests

class TestCommand(Command):
    """Runs the unit tests"""

    user_options = []
    def initialize_options(self):
        pass
    def finalize_options(self):
        pass
    def run(self):
        sys.path.insert(0, os.path.join(root_dir, package_dir))
        sys.path.insert(0, os.path.join(root_dir, test_dir))
        import all_tests
        t = TextTestRunner(verbosity=2)
        t.run(all_tests.suite())

commands['test'] = TestCommand

##################################################
# Build documentation if sphinx is installed

try:

    from sphinx.setup_command import BuildDoc as SphinxBuildDoc

    class BuildDoc(SphinxBuildDoc):
        def run(self):
            # Make sure the python path is pointing to the newly built
            # code so that the documentation is built on this and not a
            # previously installed version.
            build = self.get_finalized_command('build')
            sys.path.insert(0, os.path.abspath(build.build_lib))
            SphinxBuildDoc.run(self)
            sys.path.pop(0)

    commands['build_sphinx'] = BuildDoc

except ImportError:
    pass

##################################################


with open('README.txt') as file:
    long_description = file.read()

setup(name='pyMosaic',
      version='0.9',
      description='MOlecular SimulAtion Interchange Conventions',
      long_description=long_description,
      author='Konrad Hinsen',
      author_email='research@khinsen.fastmail.net',
      url='https://bitbucket.org/molsim/mosaic',
      package_dir = {'': package_dir},
      packages=['mosaic', 'mosaic_pdb'],
      scripts=[os.path.join(script_dir, s) for s in os.listdir(script_dir)],
      cmdclass = commands,
     )
