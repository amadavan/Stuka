# Copyright 2015 The TensorFlow Authors. All Rights Reserved.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# ==============================================================================
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import fnmatch
import os
import re
import sys

from setuptools import Command
from setuptools import find_packages
from setuptools import setup
from setuptools.command.install import install as InstallCommandBase
from setuptools.dist import Distribution

REQUIRED_PACKAGES = [
    'numpy >= 1.14.5, < 2.0',
    'scipy',
]


class BinaryDistribution(Distribution):

    def has_ext_modules(self):
        return True


class InstallCommand(InstallCommandBase):
    """Override the dir where the headers go."""

    def finalize_options(self):
        ret = InstallCommandBase.finalize_options(self)
        self.install_headers = os.path.join(self.install_purelib, 'include')
        self.install_lib = self.install_platlib
        return ret


class InstallHeaders(Command):
    """Override how headers are copied.
    The install_headers that comes with setuptools copies all files to
    the same directory. But we need the files to be in a specific directory
    hierarchy for -I <include_dir> to work correctly.
    """
    description = 'install C/C++ header files'

    user_options = [('install-dir=', 'd',
                     'directory to install header files to'),
                    ('force', 'f',
                     'force installation (overwrite existing files)'),
                    ]

    boolean_options = ['force']

    def initialize_options(self):
        self.install_dir = None
        self.force = 0
        self.outfiles = []

    def finalize_options(self):
        self.set_undefined_options('install',
                                   ('install_headers', 'install_dir'),
                                   ('force', 'force'))

    def mkdir_and_copy_file(self, header):
        install_dir = os.path.join(self.install_dir, os.path.dirname(header))

        # Copy external code headers into tensorflow_core/include.
        # A symlink would do, but the wheel file that gets created ignores
        # symlink within the directory hierarchy.
        # NOTE(keveman): Figure out how to customize bdist_wheel package so
        # we can do the symlink.
        external_header_locations = [
            'include/stuka',
        ]
        for location in external_header_locations:
            if location in install_dir:
                extra_dir = install_dir.replace(location, '')
                if not os.path.exists(extra_dir):
                    self.mkpath(extra_dir)
                self.copy_file(header, extra_dir)

        if not os.path.exists(install_dir):
            self.mkpath(install_dir)
        return self.copy_file(header, install_dir)

    def run(self):
        hdrs = self.distribution.headers
        if not hdrs:
            return

        self.mkpath(self.install_dir)
        for header in hdrs:
            (out, _) = self.mkdir_and_copy_file(header)
            self.outfiles.append(out)

    def get_inputs(self):
        return self.distribution.headers or []

    def get_outputs(self):
        return self.outfiles


def find_files(pattern, root):
    """Return all the files matching pattern below root dir."""
    for dirpath, _, files in os.walk(root):
        for filename in fnmatch.filter(files, pattern):
            yield os.path.join(dirpath, filename)


headers = (
    list(find_files('*.h', 'include/stuka'))
)


setup(
    name='stuka',
    version='1.0.0',
    description='An optimization package leveraging fast C++ computation.',
    long_description='Exposes C++ library to solve optimization problems such as linear programs using a variety of '
                     'algorithms. It utilizes the Eigen library for fast computation allowing for rapid solution of '
                     'even large-scale problems.',
    # url='https://www.tensorflow.org/',
    # download_url='https://github.com/tensorflow/tensorflow/tags',
    author='Avinash Madavan',
    author_email='avinash.madavan@gmail.com',
    # Contained modules and scripts.
    packages=find_packages(),
    headers=headers,
    install_requires=REQUIRED_PACKAGES,
    # tests_require=REQUIRED_PACKAGES + TEST_PACKAGES,
    # Add in any packaged data.
    include_package_data=True,
    package_data={
        'stuka': ['stuka/stuka_extension.so'],
    },
    zip_safe=False,
    distclass=BinaryDistribution,
    cmdclass={
        'install_headers': InstallHeaders,
        'install': InstallCommand,
    },
    # PyPI package information.
    classifiers=[
        'Intended Audience :: Developers',
        'Intended Audience :: Education',
        'Intended Audience :: Science/Research',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Mathematics',
        'Topic :: Scientific/Engineering :: Artificial Intelligence',
        'Topic :: Software Development',
        'Topic :: Software Development :: Libraries',
        'Topic :: Software Development :: Libraries :: Python Modules',
    ],
    keywords='stuka optimization',
)