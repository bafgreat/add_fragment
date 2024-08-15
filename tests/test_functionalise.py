#!/usr/bin/python
from __future__ import print_function
__author__ = "Dr. Dinga Wonanke and Dr. Shabhajit Das"
__status__ = "production"

from add_fragment.io import load_system
from add_fragment.insert_fragment import functionalise_host


host_system = load_system('test_data/B-ditjoh-H.xyz')
sub_fragment1 = load_system('test_data/ac-R.xyz')


def test_functions(host_system, sub_fragment1):
    '''
    This function tests the functionalise_host function
    by loading a host and a sub-fragment, and then
    writing the aligned system to a file.
    '''
    functionalised_system = functionalise_host(host_system, sub_fragment1)

    assert len(functionalised_system) == 828


    # functionalised_system.write('combined_system.xyz')
# test_functions(host_system, sub_fragment1)