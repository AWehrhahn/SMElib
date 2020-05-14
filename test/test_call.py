from os.path import dirname, join

import pytest

from sme_synth import SME_DLL
from cwrapper import get_lib_name

@pytest.fixture
def libfile():
    return join(dirname(__file__), "../lib/", get_lib_name())

@pytest.fixture
def dll(libfile):
    return SME_DLL(libfile)

def test_simple_call(dll):
    version = dll.SMELibraryVersion()
    assert isinstance(version, str)

def test_call_with_input(dll):
    dll.InputWaveRange(5000, 6000)
