from os.path import dirname, join

import pytest

from sme_synth import SME_DLL
from cwrapper import get_lib_name

@pytest.fixture
def libfile():
    return join(dirname(__file__), "../lib/", get_lib_name())

@pytest.fixture
def datadir():
    return join(dirname(__file__), "../share/libsme")

@pytest.fixture
def dll(libfile, datadir):
    # We set the data directory explicitly in case
    # the library file has been loaded from a different location
    return SME_DLL(libfile, datadir)

def test_simple_call(dll):
    version = dll.SMELibraryVersion()
    assert isinstance(version, str)

def test_call_with_input(dll):
    dll.InputWaveRange(5000, 6000)

if __name__ == "__main__":
    libfile = join(dirname(__file__), "../lib/", "libsme-5.dll")
    datadir = join(dirname(__file__), "../share/libsme")
    dll = SME_DLL(libfile, datadir)
    version = dll.SMELibraryVersion()
    print(version)