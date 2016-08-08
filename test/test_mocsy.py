from __future__ import (absolute_import, division, print_function)

import pytest
import numpy as np

import mocsy


@pytest.fixture
def scalar_variables():
    """
    Functions return 1-dimensional numpy arrays.
    Scalar inputs return length-1 arrays.
    DATA input: DIC and ALk in mol/kg, in situ temperature, pressure.
    """

    return mocsy.mvars(temp=18,
                       sal=35,
                       alk=2300.e-6,
                       dic=2000.e-6,
                       sil=0,
                       phos=0,
                       patm=1,
                       depth=100,
                       lat=0,
                       optcon='mol/kg',
                       optt='Tinsitu',
                       optp='db',
                       optb='u74',
                       optk1k2='l',
                       optkf='dg',
                       optgas='Pinsitu')

def test_return_12():
    ret = scalar_variables()
    assert len(ret) == 12

def test_return_scalar():
    ret = scalar_variables()
    for var in ret:
        assert len(var) == 1

def test_return_real():
    ret = scalar_variables()
    for var in ret:
        assert np.isreal(var)

def test_known_values():
    ret = scalar_variables()
    pH = 8.14892578
    pco2 = 312.28662109
    fco2 = 300.68057251
    co2 = 1.01729711e-05
    hco3 = 0.00177952
    co3 = 0.00021031
    OmegaA = 3.19940853
    OmegaC = 4.94189167
    BetaD = 9.68977737
    DENis = 1025.71105957
    p = 100.0
    Tis = 18.0
    known = pH, pco2, fco2, co2, hco3, co3, OmegaA, OmegaC, BetaD, DENis, p, Tis

    np.testing.assert_allclose(known, np.array(ret).ravel(), rtol=1e-05)
