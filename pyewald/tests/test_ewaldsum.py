from .. import ewaldsum

from numpy import array, diag, sqrt
from numpy.linalg import det, inv

from numpy.testing import assert_almost_equal, assert_no_warnings


def test_madelung_NaCL():
    """
    NaCl
    """
    cell = diag((1., 1., 1.))
    invcell = inv(cell).T
    r = array([[0., 0., 0.],
               [.5, .5, .0],
               [.5, .0, .5],
               [.0, .5, .5],
               [.5, 0., 0.],
               [.0, .5, .0],
               [.0, .0, .5],
               [.5, .5, .5]])
    q = array([1, 1, 1, 1, -1, -1, -1, -1])
    i = 0
    area = abs(det(cell))

    V = ewaldsum.potential(i, r, q, cell, area, invcell)
    V = - 0.5 * V
    refV = 1.7475645946331822
    assert_almost_equal(refV, V, err_msg='problem with NaCl')


def test_madelung_CsCL():
    """
    CsCl
    """
    cell = diag((1., 1., 1.))
    r = array([[0., 0., 0.],
              [.5,  .5,  .5]])
    invcell = inv(cell).T
    i = 0
    q = array([1, -1])
    area = abs(det(cell))
    V = ewaldsum.potential(i, r, q, cell, area, invcell)
    V *= - .5 * sqrt(3.)
    refV = 1.76267477307099
    assert_almost_equal(refV, V, err_msg='problem with CsCl')


def test_madelung_ZnS():
    """
    ZnS
    """
    cell = array([[.0, .5, .5],
                  [.5, .0, .5],
                  [.5, .5, .0]])
    invcell = inv(cell).T
    r = array([[0., 0., 0.],
              [.25, .25, .25]])
    i = 0
    alpha = 5./sqrt((cell**2).sum(1).max())
    area = abs(det(cell))
    q = array([1, -1])
    V = ewaldsum.potential(i, r, q, cell, area, invcell, alpha)
    V = - 0.25 * sqrt(3) * V
    refV = 1.63805505338879
    assert_almost_equal(refV, V, err_msg='problem with ZnS')


def test_madelung_ZnSB4():
    """
    ZnSB4
    """
    a = 1.
    c = sqrt(8./3.)*a
    u = 3./8.
    cell = array([[.5*a, -0.5*sqrt(3)*a, 0],
                  [.5*a, 0.5*sqrt(3)*a, 0],
                  [0., 0., c]])
    invcell = inv(cell).T
    area = abs(det(cell))
    alpha = 5./sqrt((cell**2).sum(1).max())
    r = array([[.5*a, .5/sqrt(3)*a, 0.],
               [.5*a, -.5/sqrt(3)*a, 0.5*c],
               [.5*a, .5/sqrt(3)*a, u*c],
               [.5*a, -.5/sqrt(3)*a, (.5+u)*c]])
    i = 0
    q = array([1., 1., -1, -1])
    V = ewaldsum.potential(i, r, q, cell, area, invcell, alpha)
    V = - 1. * u * c * V
    refV = 1.64132162737
    assert_almost_equal(refV, V, err_msg='problem with ZnSB4')


def test_madelung_random():
    """
    random
    """
    cell = diag((1., 1., 1.))
    r = array([[0., 0., 0.],
              [.234, .978789, .05]])
    invcell = inv(cell).T
    i = 0
    q = array([1, -1])
    area = abs(det(cell))
    assert_no_warnings(ewaldsum.potential, i, r, q, cell, area, invcell)
