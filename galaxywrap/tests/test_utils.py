from galaxywrap import utils


def test_isiterable_True():
    assert utils.isiterable((3, 4))


def test_isiterable_False():
    assert not utils.isiterable(3)
