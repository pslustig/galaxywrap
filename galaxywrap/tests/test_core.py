from galaxywrap import core


def test_keywordtranslator():
    assert core.keywordtranslator.to_python('RE') == 'r'
    assert core.keywordtranslator.to_python('RA') == 'ratio'
    assert core.keywordtranslator.to_python('MAG') == 'mag'

    assert core.keywordtranslator.to_galfit('r') == 'RE'
    assert core.keywordtranslator.to_galfit('ratio') == 'RA'
    assert core.keywordtranslator.to_galfit('mag') == 'MAG'


def test_keywordtranslator_instance():
    trans = core.keywordtranslator()
    assert trans.to_python('RE') == 'r'
    assert trans.to_python('RA') == 'ratio'
    assert trans.to_python('MAG') == 'mag'

    assert trans.to_galfit('r') == 'RE'
    assert trans.to_galfit('ratio') == 'RA'
    assert trans.to_galfit('mag') == 'MAG'
