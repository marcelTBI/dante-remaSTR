import dante_remastr_standalone
from dante_remastr_standalone import Motif


def test_augment_nomenclature() -> None:
    print()

    motif: Motif
    motif = Motif("chr4:g.3074878_3074973AGC[21]CGC[11]", "HD")
    motif.monoallelic = False

    hapl: list[str] = ['17', '11']
    nomenclatures: list[list[tuple[int, int, str]]] = [
        [
            (17, 4, 'AGC[15]AAC[1]AGC[1]'),
            (17, 1, 'AGC[4]ATC[1]AGC[7]AAC[1]AGC[2]AAC[1]AGC[1]')
        ], [
            (10, 5, 'CGC[1]CAC[1]CGC[7]CTCCTC[1]'),
            (9, 1, 'CGC[1]CATCGT[1]CGC[6]CTCCTT[1]')
        ]
    ]

    result = dante_remastr_standalone.augment_nomenclature(motif, hapl, nomenclatures, 0.75)
    aug_nom, nom_len, nomenclatures, errors = result
    print(aug_nom)
    print(nom_len)
    print(nomenclatures)
    print(errors)

    # expected results
    # ['AGC[15]AAC[1]AGC[1]', 'CGC[1]CAC[1]CGC[7]CTCCTC[1]']
    # 84
    # [[(17, 1, 'AGC[15]AAC[1]AGC[1]'), (17, 1, 'AGC[4]ATC[1]AGC[7]AAC[1]AGC[2]AAC[1]AGC[1]')], [(10, 1, 'CGC[1]CAC[1]CGC[7]CTCCTC[1]'), (9, 1, 'CGC[1]CATCGT[1]CGC[6]CTCCTT[1]')]]
    # []
