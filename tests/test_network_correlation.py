from __future__ import annotations

import numpy as np

from scripts.run_bo2023_network_correlation_validation import select_group_discriminative_genes


def test_network_discriminative_gene_selection_identifies_group_signal():
    values = np.asarray(
        [
            [9.0, 8.5, 1.0, 1.5],
            [3.0, 7.0, 3.0, 7.0],
            [1.0, 1.5, 8.5, 9.0],
        ]
    )
    selected, _ = select_group_discriminative_genes(
        values,
        ["A", "B"],
        {"A": np.asarray([0, 1]), "B": np.asarray([2, 3])},
        top_n=2,
    )
    assert set(selected) == {0, 2}
