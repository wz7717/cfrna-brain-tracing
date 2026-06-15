from __future__ import annotations

import numpy as np

from scripts.run_bo2023_discriminative_correlation_validation import (
    select_fold_discriminative_genes,
)


def test_select_fold_discriminative_genes_uses_region_stable_signal():
    values = np.asarray(
        [
            [8.0, 8.1, 7.9, 1.0, 1.1, 0.9],
            [4.0, 1.0, 7.0, 4.0, 1.0, 7.0],
            [1.0, 1.1, 0.9, 8.0, 8.1, 7.9],
        ]
    )
    selected, audit = select_fold_discriminative_genes(
        values,
        ["R1", "R2"],
        {"R1": np.asarray([0, 1, 2]), "R2": np.asarray([3, 4, 5])},
        heldout_idx=-1,
        top_n=2,
        min_region_train_samples=3,
    )
    assert set(selected) == {0, 2}
    assert float(audit["fisher_score"].min()) > 1.0
