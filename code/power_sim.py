"""Monte-Carlo power for a sorted pooled CRISPR screen.

Noise model is anchored on real replicate data (GEO GSE226977, genome-wide
KO screen, 3 replicates/arm). Total count overdispersion decomposes as

    phi_total = phi_bio + 1/cells_per_guide + 1/reads_per_guide

where phi_bio is what remains of the measured endpoint dispersion after
removing the sequencing (Poisson) and cell-sampling terms of that experiment.
"""
import numpy as np

PHI_BIO = 0.288          # from GSE226977 endpoint arm, see notebook cell
EFF_A, EFF_B = 2.0, 1.33  # Beta(a,b) guide-efficacy multiplier, mean ~0.6


def phi_total(cells, reads, phi_bio=None):
    # read the module global at call time so sensitivity sweeps take effect
    pb = PHI_BIO if phi_bio is None else phi_bio
    return pb + 1.0 / cells + 1.0 / reads


def nb_counts(rng, mean, phi, size):
    """Negative-binomial draws with mean `mean` and dispersion `phi`."""
    if phi <= 0:
        return rng.poisson(mean, size)
    r = 1.0 / phi
    p = r / (r + mean)
    return rng.negative_binomial(r, p, size)


def simulate(rng, n_genes, guides, reps, cells, reads, lfc):
    """Return gene-level mean log2 fold change for `n_genes` genes."""
    phi = phi_total(cells, reads)
    eff = rng.beta(EFF_A, EFF_B, size=(n_genes, guides, 1))
    ctrl = nb_counts(rng, reads, phi, (n_genes, guides, reps)).astype(float)
    treat_mean = reads * np.power(2.0, -lfc * eff)
    treat = nb_counts(rng, treat_mean, phi, (n_genes, guides, reps)).astype(float)
    l2 = np.log2((treat + 1) / (ctrl + 1))
    return l2.mean(axis=(1, 2))


def power_at_fdr(rng, guides, reps, cells, reads, lfc, n_null=4000,
                 n_alt=2000, fdr=0.05):
    """Empirical power: null genes set the threshold, effect genes are scored."""
    null = simulate(rng, n_null, guides, reps, cells, reads, 0.0)
    alt = simulate(rng, n_alt, guides, reps, cells, reads, lfc)
    # one-sided (depletion from the high bin): p from the null ECDF
    order = np.sort(null)
    p_alt = np.searchsorted(order, alt, "right") / n_null
    p_null = np.searchsorted(order, null, "right") / n_null
    allp = np.concatenate([p_alt, p_null])
    lab = np.concatenate([np.ones(n_alt, bool), np.zeros(n_null, bool)])
    o = np.argsort(allp)
    allp, lab = allp[o], lab[o]
    m = allp.size
    bh = allp * m / np.arange(1, m + 1)
    bh = np.minimum.accumulate(bh[::-1])[::-1]
    sig = bh <= fdr
    return lab[sig].sum() / n_alt if sig.any() else 0.0
