# Federated k-Means over Networks via GTVMin — Numerical Experiments

Self-contained companion code for

> X. Yang, S. Rastelli, and A. Jung,
> **"Federated k-Means over Networks via Generalized Total Variation
> Minimization,"** submitted to *IEEE Transactions on Signal
> Processing*, 2026.

The paper formulates federated (hard) k-means clustering over a network
of nodes — each holding a private local dataset — as an instance of
generalized total variation minimization (GTVMin), with two
permutation-invariant couplings of the local cluster centroids:
nearest-centroid assignment (**FedKM-NC**) and optimal transport
assignment (**FedKM-OT**).

## Scripts

Each script is **fully self-contained** (algorithms, baselines, metrics
and data generation are embedded; nothing is imported across files) and
reproduces one experiment of Section V:

| Script | Paper artifact | Full runtime* |
|---|---|---|
| `exp1_effect_of_alpha.py` | Fig. 1, Table III — effect of the GTV weight α | 1–2 h |
| `exp2_samplesize_connectivity.py` | Sec. V-C — sample-size and connectivity sweeps | 2–4 h |
| `exp3_regime_check.py` | Table I — guarantee-covered update regimes of Prop. 1 | 30–60 min |
| `exp4_small_scale_comparison.py` | Table II — comparison with Distributed k-means and DGC-KM | 1–2 h |
| `exp5_communities.py` | Table IV — two weakly linked node communities | 20–40 min |
| `exp6_cifar10.py` | Fig. 2, Table V — CIFAR-10: convergence, communication, runtime | 3–5 h |

*single-threaded (the scripts pin BLAS to one thread so the wall-clock
numbers are comparable across methods, as in the paper).

Every script accepts:

```
--outdir DIR    output directory (default: outputs/)
--seeds N       number of random seeds (default: 10, as in the paper)
--quick         minutes-scale smoke test (small T, fewer seeds)
```

Example:

```bash
pip install numpy scipy scikit-learn networkx pandas
python3 exp5_communities.py --quick          # smoke test
python3 exp5_communities.py                  # full run (paper settings)
```

Outputs are CSV files (per-seed results and mean ± 95% CI aggregates)
plus, where applicable, the LaTeX table bodies used in the paper.

## Baselines implemented

* **Local / Central k-means** — scikit-learn KMeans per node / on the
  pooled data (Central defines the reference W* of the ξ₁ metric).
* **Distributed k-means** — Forero, Cano, Giannakis (IEEE JSTSP 2011):
  ADMM on the label-wise consensus formulation (optional
  shared-initialization variant).
* **DGC-KM_ρ** — the k-means instance of the distributed gradient
  clustering framework of Armacki, Bajović, Jakovetić, Kar
  (IEEE TSP 2025), implemented from their Algorithm 1 / eq. (7) with
  the Theorem-1-compatible step size.

## Notes

* `exp6_cifar10.py` downloads the CIFAR-10 python tarball (~170 MB)
  to `~/.cache/cifar10_repro` on first use.
* The 10 seeds are fixed inside the scripts
  (`np.random.default_rng(42).integers(0, 100000, 10)`) — identical to
  the runs reported in the paper.
* Metrics: Global Centroid Deviation ξ₁ and Consensus Variation ξ₂,
  both Hungarian-aligned (permutation-invariant), as defined in
  Sec. V-B of the paper.
