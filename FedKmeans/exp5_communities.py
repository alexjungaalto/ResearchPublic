#!/usr/bin/env python3
"""Two weakly linked node communities (Table IV of the paper).

Companion code for:
  X. Yang, S. Rastelli, and A. Jung, "Federated k-Means over
  Networks via Generalized Total Variation Minimization,"
  submitted to IEEE Transactions on Signal Processing, 2026.

Two communities of 5 nodes with translated data mixtures,
SBM graph (p_in=0.8, p_out in {0, 0.1}); deviation
xi_1^comm from the community-wise centralized solutions for
FedKM-NC/OT (alpha in {0.01,0.1,1}), Distributed k-means,
DGC-KM_1 and Local k-means; T=200, P=10, m_i=50.

Approximate full-settings runtime: 20-40 min single-threaded.
Use --quick for a minutes-scale smoke test.

Dependencies: numpy, scipy, scikit-learn, networkx, pandas
    pip install numpy scipy scikit-learn networkx pandas
"""

import os

# Pin BLAS to a single thread, as in Xu's notebooks: (a) the timing
# numbers in Table I of the paper were measured single-threaded, and
# (b) sklearn's KMeans warns about oversubscription with joblib.
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import argparse
import pickle
import tarfile
import time
import urllib.request
from pathlib import Path

import networkx as nx
import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment
from scipy.stats import t as t_dist
from sklearn.cluster import KMeans
from sklearn.datasets import make_blobs
from sklearn.decomposition import PCA
from sklearn.metrics import pairwise_distances_argmin

# The 10 seeds used throughout Sec. V, generated exactly as in
# NumExpXu/FedKMeans_Synthetic.ipynb:
#   [8925, 77395, 65457, 43887, 43301, 85859, 8594, 69736, 20146, 9417]
# Seed 20146 (index 8) is the run shown in Fig. centroid_traj.
SEEDS = np.random.default_rng(42).integers(0, 100000, 10)

# epsilon in the denominators of steps Alg1-M / Alg2-M.  The paper
# requires epsilon > 0 as an input of Algorithms 1-2 (Proposition 1 is
# stated for the limit epsilon = 0 with a keep-previous-value
# convention for empty clusters; the code realizes that convention
# explicitly via the `mask` below, so epsilon only guards against
# division by tiny positive denominators).
EPS = 1e-12

# Stopping threshold delta of Algorithms 1-2: leave the inner loop as
# soon as one pass decreases the local objective by less than delta
# (matches tol=1e-6 in the notebooks).
DELTA = 1e-6


# ----------------------------------------------------------------------
# Data and graph generation (conventions of FedKMeans_Synthetic.ipynb)
# ----------------------------------------------------------------------

def generate_centers(K=3, dim=2, center_box=(-10, 10), min_sep=5.0, seed=None):
    """Random cluster centers with pairwise separation >= min_sep.

    Sec. V-C: "The minimum separation between cluster centroids is set
    to 5."  Rejection sampling from the uniform box, identical to the
    notebook helper (uses the legacy RandomState so that the same seed
    reproduces the notebook's centers bit-for-bit).
    """
    rng = np.random.RandomState(seed)
    centers = []
    while len(centers) < K:
        cand = rng.uniform(center_box[0], center_box[1], size=dim)
        if all(np.linalg.norm(cand - c) >= min_sep for c in centers):
            centers.append(cand)
    return np.array(centers)


def gen_synthetic(name, n_samples, K=3, seed=None):
    """One of the three 2-D geometries of Sec. V-C.

    iso    D_iso:   isotropic blobs, sigma = 1.0 for every cluster.
    aniso  D_aniso: the same blobs pushed through the fixed
           rotation-scaling A = [[0.6,-0.6],[-0.4,0.8]] (elliptical
           clusters; stresses the isotropic k-means loss).
    var    D_var:   cluster-specific sigmas {1.0, 2.5, 0.5} (unequal
           spreads; stresses the equal-variance assumption).

    Returns (X, y) with y the ground-truth component labels (used only
    for evaluation, never by the algorithms).
    """
    centers = generate_centers(K=K, dim=2, min_sep=5.0, seed=seed)
    if name == "iso":
        X, y = make_blobs(n_samples=n_samples, centers=centers,
                          cluster_std=1.0, random_state=seed)
    elif name == "aniso":
        X, y = make_blobs(n_samples=n_samples, centers=centers,
                          cluster_std=1.0, random_state=seed)
        X = X.dot(np.array([[0.6, -0.6], [-0.4, 0.8]]))
    elif name == "var":
        X, y = make_blobs(n_samples=n_samples, centers=centers,
                          cluster_std=[1.0, 2.5, 0.5][:K], random_state=seed)
    else:
        raise ValueError(name)
    return X, y


def partition_uniform(X, num_nodes, seed=None):
    """Uniformly random split of X into num_nodes local datasets.

    Sec. V-C: "The datasets are uniformly split across n = 10 nodes."
    A random permutation chopped into near-equal parts; local sample
    sizes m_i differ by at most 1, so the paper's equal-sample-size
    equivalence (Sec. II) holds essentially exactly.
    """
    rng = np.random.default_rng(seed)
    idx = np.arange(X.shape[0])
    rng.shuffle(idx)
    return [X[p] for p in np.array_split(idx, num_nodes)]


def build_er_graph(num_nodes, p, seed=None):
    """Erdos-Renyi G(n, p) adjacency matrix (0/1, symmetric, zero
    diagonal).  This is the FL network of the synthetic and CIFAR
    experiments; p is the edge probability swept in Fig. ablation
    row 3.  Note p=0 yields the empty graph (no communication) and the
    graph is NOT forced to be connected -- isolated nodes are handled
    by the n_eff convention of xi_2, exactly as in the paper.
    """
    rng = np.random.default_rng(seed)
    upper = rng.random((num_nodes, num_nodes)) < p
    A = np.triu(upper, k=1).astype(float)
    return A + A.T


def build_rgg_graph(num_nodes=30, target_degree=4.867, seed=0):
    """Connected random geometric graph statistically matched to the
    benchmark graph of Forero et al., Sec. V-A (n=30 nodes, average
    degree 4.867, algebraic connectivity 0.07057).

    Their exact graph realization is not published, so we sweep the
    connection radius upward and return the first connected RGG whose
    average degree reaches the target -- a graph with the same
    n/degree profile, not the identical instance.  The realized degree
    is printed by the benchmark experiment for transparency.
    """
    for radius in np.linspace(0.15, 0.6, 200):
        G = nx.random_geometric_graph(num_nodes, radius, seed=seed)
        deg = 2 * G.number_of_edges() / num_nodes
        if nx.is_connected(G) and deg >= target_degree:
            return nx.to_numpy_array(G), deg
    raise RuntimeError("no connected RGG with the target degree found")


# ----------------------------------------------------------------------
# Metrics (Hungarian-aligned, eqs. central_error / nodes_error / local_cost)
# ----------------------------------------------------------------------
# The paper (end of Sec. V-B): "we always use the Hungarian method to
# permute the cluster centroid labels ... before computing xi_1 and
# xi_2".  Without this alignment, a mere relabeling of centroids would
# inflate both metrics; with it, they measure genuine geometric
# deviation.

def match_clusters(mu_ref, mu_other):
    """Align mu_other to mu_ref by the Hungarian method (Kuhn 1955).

    Solves the linear assignment problem on the K x K matrix of
    squared Euclidean distances between centroids and returns
    (mu_other reordered to best match mu_ref row-wise, permutation).
    """
    D = ((mu_ref[:, None, :] - mu_other[None, :, :]) ** 2).sum(axis=2)
    # Guard against NaN/inf centroids (e.g. a diverged baseline run):
    # replace by a large finite cost so the assignment stays defined.
    D = np.nan_to_num(D, nan=1e15, posinf=1e15)
    row, col = linear_sum_assignment(D)
    perm = np.empty(mu_ref.shape[0], dtype=int)
    perm[row] = col
    return mu_other[perm], perm


def gcd(mu_central, Ws):
    """xi_1: Global Centroid Deviation, eq. (central_error):

        xi_1 = (1 / (n K)) sum_i || W^(i) - W* ||_F^2

    with W* the centralized k-means solution on the pooled dataset and
    each W^(i) Hungarian-aligned to W* first.  Lower = closer to the
    ideal centralized clustering.
    """
    n, K = len(Ws), Ws[0].shape[0]
    total = 0.0
    for W in Ws:
        aligned, _ = match_clusters(mu_central, W)
        total += np.sum((mu_central - aligned) ** 2)
    return total / (n * K)


def cv(Ws, A):
    """xi_2: Consensus Variation, eq. (nodes_error):

        xi_2 = (1 / (K n_eff)) sum_{i: |N(i)|>0} (1/|N(i)|)
               sum_{i' in N(i)} || W^(i) - W^(i') ||_F^2

    where n_eff counts nodes with at least one neighbor (isolated
    nodes are skipped, so xi_2 is well-defined for sparse graphs, e.g.
    the p -> 0 end of the connectivity sweep).  Serves as an empirical
    proxy for the GTV term; lower = stronger cross-node consensus.
    """
    n, K = len(Ws), Ws[0].shape[0]
    total, n_eff = 0.0, 0
    for i in range(n):
        neigh = np.where(A[i] != 0)[0]
        if len(neigh) == 0:
            continue
        n_eff += 1
        d_i = 0.0
        for j in neigh:
            aligned, _ = match_clusters(Ws[i], Ws[j])
            d_i += np.sum((Ws[i] - aligned) ** 2)
        total += d_i / len(neigh)
    return total / (K * n_eff) if n_eff else 0.0


def local_distortion(node_data, Ws):
    """LD, eq. (local_cost): the metric of Forero et al. used in the
    benchmark experiment (Figs. kmpp_init*),

        LD = (1/2) sum_i sum_r min_c || x^(i,r) - w^(i,c) ||_2^2 .

    Sum of squared distances from each datapoint to its nearest LOCAL
    centroid.  Note the paper's caveat (Sec. V-B): because each node
    may overfit its own local dataset, federated methods can undercut
    the centralized LD -- which is why xi_1/xi_2 are preferred.
    """
    total = 0.0
    for X, W in zip(node_data, Ws):
        D = ((X[:, None, :] - W[None, :, :]) ** 2).sum(axis=2)
        total += D.min(axis=1).sum()
    return 0.5 * total


def gtv_discrepancy(Wi, Wj, ell):
    """The edge discrepancy d^(i,i') between two centroid matrices.

    ell = 1, eq. (equ-def_gtv), nearest-centroid assignment:
        component I  = sum over rows of Wi of the distance to the
                       nearest row of Wj      (D.min(axis=1))
        component II = sum over rows of Wj of the distance to the
                       nearest row of Wi      (D.min(axis=0))
        d = I + II.  Non-bijective: several centroids may share the
        same nearest counterpart.

    ell = 2, eq. (equ-def_gtv_perm), optimal transport assignment:
        d = min over permutations s of sum_c || Wi[c] - Wj[s(c)] ||^2,
        solved exactly by the Hungarian method.  Bijective by
        construction.

    Both vanish iff the two centroid SETS coincide (regardless of
    ordering) -- this permutation invariance is the paper's central
    design point.
    """
    D = ((Wi[:, None, :] - Wj[None, :, :]) ** 2).sum(axis=2)
    if ell == 1:
        return D.min(axis=1).sum() + D.min(axis=0).sum()
    row, col = linear_sum_assignment(D)
    return D[row, col].sum()


def gtvmin_objective(node_data, Ws, A, alpha, ell):
    """f_ell, the GTVMin objective eq. (equ_fedkmeans_gtvmin):

        f_ell(W^(1),...,W^(n)) = sum_i L_i(W^(i))
                                 + alpha * sum_{edges {i,i'}} d^(i,i')

    with the 1/m_i-normalized local k-means losses of
    eq. (equ_def_local_loss).  Used for the descent logging of
    experiment spec E2: Proposition 1 guarantees this value never
    increases across outer iterations (regime (i): P=1, any gamma;
    regime (ii): any P with the safeguard, gamma <= 1/2).
    """
    obj = 0.0
    for X, W in zip(node_data, Ws):
        D = ((X[:, None, :] - W[None, :, :]) ** 2).sum(axis=2)
        obj += D.min(axis=1).sum() / len(X)
    n = len(Ws)
    # Each undirected edge {i, i'} counted once, as in the paper.
    for i in range(n):
        for j in range(i + 1, n):
            if A[i, j] != 0:
                obj += alpha * gtv_discrepancy(Ws[i], Ws[j], ell)
    return obj


# ----------------------------------------------------------------------
# Algorithm 1: local solver via nearest-centroid assignment (FedKM-NC)
# ----------------------------------------------------------------------
# Solves eq. (equ_def_Jacobi_GTVMin) at device i: minimize the local
# k-means loss plus alpha times the ell=1 GTV terms to the (frozen)
# neighbor centroids, by a Lloyd-style alternation of assignment steps
# (E1, E2, E3) and a closed-form centroid update (M).

def _surrogate_F(W, X, neigh_mus, assign, fwd_list, rev_list, alpha):
    """F(W; pi_1), eq. (equ_surrogate_F): the convex surrogate of
    eq. (equ_def_Jacobi_GTVMin) obtained by FREEZING all assignment
    variables at pi_1 = ({khat_r}, {khat_{i',c'}}, {b_{i',c}}).

    Term by term:
      * (1/m_i) sum_r || W[khat_r] - x_r ||^2
            local loss with frozen datapoint assignments,
      * alpha * sum_{i'} sum_{c'} || W[khat_{i',c'}] - w^(i',c') ||^2
            component I with frozen reverse assignments (rev),
      * alpha * sum_{i'} sum_c  || W[c] - w^(i', b_{i',c}) ||^2
            component II with frozen forward assignments (fwd).

    F is a quadratic (hence convex) function of W that upper-bounds
    the true local objective and is tight at the point where pi_1 was
    computed -- the majorization used in the proof of Proposition 1.
    """
    val = np.sum((W[assign] - X) ** 2) / len(X)
    for mu_j, fwd, rev in zip(neigh_mus, fwd_list, rev_list):
        val += alpha * np.sum((W[rev] - mu_j) ** 2)        # component I
        val += alpha * np.sum((W - mu_j[fwd]) ** 2)        # component II
    return val


def local_update_nc(X, mu, neigh_mus, alpha, P=10, delta=DELTA, eps=EPS):
    """Algorithm 1 of the paper.  Returns (W, safeguard_triggered).

    Inputs mirror the algorithm's Require line: local dataset X
    (m_i x d), current centroids mu (K x d), the neighbors' centroid
    matrices neigh_mus (list of K x d), GTV weight alpha > 0, max
    number of passes P, stopping threshold delta, epsilon.
    """
    K, d = mu.shape
    m = max(len(X), 1)
    W_in = mu.copy()          # W_in <- W^(i): the round's input, kept
    W = mu.copy()             # for the safeguard comparison below
    pi1 = W1 = None           # pass-1 assignments / iterate storage
    prev_obj = None
    for tau in range(1, P + 1):
        # ---- (Alg1-E1) local assignment step -------------------------
        # khat_r <- argmin_c || w^(i,c) - x^(i,r) ||^2 for every local
        # datapoint r.
        assign = pairwise_distances_argmin(X, W)

        # Accumulators for the closed-form update (Alg1-M):
        #   num[c] = (1/m_i) sum_{r: khat_r = c} x_r
        #            + alpha * sum_{(i',c') in U_c} w^(i',c')   (comp. I)
        #            + alpha * sum_{i'} w^(i', b_{i',c})        (comp. II)
        #   den[c] = |{r: khat_r = c}| / m_i
        #            + alpha * (|U_c| + |N(i)|)
        fwd_list, rev_list = [], []
        num = np.zeros((K, d))
        den = np.zeros(K)
        n_c = np.bincount(assign, minlength=K).astype(float)
        for c in range(K):
            if n_c[c] > 0:
                num[c] = X[assign == c].sum(axis=0)
        num /= m                      # the 1/m_i normalization of L_i
        den += n_c / m
        for mu_j in neigh_mus:
            # K x K squared distances between own and neighbor centroids
            D = ((W[:, None, :] - mu_j[None, :, :]) ** 2).sum(axis=2)
            # ---- (Alg1-E3) forward map: for each LOCAL centroid c,
            # b_{i',c} = index of the nearest NEIGHBOR centroid
            # (solves component II of eq. equ_def_Jacobi_GTVMin).
            fwd = D.argmin(axis=1)
            # ---- (Alg1-E2) reverse map: for each NEIGHBOR centroid c',
            # khat_{i',c'} = index of the nearest LOCAL centroid
            # (solves component I).  Non-bijective: several neighbor
            # centroids may map to the same local centroid, which is
            # why |U_c| varies across c (contrast with Algorithm 2).
            rev = D.argmin(axis=0)
            fwd_list.append(fwd)
            rev_list.append(rev)
            # component II contributions: exactly one matched neighbor
            # centroid per local c per neighbor.
            num += alpha * mu_j[fwd]
            den += alpha
            # component I contributions: scatter each neighbor centroid
            # onto its nearest local centroid (the set U_c of Alg1-M).
            np.add.at(num, rev, alpha * mu_j)
            np.add.at(den, rev, alpha)
        # ---- (Alg1-M) closed-form centroid update --------------------
        # w~^(i,c) <- num / (den + eps).  A vanishing denominator means
        # cluster c received neither datapoints nor matched neighbor
        # centroids; per the convention of Proposition 1, such a
        # centroid retains its previous value (the `mask`).
        W_new = W.copy()
        mask = den > 0
        W_new[mask] = num[mask] / (den[mask, None] + eps)
        W = W_new
        if tau == 1:
            # Store the pass-1 assignments pi_1 and iterate W_1 for the
            # safeguard (Algorithm 1, line after the M step).
            pi1 = (assign, fwd_list, rev_list)
            W1 = W.copy()
        # ---- delta stopping rule -------------------------------------
        # Leave the loop early once a pass decreases the objective of
        # eq. (equ_def_Jacobi_GTVMin) by less than delta.
        obj = _local_obj_nc(X, W, neigh_mus, alpha)
        if prev_obj is not None and prev_obj - obj < delta:
            break
        prev_obj = obj
    # ---- Safeguard (regime (ii) of Proposition 1) --------------------
    # Accept the multi-pass iterate only if it does not increase the
    # pass-1 surrogate relative to the round's input:
    #     F(W ; pi_1) <= F(W_in ; pi_1),
    # else fall back to the pass-1 iterate W_1 (which satisfies the
    # inequality by construction, being the minimizer of F(.; pi_1)).
    # Locally checkable: needs only X, the neighbors' centroids, pi_1.
    assign, fwd_list, rev_list = pi1
    triggered = (_surrogate_F(W, X, neigh_mus, assign, fwd_list, rev_list, alpha)
                 > _surrogate_F(W_in, X, neigh_mus, assign, fwd_list, rev_list, alpha))
    if triggered:
        W = W1
    return W, triggered


def _local_obj_nc(X, W, neigh_mus, alpha):
    """The TRUE (assignment-minimized) objective of
    eq. (equ_def_Jacobi_GTVMin): local loss L_i(W) plus alpha times the
    ell=1 discrepancies to the frozen neighbor centroids.  Used only
    for the delta stopping rule -- distinct from the frozen-assignment
    surrogate F above.
    """
    D = ((X[:, None, :] - W[None, :, :]) ** 2).sum(axis=2)
    val = D.min(axis=1).sum() / max(len(X), 1)
    for mu_j in neigh_mus:
        val += alpha * gtv_discrepancy(W, mu_j, ell=1)
    return val


# ----------------------------------------------------------------------
# Algorithm 2: local solver via optimal transport (FedKM-OT)
# ----------------------------------------------------------------------
# Solves eq. (equ_def_Jacobi_GTVMin_perm): same structure as
# Algorithm 1, but the two nearest-centroid steps (E2)+(E3) are
# replaced by ONE bijective assignment per neighbor (the Hungarian
# permutation s_{i,i'}), so each local centroid is matched to exactly
# one neighbor centroid per neighbor and den is constant across c.

def _surrogate_G(W, X, neigh_mus, assign, perm_list, alpha):
    """G(W; theta_1), eq. (equ_upper_bound_G): the convex surrogate of
    eq. (equ_def_Jacobi_GTVMin_perm) with all assignment variables
    frozen at theta_1 = ({khat_r}, {s_{i,i'}}):

        (1/m_i) sum_r || W[khat_r] - x_r ||^2
        + alpha * sum_{i'} sum_c || W[c] - w^(i', s_{i,i'}(c)) ||^2 .
    """
    val = np.sum((W[assign] - X) ** 2) / len(X)
    for mu_j, perm in zip(neigh_mus, perm_list):
        val += alpha * np.sum((W - mu_j[perm]) ** 2)
    return val


def local_update_ot(X, mu, neigh_mus, alpha, P=10, delta=DELTA, eps=EPS):
    """Algorithm 2 of the paper.  Returns (W, safeguard_triggered).
    Same interface as local_update_nc."""
    K, d = mu.shape
    m = max(len(X), 1)
    W_in = mu.copy()
    W = mu.copy()
    theta1 = W1 = None
    prev_obj = None
    for tau in range(1, P + 1):
        # ---- (Alg2-E1) local assignment step (identical to Alg1-E1) --
        assign = pairwise_distances_argmin(X, W)
        # Accumulators for (Alg2-M):
        #   num[c] = (1/m_i) sum_{r: khat_r=c} x_r
        #            + alpha * sum_{i'} w^(i', s_{i,i'}(c))
        #   den[c] = |{r: khat_r=c}| / m_i + alpha * |N(i)|
        # (den's neighbor part is constant in c -- the bijectivity
        # advantage discussed after Algorithm 2 in the paper.)
        perm_list = []
        num = np.zeros((K, d))
        den = np.zeros(K)
        n_c = np.bincount(assign, minlength=K).astype(float)
        for c in range(K):
            if n_c[c] > 0:
                num[c] = X[assign == c].sum(axis=0)
        num /= m
        den += n_c / m
        for mu_j in neigh_mus:
            # ---- (Alg2-E2) bijective assignment via the Hungarian
            # method: s_{i,i'} = argmin over permutations s of
            # sum_c || w^(i,c) - w^(i',s(c)) ||^2.  This single step
            # replaces (Alg1-E2) + (Alg1-E3).  Cost K^3 per neighbor
            # vs 2K^2 for Algorithm 1 (paper, Sec. V-D runtime remark).
            D = ((W[:, None, :] - mu_j[None, :, :]) ** 2).sum(axis=2)
            row, col = linear_sum_assignment(D)
            perm = np.empty(K, dtype=int)
            perm[row] = col                      # perm[c] = s_{i,i'}(c)
            perm_list.append(perm)
            num += alpha * mu_j[perm]
            den += alpha
        # ---- (Alg2-M) closed-form centroid update (same empty-cluster
        # convention as Alg1-M) -----------------------------------------
        W_new = W.copy()
        mask = den > 0
        W_new[mask] = num[mask] / (den[mask, None] + eps)
        W = W_new
        if tau == 1:
            # Store theta_1 and W_1 for the safeguard.
            theta1 = (assign, perm_list)
            W1 = W.copy()
        # ---- delta stopping rule --------------------------------------
        obj = _local_obj_ot(X, W, neigh_mus, alpha)
        if prev_obj is not None and prev_obj - obj < delta:
            break
        prev_obj = obj
    # ---- Safeguard: G(W; theta_1) <= G(W_in; theta_1) or fall back ----
    assign, perm_list = theta1
    triggered = (_surrogate_G(W, X, neigh_mus, assign, perm_list, alpha)
                 > _surrogate_G(W_in, X, neigh_mus, assign, perm_list, alpha))
    if triggered:
        W = W1
    return W, triggered


def _local_obj_ot(X, W, neigh_mus, alpha):
    """True objective of eq. (equ_def_Jacobi_GTVMin_perm), for the
    delta stopping rule (ell=2 discrepancies re-minimized over s)."""
    D = ((X[:, None, :] - W[None, :, :]) ** 2).sum(axis=2)
    val = D.min(axis=1).sum() / max(len(X), 1)
    for mu_j in neigh_mus:
        val += alpha * gtv_discrepancy(W, mu_j, ell=2)
    return val


# ----------------------------------------------------------------------
# Algorithm 3: outer parallel Jacobi iteration
# ----------------------------------------------------------------------

def kmeanspp_init(node_data, K, seed):
    """Independent k-means++ per node (Algorithm 3, first step /
    Sec. V: "apply k-means++ initialization independently on each
    node").  Independence across nodes is deliberate: it creates the
    label-mismatch scenario that breaks label-wise couplings (Forero,
    DGC-KM) but not our permutation-invariant GTV measures.  n_init=10
    restarts, as in the notebooks.
    """
    return [KMeans(n_clusters=K, init="k-means++", n_init=10,
                   random_state=seed).fit(X).cluster_centers_
            for X in node_data]


def fed_kmeans(node_data, A, K, ell, alpha=1.0, T=300, P=10, gamma=1.0,
               seed=0, callback=None, log_objective=False):
    """Algorithm 3 with local solver Algorithm 1 (ell=1, FedKM-NC) or
    Algorithm 2 (ell=2, FedKM-OT).

    Parameters mirror the Require line of Algorithm 3: local datasets,
    graph (adjacency A), GTV selector ell, max iterations T, damping
    gamma in (0,1] with default 1 (plain Jacobi), plus alpha and P
    forwarded to the local solver.

    callback(t, Ws) is invoked once before the loop (t=-1, the
    k-means++ initialization) and after every outer iteration -- used
    by the experiments to record xi_1/xi_2 trajectories.

    Returns (Ws, info) with
      info['safeguard_count']  total number of local safeguard
                               fallbacks across all devices/rounds
                               (experiment spec E2 asks how often it
                               triggers; empirically: rarely), and
      info['objective']        f_ell after each iteration when
                               log_objective=True (descent check).
    """
    local = local_update_nc if ell == 1 else local_update_ot
    n = len(node_data)
    G = nx.from_numpy_array(A)
    Ws = kmeanspp_init(node_data, K, seed)
    safeguard_count = 0
    obj_hist = []
    if callback is not None:
        callback(-1, Ws)
    for t in range(T):
        # Snapshot = the Jacobi (parallel) schedule: every device
        # computes its candidate from the neighbors' centroids at the
        # START of the round.  (A Gauss-Seidel variant would read
        # already-updated neighbors -- not what Algorithm 3 does.)
        snapshot = [W.copy() for W in Ws]
        new_Ws = []
        for i in range(n):
            neigh_mus = [snapshot[j] for j in G.neighbors(i)]
            W_tilde, trig = local(node_data[i], snapshot[i], neigh_mus,
                                  alpha=alpha, P=P)
            safeguard_count += trig
            # eq. (equ_damped_update):
            #   W^(i) <- (1-gamma) W^(i) + gamma W~^(i);
            # gamma=1 is the plain Jacobi update used in Sec. V.
            new_Ws.append((1.0 - gamma) * snapshot[i] + gamma * W_tilde)
        Ws = new_Ws
        if log_objective:
            obj_hist.append(gtvmin_objective(node_data, Ws, A, alpha, ell))
        if callback is not None:
            callback(t, Ws)
    return Ws, {"safeguard_count": safeguard_count,
                "objective": np.array(obj_hist)}


# ----------------------------------------------------------------------
# Baselines
# ----------------------------------------------------------------------

def local_kmeans(node_data, K, seed):
    """Local k-means baseline (Sec. V-A): every node clusters its own
    local dataset independently, no communication.  Lower bound on
    collaboration; its xi_1 gap to the federated methods quantifies
    the benefit of the GTV coupling."""
    return [KMeans(n_clusters=K, init="k-means++", n_init=10,
                   random_state=seed).fit(X).cluster_centers_
            for X in node_data]


def central_kmeans(node_data, K, seed):
    """Central k-means baseline (Sec. V-A): standard k-means on the
    pooled data of all nodes -- the ideal reference W* used by xi_1."""
    X = np.vstack(node_data)
    return KMeans(n_clusters=K, init="k-means++", n_init=10,
                  random_state=seed).fit(X).cluster_centers_


def distributed_kmeans(node_data, A, K, eta=1.0, T=300, seed=0,
                       shared_init=False, callback=None):
    """Distributed k-means (Forero et al., JSTSP 2011): ADMM on the
    label-wise consensus formulation, eq. (eq:dist_kmeans) of the
    draft.  Implementation follows Xu's notebook, which follows the
    closed-form updates of their paper (their eqs. (11a)-(11c)):

      per round t, per node j:
        1. assign local datapoints to the nearest local centroid
           (binary memberships mu),
        2. centroid update (primal):
             w_j,c <- [ sum_{r in cluster c} x_r  -  2 lambda_j,c
                        + eta (|N_j| w_j,c + sum_{j' in N_j} w_j',c) ]
                      / [ |cluster c| + 2 eta |N_j| ]
      then, after ALL nodes updated (synchronous broadcast):
        3. dual update:
             lambda_j,c += (eta/2) (|N_j| w_j,c - sum_{j' in N_j} w_j',c)

    eta is the augmented-Lagrangian penalty of the ADMM solver -- it
    tunes convergence speed, NOT the consensus itself, which remains a
    hard constraint (Sec. V-A of the draft).  Crucially, centroids are
    coupled BY INDEX c: with independent per-node initializations, the
    constraints may tie together centroids of different clusters --
    the label-mismatch failure mode shown in Fig. centroid_traj.

    shared_init=True gives the 'Distributed k-means (shared)' variant
    of Sec. V-D: every node starts from node 0's k-means++ centroids,
    which removes label mismatch by construction (at the price of one
    extra coordination round, not counted here).
    """
    n = len(node_data)
    d = node_data[0].shape[1]
    G = nx.from_numpy_array(A)
    W = np.zeros((n, K, d))
    if shared_init:
        km = KMeans(n_clusters=K, init="k-means++", n_init=10,
                    random_state=seed).fit(node_data[0])
        W[:] = km.cluster_centers_
    else:
        # Independent inits; seed+j gives each node a different
        # k-means++ draw (matches the notebook convention).
        for j in range(n):
            km = KMeans(n_clusters=K, init="k-means++", n_init=10,
                        random_state=seed + j).fit(node_data[j])
            W[j] = km.cluster_centers_
    lam = np.zeros((n, K, d))          # dual variables lambda_j,c
    W_new = np.zeros_like(W)
    if callback is not None:
        callback(-1, list(W))
    for t in range(T):
        for j in range(n):
            Xj = node_data[j]
            assign = pairwise_distances_argmin(Xj, W[j])
            neigh = list(G.neighbors(j))
            for c in range(K):
                members = Xj[assign == c]
                s = members.sum(axis=0) if len(members) else np.zeros(d)
                nom = (s - 2 * lam[j, c]
                       + eta * (len(neigh) * W[j, c] + W[neigh, c].sum(axis=0)))
                den = len(members) + 2 * eta * len(neigh)
                W_new[j, c] = nom / den if den > 0 else W[j, c]
        # Synchronous broadcast: duals are updated only after every
        # node has computed its new centroids.
        W = W_new.copy()
        for j in range(n):
            neigh = list(G.neighbors(j))
            for c in range(K):
                lam[j, c] += eta / 2 * (len(neigh) * W[j, c]
                                        - W[neigh, c].sum(axis=0))
        if callback is not None:
            callback(t, list(W))
    return list(W)


# ----------------------------------------------------------------------
# Aggregation helpers
# ----------------------------------------------------------------------

def mean_ci(values, conf=0.95):
    """Mean and 95% CI half-width via the Student t-distribution --
    the convention of the notebooks and of every 'mean +/- 95% CI'
    entry in the paper's tables.  With fewer than two values the CI is
    reported as 0."""
    v = np.asarray(values, dtype=float)
    m = v.mean()
    if len(v) < 2:
        return m, 0.0
    half = v.std(ddof=1) / np.sqrt(len(v)) * t_dist.ppf((1 + conf) / 2,
                                                        len(v) - 1)
    return m, half


class MetricRecorder:
    """Callback for fed_kmeans / distributed_kmeans recording the full
    xi_1 / xi_2 trajectories (one value per outer iteration, including
    the initialization at t=-1)."""

    def __init__(self, mu_central, A):
        self.mu_central = mu_central
        self.A = A
        self.gcd, self.cv = [], []

    def __call__(self, t, Ws):
        self.gcd.append(gcd(self.mu_central, Ws))
        self.cv.append(cv(Ws, self.A))


# ----------------------------------------------------------------------
# Experiments
# ----------------------------------------------------------------------
# Parameter provenance: every constant below is stated either in
# Sec. V of the draft or in Xu's notebooks (where the paper is silent,
# e.g. the ER seed convention).  --quick shrinks T / sweep grids for a
# smoke test; it does NOT change the algorithms.

DATASETS = ["iso", "aniso", "var"]


# ----------------------------------------------------------------------

# DGC-KM_rho (Armacki et al., IEEE TSP 2025, Algorithm 1 + eq. (7))
# ----------------------------------------------------------------------
# Their framework minimizes the rho-relaxed objective (their eq. (4))
#
#     J_rho(x, C) = (1/rho) J(x, C) + (1/2) <x, L x>,
#
# where x stacks all users' center estimates, J is the (weighted)
# clustering loss and L is the graph Laplacian.  Note the consensus
# term couples centers BY INDEX k -- there is no matching step -- which
# is exactly the property our permutation-invariant GTV measures avoid.
#
# Their Theorem 1: for step size alpha < (beta/rho + lambda_max(L))^-1
# the center sequence converges to a fixed point (clusters optimal for
# the centers, centers stationary for J_rho); Theorem 2: as
# rho -> infty the fixed points approach consensus fixed points, which
# for the k-means loss are Lloyd points of the pooled data (Lemma 1).

def dgc_km(node_data, A, K, rho=10.0, B=3, T=200, seed=0, callback=None):
    """DGC-KM_rho with uniform sample weights w_{i,r} = 1/N (as in
    their Sec. V experiments) and independent per-node k-means++
    initialization (their Algorithm 1 allows arbitrary, unsynchronized
    initialization -- highlighted as a feature after Theorem 1).

    Per round t (their Algorithm 1):
      steps 2-5: every node assigns each local datapoint to the
                 nearest of its OWN centers (cluster update; the
                 distance g is Euclidean per their Assumption 6),
      steps 7-10: B gradient updates of the centers, their eq. (6)/(7),
                 with a neighbor exchange of centers before each --
                 so one round costs B communications (their Remark 8).

    The k-means center update, their eq. (7), reads per center k:

      x_i(k) <- (1 - alpha [ W_ik / rho + |N_i| ]) x_i(k)
                + (alpha/rho) sum_{r in C_i(k)} w_r y_r
                + alpha sum_{j in N_i} x_j(k)

    with W_ik = sum_{r in C_i(k)} w_r the total weight assigned to
    center k at node i.  The first line is the consensus part
    (gradient of the Laplacian term), the second the innovation part
    (gradient of the 1/rho-weighted loss).  Note the neighbor sum runs
    over the SAME index k at all nodes -- the label-wise coupling.

    Step size: their Theorem 1 requires
    alpha < (beta/rho + lambda_max(L))^{-1}; for the k-means loss with
    uniform weights, the smoothness constant is beta = max_i m_i / N.
    We take 90% of the bound.  This is the theoretically safe (and
    conservative) choice; no tuning is attempted.
    """
    n = len(node_data)
    d = node_data[0].shape[1]
    N = sum(len(X) for X in node_data)
    w = 1.0 / N                      # uniform sample weights, their Sec. V
    G = nx.from_numpy_array(A)
    L = nx.laplacian_matrix(G).toarray()
    lam_max = np.linalg.eigvalsh(L).max() if L.any() else 1.0
    beta = max(len(X) for X in node_data) / N
    alpha = 0.9 / (beta / rho + lam_max)   # 90% of the Thm. 1 bound

    # Independent per-node k-means++ (seed+i so each node draws its
    # own initialization -- the label-mismatch regime).
    X_centers = np.zeros((n, K, d))
    for i in range(n):
        km = KMeans(n_clusters=K, init="k-means++", n_init=10,
                    random_state=seed + i).fit(node_data[i])
        X_centers[i] = km.cluster_centers_

    if callback is not None:
        callback(-1, list(X_centers))
    for t in range(T):
        # ---- cluster step (their steps 2-5) --------------------------
        # Assignments are fixed for the whole round; precompute, per
        # node and per center k, the weighted sum of assigned points
        # (sums) and the total assigned weight W_ik (wsum) that enter
        # the innovation part of eq. (7).
        sums = np.zeros((n, K, d))
        wsum = np.zeros((n, K))
        for i in range(n):
            assign = pairwise_distances_argmin(node_data[i], X_centers[i])
            for c in range(K):
                members = node_data[i][assign == c]
                if len(members):
                    sums[i, c] = w * members.sum(axis=0)
                    wsum[i, c] = w * len(members)
        # ---- B center updates of eq. (7) (their steps 7-10) ----------
        # Snapshot before each update = synchronous neighbor exchange:
        # all nodes read their neighbors' centers from the start of
        # the inner step (parallel schedule, as in their Algorithm 1).
        for b in range(B):
            snapshot = X_centers.copy()
            for i in range(n):
                neigh = list(G.neighbors(i))
                deg = len(neigh)
                # sum_{j in N_i} x_j(k) for all k at once (K x d).
                nb_sum = snapshot[neigh].sum(axis=0) if deg else 0.0
                X_centers[i] = ((1.0 - alpha * (wsum[i, :, None] / rho + deg))
                                * snapshot[i]
                                + (alpha / rho) * sums[i]
                                + alpha * nb_sum)
        if callback is not None:
            callback(t, list(X_centers))
    return list(X_centers)


def gen_two_communities(n_per_comm, m, K, shift, seed):
    """Two communities of n_per_comm nodes; community A samples from a
    K-component isotropic mixture with centers C, community B from
    C + shift.  Returns (node_data, comm) with comm[i] in {0, 1}."""
    centers = generate_centers(K=K, dim=2, min_sep=5.0, seed=seed)
    node_data, comm = [], []
    rng = np.random.default_rng(seed)
    for c, offset in [(0, 0.0), (1, shift)]:
        C = centers + offset
        for _ in range(n_per_comm):
            X, _ = make_blobs(n_samples=m, centers=C, cluster_std=1.0,
                              random_state=int(rng.integers(0, 2**31 - 1)))
            node_data.append(X)
            comm.append(c)
    return node_data, np.array(comm)


def build_sbm_graph(comm, p_in, p_out, seed):
    """SBM adjacency; resampled until every community subgraph is
    connected (so the intra-community GTV coupling can spread)."""
    n = len(comm)
    rng = np.random.default_rng(seed)
    for _ in range(1000):
        A = np.zeros((n, n))
        for i in range(n):
            for j in range(i + 1, n):
                p = p_in if comm[i] == comm[j] else p_out
                if rng.random() < p:
                    A[i, j] = A[j, i] = 1.0
        ok = all(nx.is_connected(nx.from_numpy_array(
                 A[np.ix_(comm == c, comm == c)])) for c in (0, 1))
        if ok:
            return A
    raise RuntimeError("could not sample a per-community-connected SBM")


def gcd_communities(Ws, comm, refs):
    """xi_1^comm: per-node deviation from the node's community-wise
    central k-means solution, Hungarian-aligned (reuses gcd per
    community and averages with node weights)."""
    n = len(Ws)
    total = 0.0
    for c in (0, 1):
        idx = np.where(comm == c)[0]
        total += gcd(refs[c], [Ws[i] for i in idx]) * len(idx)
    return total / n

COMM_ALPHAS = [0.01, 0.1, 1.0]


def exp_communities(outdir, seeds, quick=False):
    T = 50 if quick else 200
    n_per, m, K, shift = 5, 50, 3, 15.0
    rows = []
    for p_out in [0.0, 0.1]:
        for seed in seeds:
            seed = int(seed)
            node_data, comm = gen_two_communities(n_per, m, K, shift, seed)
            A = build_sbm_graph(comm, p_in=0.8, p_out=p_out, seed=seed)
            n_edges_x = int(A[np.ix_(comm == 0, comm == 1)].sum())
            # community-wise and pooled central references
            refs = {c: central_kmeans([node_data[i]
                                       for i in np.where(comm == c)[0]],
                                      K, seed) for c in (0, 1)}
            mu_glob = central_kmeans(node_data, K, seed)

            results = {}
            results["Local $k$-means"] = local_kmeans(node_data, K, seed)
            for alpha in COMM_ALPHAS:
                for ell, meth in [(1, "FedKM-NC"), (2, "FedKM-OT")]:
                    Ws, _ = fed_kmeans(node_data, A, K, ell=ell,
                                       alpha=alpha, T=T, P=10, seed=seed)
                    results[f"{meth} ($\\regparam={alpha:g}$)"] = Ws
            results["Distributed $k$-means"] = distributed_kmeans(
                node_data, A, K, eta=1.0, T=T, seed=seed)
            results["$\\text{DGC-KM}_{1}$"] = dgc_km(
                node_data, A, K, rho=1.0, B=3, T=T, seed=seed)

            for meth, Ws in results.items():
                rows.append({
                    "p_out": p_out, "seed": seed, "method": meth,
                    "cross_edges": n_edges_x,
                    "gcd_comm": gcd_communities(Ws, comm, refs),
                    "gcd_glob": gcd(mu_glob, Ws),
                })
            print(f"  [communities] p_out={p_out} seed {seed} done "
                  f"({n_edges_x} cross edges)")
    df = pd.DataFrame(rows)
    df.to_csv(outdir / "communities.csv", index=False)
    return df


def summarize_communities(df):
    print("\n== Experiment C: two weakly linked communities ==")
    for (p_out, meth), g in df.groupby(["p_out", "method"], sort=False):
        mc, cc = mean_ci(g["gcd_comm"].to_numpy())
        mg, cg = mean_ci(g["gcd_glob"].to_numpy())
        print(f"  p_out={p_out}  {meth:22s}  "
              f"xi1_comm {mc:9.3e} +/- {cc:.2e}   "
              f"xi1_glob {mg:9.3e} +/- {cg:.2e}")


def write_communities_table(df, path):
    """LaTeX table body (single-column layout): rows = method,
    columns = xi1_comm for p_out in {0, 0.1}; xi1_glob is summarized
    in the text."""

    def fmt(vals):
        mu, ci = mean_ci(vals)
        if mu != 0 and (abs(mu) < 1e-2 or abs(mu) >= 1e3):
            e = int(np.floor(np.log10(abs(mu))))
            return (f"$({mu / 10**e:.2f}\\pm{ci / 10**e:.2f})"
                    f"{{\\times}}10^{{{e}}}$")
        return f"${mu:.2f}\\pm{ci:.2f}$"

    methods = ([f"{m} ($\\regparam={a:g}$)"
                for m in ("FedKM-NC", "FedKM-OT") for a in COMM_ALPHAS]
               + ["Distributed $k$-means", "$\\text{DGC-KM}_{1}$",
                  "Local $k$-means"])
    with open(path, "w") as f:
        f.write("% Auto-generated by NumExpXu/regime_and_communities.py"
                " (exp_communities).\n")
        f.write("\\begin{tabular}{l cc}\n\\toprule\n")
        f.write("& \\multicolumn{2}{c}{$\\GCD^{\\rm comm}$} \\\\\n")
        f.write("\\cmidrule(lr){2-3}\n")
        f.write("Method & $p_{\\rm out}=0$ & $p_{\\rm out}=0.1$ \\\\\n")
        f.write("\\midrule\n")
        for meth in methods:
            cells = []
            for p_out in [0.0, 0.1]:
                g = df[(df.p_out == p_out) & (df.method == meth)]
                cells.append(fmt(g["gcd_comm"].to_numpy()))
            f.write(meth + " & " + " & ".join(cells) + " \\\\\n")
        f.write("\\bottomrule\n\\end{tabular}\n")
    print(f"  wrote {path}")


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--outdir", default="outputs")
    ap.add_argument("--seeds", type=int, default=10,
                    help="number of seeds (max 10, paper uses 10)")
    ap.add_argument("--quick", action="store_true",
                    help="tiny settings for a smoke test")
    args = ap.parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    seeds = SEEDS[: (3 if args.quick else args.seeds)]
    df = exp_communities(outdir, seeds, quick=args.quick)
    summarize_communities(df)
    write_communities_table(df, outdir / "communities_table.tex")
    print("done; results in", outdir.resolve())


if __name__ == "__main__":
    main()
