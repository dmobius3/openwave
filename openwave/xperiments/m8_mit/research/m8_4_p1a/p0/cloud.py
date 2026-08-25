"""Orbit cloud on S³/2I for the RBF-FD quotient operator.

Generates seed points on S³, expands each under the 120-element 2I action,
and returns the full cloud with orbit and group-element indices.
"""

import numpy as np
from .algebra import qmul


def fibonacci_seeds_s3(n_seeds, rng_seed=0):
    """Generate well-separated seed points on S³ via generalized Fibonacci.

    Uses the 4D Fibonacci lattice approach: golden-ratio spacing in three
    angular coordinates of the Hopf-like parameterization.
    """
    phi = (1 + 5**0.5) / 2
    seeds = []
    for i in range(n_seeds):
        t = (i + 0.5) / n_seeds
        theta1 = np.arccos(1 - 2 * t)
        theta2 = 2 * np.pi * (i * phi) % (2 * np.pi)
        theta3 = 2 * np.pi * (i * phi**2) % (2 * np.pi)
        w = np.cos(theta1 / 2) * np.cos(theta2 / 2)
        x = np.cos(theta1 / 2) * np.sin(theta2 / 2)
        y = np.sin(theta1 / 2) * np.cos(theta3 / 2)
        z = np.sin(theta1 / 2) * np.sin(theta3 / 2)
        norm = (w**2 + x**2 + y**2 + z**2)**0.5
        seeds.append(np.array([w, x, y, z]) / norm)
    return seeds


def build_orbit_cloud(seeds, elems):
    """Expand seeds under 2I left-action to build the full orbit cloud.

    For each seed x and group element γ, the point is γ·x = qmul(γ, x).

    Returns:
        X: (N, 4) array of all cloud points
        oid: (N,) orbit index for each point
        gid: (N,) group-element index for each point
    """
    G = len(elems)
    N = len(seeds) * G
    X = np.zeros((N, 4))
    oid = np.zeros(N, dtype=int)
    gid = np.zeros(N, dtype=int)

    for si, seed in enumerate(seeds):
        for gi, gamma in enumerate(elems):
            idx = si * G + gi
            X[idx] = qmul(gamma, seed)
            oid[idx] = si
            gid[idx] = gi

    return X, oid, gid
