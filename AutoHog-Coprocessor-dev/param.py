import os
import sys


BASE_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(BASE_DIR, "lattice-estimator"))

from estimator import *


n = 711
q = 2**32
alpha_torus = 2.0**-17
sigma_zq = alpha_torus * q

Xs = ND.Binary
Xe = ND.DiscreteGaussian(stddev=sigma_zq)

def eval_for_m(m):
    params = LWE.Parameters(
        n=n,
        q=q,
        Xs=Xs,
        Xe=Xe,
        m=m,
        tag=f"lvl0_n{n}_m{m}",
    )
    print("=" * 70)
    print(f"Evaluating: n={n}, q=2^32, alpha=2^-17, m={m}")
    print("- rough()")
    LWE.estimate.rough(params)
    print("- full estimate()")
    LWE.estimate(params)

if __name__ == "__main__":
    for m in [n, 2 * n, 4 * n, 2**20]:
        eval_for_m(m)