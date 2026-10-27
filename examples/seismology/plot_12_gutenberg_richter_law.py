"""
Gutenberg and Richter's frequency-magnitude law (1944)
======================================================

Beno Gutenberg and Charles Richter found that the number of earthquakes of at
least magnitude M falls as 10^(a - bM), with b close to 1: each magnitude unit
brings ten times fewer events. A synthetic catalogue drawn from that law
recovers b with Aki's maximum-likelihood estimator.
"""

import matplotlib.pyplot as plt
import numpy as np

from earthsciencekit.seismology import gutenberg_richter

rng = np.random.default_rng(1944)
b, m_min, n = 1.0, 2.0, 20000
magnitudes = m_min + rng.exponential(1 / (b * np.log(10)), n)
a = np.log10(n) + b * m_min
b_hat = np.log10(np.e) / (magnitudes.mean() - m_min)
assert abs(b_hat - b) < 0.03
assert np.isclose(gutenberg_richter(a, b, m_min), n)

m = np.linspace(m_min, magnitudes.max(), 200)
observed = (magnitudes[None, :] >= m[:, None]).sum(axis=1)
fig, ax = plt.subplots()
ax.semilogy(m, observed, ".", ms=3, label="synthetic catalogue")
ax.semilogy(m, gutenberg_richter(a, b, m), label=f"law with b = {b}")
ax.set(xlabel="magnitude M", ylabel="number with magnitude >= M", title=f"Estimated b = {b_hat:.3f}")
ax.legend()

plt.show()
