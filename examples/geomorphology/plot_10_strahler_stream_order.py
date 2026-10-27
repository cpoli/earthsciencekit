"""
Strahler stream order (1952)
============================

Arthur Strahler simplified Horton's ordering of streams: every unbranched
headwater channel is first order, and the order rises by one only where
two streams of equal order meet; a smaller tributary leaves the order
unchanged. Unlike Horton's scheme it needs no decision about which branch
continues upstream, so it can be computed from a network's links alone.
The example orders a random synthetic network and draws each link with a
width set by its order.
"""

import matplotlib.pyplot as plt
import numpy as np

from earthsciencekit.geomorphology import strahler_order

rng = np.random.default_rng(7)
downstream, depth = [-1], [0]
frontier = [0]
while frontier:
    link = frontier.pop()
    # Branch with a probability that falls with depth, so the tree stays finite.
    if depth[link] < 9 and rng.random() < 0.95 - 0.06 * depth[link]:
        for _ in range(2):
            downstream.append(link)
            depth.append(depth[link] + 1)
            frontier.append(len(downstream) - 1)
downstream = np.array(downstream)
order = strahler_order(downstream)
assert np.all(order[downstream[downstream >= 0]] >= order[downstream >= 0]) and order[0] >= 4

# Lay the tree out: sources spread along x, each link centred over its tributaries.
children = {i: list(np.flatnonzero(downstream == i)) for i in range(downstream.size)}
x = np.zeros(downstream.size)
counter = [0]


def place(link):
    if not children[link]:
        x[link] = counter[0]
        counter[0] += 1
    else:
        for child in children[link]:
            place(child)
        x[link] = np.mean([x[c] for c in children[link]])


place(0)
fig, ax = plt.subplots(figsize=(10, 5))
colors = plt.get_cmap("viridis")(np.linspace(0, 0.9, int(order.max())))
for link in range(downstream.size):
    top = -depth[link]
    k = int(order[link]) - 1
    ax.plot([x[link], x[link]], [top, top + 1 if link else top + 0.4], color=colors[k], lw=1 + 1.5 * k)
    if children[link]:
        span = [x[c] for c in children[link]]
        ax.plot([min(span), max(span)], [top, top], color=colors[k], lw=1 + 1.5 * k)
for k in range(int(order.max())):
    ax.plot([], [], color=colors[k], lw=1 + 1.5 * k, label=f"order {k + 1}: {np.count_nonzero(order == k + 1)} links")
ax.set(xticks=[], yticks=[], title=f"Strahler: a random network of order {order[0]:.0f}, outlet at the bottom")
ax.invert_yaxis()
ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1))
fig.tight_layout()

plt.show()
