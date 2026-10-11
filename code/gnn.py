"""Two-layer mean-aggregator GraphSAGE for node classification, in NumPy.

h1 = relu( X Ws0 + (A_norm X) Wn0 )
z  = h1 Ws1 + (A_norm h1) Wn1          (logit per node)

A_norm is the row-normalised adjacency (mean aggregation over neighbours).
Trained with Adam on masked binary cross-entropy.
"""
import numpy as np


def row_normalise(A):
    d = np.asarray(A.sum(axis=1)).ravel()
    d[d == 0] = 1.0
    return A / d[:, None]


def sigmoid(x):
    return 1.0 / (1.0 + np.exp(-np.clip(x, -40, 40)))


class GraphSAGE:
    def __init__(self, n_feat, n_hidden=32, seed=0):
        r = np.random.default_rng(seed)
        s0 = np.sqrt(2.0 / (n_feat + n_hidden))
        s1 = np.sqrt(2.0 / (n_hidden + 1))
        self.W = {"Ws0": r.normal(0, s0, (n_feat, n_hidden)),
                  "Wn0": r.normal(0, s0, (n_feat, n_hidden)),
                  "Ws1": r.normal(0, s1, (n_hidden, 1)),
                  "Wn1": r.normal(0, s1, (n_hidden, 1))}
        self.m = {k: np.zeros_like(v) for k, v in self.W.items()}
        self.v = {k: np.zeros_like(v) for k, v in self.W.items()}
        self.t = 0

    def forward(self, X, An):
        AX = An @ X
        pre1 = X @ self.W["Ws0"] + AX @ self.W["Wn0"]
        h1 = np.maximum(pre1, 0)
        Ah1 = An @ h1
        z = h1 @ self.W["Ws1"] + Ah1 @ self.W["Wn1"]
        return z, dict(AX=AX, pre1=pre1, h1=h1, Ah1=Ah1)

    def loss_grad(self, X, An, y, mask, l2=1e-4):
        z, c = self.forward(X, An)
        p = sigmoid(z).ravel()
        n = max(int(mask.sum()), 1)
        eps = 1e-9
        loss = -(y[mask] * np.log(p[mask] + eps) +
                 (1 - y[mask]) * np.log(1 - p[mask] + eps)).mean()
        dz = np.zeros_like(z)
        dz[mask, 0] = (p[mask] - y[mask]) / n
        g = {}
        g["Ws1"] = c["h1"].T @ dz
        g["Wn1"] = c["Ah1"].T @ dz
        dh1 = dz @ self.W["Ws1"].T + An.T @ (dz @ self.W["Wn1"].T)
        dpre1 = dh1 * (c["pre1"] > 0)
        g["Ws0"] = X.T @ dpre1
        g["Wn0"] = c["AX"].T @ dpre1
        for k in g:
            g[k] = g[k] + l2 * self.W[k]
            loss += 0.5 * l2 * float((self.W[k] ** 2).sum())
        return loss, g, p

    def step(self, g, lr=0.01, b1=0.9, b2=0.999, eps=1e-8):
        self.t += 1
        for k in self.W:
            self.m[k] = b1 * self.m[k] + (1 - b1) * g[k]
            self.v[k] = b2 * self.v[k] + (1 - b2) * g[k] ** 2
            mh = self.m[k] / (1 - b1 ** self.t)
            vh = self.v[k] / (1 - b2 ** self.t)
            self.W[k] -= lr * mh / (np.sqrt(vh) + eps)


def check_gradients(model, X, An, y, mask, eps=1e-6):
    """Finite-difference check; returns max relative error over sampled entries."""
    _, g, _ = model.loss_grad(X, An, y, mask)
    rng = np.random.default_rng(0)
    worst = 0.0
    for k in model.W:
        for _ in range(4):
            i = rng.integers(model.W[k].shape[0])
            j = rng.integers(model.W[k].shape[1])
            old = model.W[k][i, j]
            model.W[k][i, j] = old + eps
            lp, _, _ = model.loss_grad(X, An, y, mask)
            model.W[k][i, j] = old - eps
            lm, _, _ = model.loss_grad(X, An, y, mask)
            model.W[k][i, j] = old
            num = (lp - lm) / (2 * eps)
            den = max(abs(num), abs(g[k][i, j]), 1e-8)
            worst = max(worst, abs(num - g[k][i, j]) / den)
    return worst
