"""Cross-check only (not used to form any estimate): closed-form spectrum
lambda_{k,m} = (m + a_k)(m + a_k + 1),  a_k = k pi / (2W),  k >= 1, m >= 0.
Derivation in STAGE1.md (section 'Closed-form cross-check')."""
import math

from common import WIDTHS, NEIG, save_json


def closed(W, n=NEIG):
    vals = []
    kmax = 1
    while (kmax * math.pi / (2 * W)) ** 2 < 10 * (n + 2) ** 2 + 400:
        kmax += 1
    for k in range(1, kmax + 1):
        a = k * math.pi / (2 * W)
        for m in range(n + 1):
            vals.append(((m + a) * (m + a + 1), k, m))
    vals.sort()
    return vals[:n]


def main():
    out = {}
    for label, W in WIDTHS:
        c = closed(W)
        out[label] = {"W": W, "lowest": [v[0] for v in c], "labels_k_m": [[v[1], v[2]] for v in c]}
        print("C W=%-5s %s  labels=%s" % (label, " ".join("%.14f" % v[0] for v in c),
                                           [[v[1], v[2]] for v in c]))
    save_json("results_C.json", out)


if __name__ == "__main__":
    main()
