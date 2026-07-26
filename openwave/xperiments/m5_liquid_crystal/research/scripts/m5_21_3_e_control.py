"""M5.21.3 grind control: static (omega = 0) continuation of the P1
endpoint for the SAME cumulative depth as the omega-ladder (5 x 3000),
so the ladder's E_stat drift separates into grind vs omega-driven.

Run: python3 m5_21_3_e_control.py <tagbase> [maxit]
Out: ../data/m5_21_3_row_ctrl_<tagbase>.json
"""
import json
import os
import runpy
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
g4 = runpy.run_path(os.path.join(HERE, "m5_21_3_a_4d.py"),
                    run_name="not_main")

tagbase = sys.argv[1]
maxit = int(sys.argv[2]) if len(sys.argv) > 2 else 15000
M, cfg = g4["load_p1"](tagbase)
free = ~g4["pin_shell"](cfg["n"], cfg["h"])
M2, info = g4["fire"](M, cfg, free, maxit, tag=f"ctrl_{tagbase}")
eu, ev = g4["e_parts"](M2, cfg)
row = {"tag": f"ctrl_{tagbase}", "maxit": maxit,
       "E_end": float(eu + ev), "E_u": float(eu), "E_v": float(ev),
       "stop": info["stop"], "trace": info["trace"][-5:]}
np.savez_compressed(os.path.join(DATA, f"m5_21_3_ctrl_{tagbase}.npz"),
                    M=M2.astype(np.float32), s=cfg["s"],
                    delta=cfg["delta"], h=cfg["h"])
with open(os.path.join(DATA, f"m5_21_3_row_ctrl_{tagbase}.json"),
          "w") as f:
    json.dump(row, f, indent=1)
print(json.dumps({k: row[k] for k in ("tag", "E_end", "E_u", "E_v",
                                      "stop")}))
