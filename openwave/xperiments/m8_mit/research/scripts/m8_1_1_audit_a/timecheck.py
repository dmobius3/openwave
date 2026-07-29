import time
from driver import analyse
from part1 import build_groups
gd={n:(g,o) for n,g,o in build_groups()}
for nm in ("C_20","BD_12","2O","2I"):
    t=time.time(); analyse(nm,*gd[nm]); print("   -> %.1fs"%(time.time()-t))
