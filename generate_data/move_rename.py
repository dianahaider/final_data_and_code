#!/usr/bin/env python3
import shutil
import os
import sys

dr = sys.argv[1]

for root, dirs, files in os.walk(dr):
    for file in files:
        if file == "taxonomy.tsv":
            spl = root.split("/"); newname = spl[-5]; sup = ("/").join(spl[:-5])
            shutil.move(root+"/"+file, sup+"/"+newname+".tsv"); shutil.rmtree(root)
