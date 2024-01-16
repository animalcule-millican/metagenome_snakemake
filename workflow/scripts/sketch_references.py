#!/usr/bin/env python3
import os
import sys
import subprocess 

dir_path = sys.argv[1]
output = sys.argv[2]
sketchout = sys.argv[3]
cpus = sys.argv[4]

magfiles = os.listdir(dir_path)
magnames = [x.replace(".fna.gz", "") for x in magfiles]
magpaths = [os.path.join(dir_path, x) for x in magfiles]
with open(output, "w") as f:
    f.write("name,genome_filename,protein_filename\n")
    for i in range(len(magfiles)):
        f.write(magnames[i] + "," + magpaths[i] + "," + "\n")

subprocess.run(["sourmash", "scripts", "manysketch", "-o", sketchout, "-c", cpus, "-p", "k=21,k=31,k=51,scaled=1000,abund", output])