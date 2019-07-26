#!/usr/bin/env python

import os
import sys
import subprocess

os.system("rm -rf src/resources/*")

def xxd(src, dst):
    os.system("xxd -i " + src + " " + " " + dst)

data  = [ "data/manuals/anaquin.txt",
          "data/manuals/rna.txt",
          "data/manuals/meta.txt",
          "data/manuals/split.txt",
          "data/manuals/somatic.txt",
          "data/manuals/germline.txt",
          "data/manuals/calibrate.txt",

          "scripts/report.py",
          "scripts/template/calibrate.html",          
          "scripts/template/gSplit.html",
          "scripts/template/rSplit.html",
          "scripts/template/mSplit.html",
          "scripts/template/somatic.html",          
          "scripts/template/germline.html",
          "scripts/template/style.css",

          "src/r/plotGene.R",
          "src/r/plotInsert.R",
          "src/r/plotLinear.R",
          "src/r/plotAllele.R",
          "src/r/plotKSomatic.R",
          "src/r/plotWhisker.R",
          "src/r/plotLDensity.R",
    	  "src/r/plotLogistic.R",
          "src/r/plotQualFilter.R",
          "src/r/plotStrelkaROC.R",
          "src/r/plotKSynthetic.R"
        ]
tests = [ ]

for file in os.listdir("data/manuals/"):
    path = os.path.join("data/manuals/", file)
    if os.path.isfile(path):
        os.system("expand -t 4 " + path + " > /tmp/tmp.txt")
        os.system("mv /tmp/tmp.txt " + path)

def checkEOF(file):
    with open(file) as r:
        x = r.read()
        if "<<@@@@>>" not in x:
            x = x + "<<@@@@>>"
            #raise Exception("No <<@@@@>> in " + file)
    with open(file, "w") as w:
        w.write(x)

for i in range(0, len(data)):
    file = os.path.basename(data[i])
    dst = "src/resources/" + file
    checkEOF(data[i])    
    xxd(data[i], dst)

os.system("rm src/data/resources.o")
