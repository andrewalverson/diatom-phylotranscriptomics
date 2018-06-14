#! /share/apps/Python/3.5.2/bin/python3.5

import sys
import os
import re

# get directory contents
folders = os.listdir(".")

for dir in folders:
    ids = dir.split("_")
    print("ln -s " + dir + "/*_1.fq.gz " + ids[0] + "_1.fq.gz")
    print("ln -s " + dir + "/*_2.fq.gz " + ids[0] + "_2.fq.gz")

