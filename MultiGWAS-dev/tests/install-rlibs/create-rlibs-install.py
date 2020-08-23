#!/usr/bin/python

import sys

args = sys.argv
REPOS="http://cran.r-project.org"

infile = args [1]
libs = open (infile).readlines ()
for l in libs:
    cmm = 'install.packages ("%s", dependencies=T, repos="%s")' % (l.strip(), REPOS)
    print (cmm)
