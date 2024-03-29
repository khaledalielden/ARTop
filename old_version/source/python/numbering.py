#!/usr/bin/env python
# coding: utf-8

from __future__ import division, print_function
import os
import sys
import datetime


def number_list(stime, etime, output_dir):
    expct_filesNo = (etime - stime).total_seconds()/720
#    expct_filesNo += 1 # start counting from zero
    add = open(output_dir+'/specifications.txt', 'a')
    add.write(str(int(expct_filesNo))+'\n')
    add.close()
    return expct_filesNo

sYear = sys.argv[1]
sMonth = sys.argv[2]
sDay = sys.argv[3]
sHour = sys.argv[4]

eYear = sys.argv[5]
eMonth = sys.argv[6]
eDay = sys.argv[7]
eHour = sys.argv[8]

outputdir = sys.argv[9]

start_time = datetime.datetime(int(sYear),int(sMonth),int(sDay),int(sHour),int(00),int(00))
end_time = datetime.datetime(int(eYear),int(eMonth),int(eDay),int(eHour),int(00),int(00))

nof = number_list(start_time, end_time, outputdir)

print('number of files = ', int(nof))

print('time span = ', end_time - start_time)


savefile= sys.argv[11]+"/savedpaths.dat"
if os.path.isfile(savefile):
    savepaths=open(savefile,'a')
else:
    savepaths=open(savefile,'w')
note={sys.argv[12]:[sys.argv[10],outputdir]}
for AR in note.keys():
    savepaths.write(str(AR)+':'+str(note[AR])+"\n")
savepaths.close()
