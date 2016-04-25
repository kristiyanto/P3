# THIS SCRIPT WILL SCAN FOR THE FILES IN A FOLDER
# EXTRACT, IF NECESSARY, AND RUN MGSF PLUS TO 
# NORMALIZE THE DATA
# DANIEL.KRISTIYANTO@PNNL.GOV
# Docker container: msgf

######################## ENV ###########################
import gzip
import os
import glob
import os.path
import subprocess
import re
import csv
from ftplib import FTP
from _functions import *
from datetime import timedelta
from datetime import datetime
import time

working_dir = os.getcwd() +"/data/"

spectrum 	= []
out 		= []
db 			= None
module 		= None
input_csv	= None
start_time = time.time()
print("Start at:", str(datetime.now()))

db, input_csv 			= scan_dir(working_dir)
######################## DOWNLOAD FTP ###########################
if (input_csv != None):
	with open(input_csv) as f:
	    reader = csv.reader(f)
	    header = next(reader)
	    for row in reader:
	        for src in row:
	        	print(src)
	        	if (check_url(src)):
	        		get_ftp(src)

######################## DECOMPRESS ###########################
spectrum_tmp = []
for src_name in glob.glob(os.path.join(working_dir, '*.gz')):
    base = os.path.basename(src_name)
    #print("Extracting", src_name)
    dest_name = os.path.join(working_dir, base[:-3])
    if not os.path.isfile(dest_name):
    	spectrum_tmp.append(dest_name)
    	with gzip.open(src_name, 'rb') as infile:
	        with open(dest_name, 'wb') as outfile:
	            for line in infile:
	                outfile.write(line)

spectrum 				= scan_spectrum(working_dir)
db, input_csv 			= scan_dir(working_dir)
#module 		= dir_content['module']


######################## CHECK FILES ###########################
spectrum 				= scan_spectrum(working_dir)
db, input_csv 			= scan_dir(working_dir)

for s in spectrum:
	out = os.path.splitext(s)[0]+'.mzid'
	msgf(s,db,out)
	print("Done.")

elapsed_time = time.time() - start_time
print("Finished at:", str(datetime.now()))
print("Elapsed time:", str(timedelta(seconds=elapsed_time)))
