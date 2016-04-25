# THIS SCRIPT WILL SCAN FOR THE FILES IN A FOLDER
# EXTRACT, IF NECESSARY, AND RUN MGSF PLUS TO 
# NORMALIZE THE DATA
# DANIEL.KRISTIYANTO@PNNL.GOV
# Docker Container: itraquant

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
import time

working_dir = os.getcwd() +"/data/"

spectrum 	= []
out 		= []
db 			= None
module 		= None
input_csv	= None
start_time = time.time()

######################## READ CONFIG FILES ###########################

working_dir = os.getcwd() +"/data/"
#working_dir = "/Users/Daniel/Desktop/t/"
cfg_file = os.path.join(working_dir, "p3.config")
if os.path.isfile(cfg_file):
	r_options = get_itraq_opts(cfg_file)
else:
	write_config(working_dir)
	r_options = get_itraq_opts(cfg_file)

######################## DECOMPRESS ###########################
spectrum_tmp = []
for src_name in glob.glob(os.path.join(working_dir, '*.gz')):
    base = os.path.basename(src_name)
    # print("Extracting", src_name)
    dest_name = os.path.join(working_dir, base[:-3])
    if not os.path.isfile(dest_name):
    	spectrum_tmp.append(dest_name)
    	with gzip.open(src_name, 'rb') as infile:
	        with open(dest_name, 'wb') as outfile:
	            for line in infile:
	                outfile.write(line)
	                
spectrum 				= scan_spectrum(working_dir)
db, input_csv 			= scan_dir(working_dir)

######################## DOWNLOAD FTP ###########################
if (input_csv != None):
	with open(input_csv, 'rb') as f:
	    reader = csv.reader(f)
	    header = reader.next()
	    for row in reader:
	        for src in row:
	        	print(src)
	        	if (check_url(src)):
	        		get_ftp(src)

######################## CHECK FILES ###########################
spectrum 				= scan_spectrum(working_dir)
db, input_csv 			= scan_dir(working_dir)

if len(scan_mzid(working_dir)) != 0:
	try:
		print("Filtering.")
		rscript(r_options)
		print("Done.")
	except:
		print("Filtering failed.")
else:
	print("Missing MZID Files.")

elapsed_time = time.time() - start_time
print("Elapsed time:", str(timedelta(seconds=elapsed_time)))
