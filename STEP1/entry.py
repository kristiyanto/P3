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
download_ftp(input_csv)

######################## DECOMPRESS ###########################
unzip_gz(working_dir)

spectrum 				= scan_spectrum(working_dir)
db, input_csv 			= scan_dir(working_dir)

######################## CHECK FILES ###########################
spectrum 				= scan_spectrum(working_dir)
db, input_csv 			= scan_dir(working_dir)

for s in spectrum:
	out = os.path.splitext(s)[0]+'.mzid'
	msgf(s,db,out)

elapsed_time = time.time() - start_time
print("Finished at:", str(datetime.now()))
print("Elapsed time:", str(timedelta(seconds=elapsed_time)))
