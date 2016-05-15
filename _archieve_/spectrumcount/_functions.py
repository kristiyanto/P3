# DANIEL.KRISTIYANTO@PNNL.GOV

import os
import sys
import subprocess
import re
from os.path import isfile, join
from urlparse import urlparse
from ftplib import FTP

######################## SUB PROCESSES ###########################

def msgf(spectrum, db, out):
	#print(spectrum, db, out)
	if not os.path.isfile(out):
		try:
			subprocess.call(['java', '-Xmx3500M', '-jar', 'MSGFPlus.jar', '-s', spectrum, '-d', db, '-o', out])
		except:
			print("MZID conversion failed.")

def rscript():
	print("Filtering.")
	runr = subprocess.call(['Rscript','labelfreequant.R'])

######################## SCAN FILES ###########################

def scan_spectrum(working_dir):
	spectrum = []
	for file in os.listdir(working_dir):
		if file.lower().endswith((".mzml",".mgf", ".mzxml", ".ms2", ".pkl")):
			spectrum.append(join(working_dir, file))
			if (len(spectrum)==0):
				print("Missing spectrum files.")
				exit()
			#print("Spectrum: ", file)
	return spectrum

def scan_mzid(working_dir):
	mzid = []
	for file in os.listdir(working_dir):
		if file.lower().endswith((".mzid")):
			mzid.append(join(working_dir, file))
			if (len(mzid)==0):
				print("No MZID file found.")
				exit()
			#print("MZID: ", file)
	return mzid


def scan_dir(working_dir):
	db 			= None
	module 		= None
	input_csv	= None
	for file in os.listdir(working_dir):
		if file.lower().endswith(".fasta"):
			db 			= join(working_dir, file)
			if (db == None):
				print("Missing database file.")
				exit()
			#print("Database: ", file)
		if file.lower().endswith(".mod"):
			module		= join(working_dir, file)
			if (module == None):
				print("Missing module file.")
			#print("Module:", file)
		if file.lower().endswith(".csv"):
			input_csv 	= join(working_dir, file)
			if (input_csv != None):
				print(input_csv)
			#print("CSV:", file)
	return (db, input_csv)

######################## FTP HARVESTING ###########################


def check_url(url):
	regex = re.compile(
        r'^(?:ftp)s?://' 
        r'(?:(?:[A-Z0-9](?:[A-Z0-9-]{0,61}[A-Z0-9])?\.)+(?:[A-Z]{2,6}\.?|[A-Z0-9-]{2,}\.?)|' 
        r'localhost|' 
        r'\d{1,3}\.\d{1,3}\.\d{1,3}\.\d{1,3})' 
        r'(?::\d+)?' 
        r'(?:/?|[/?]\S+)$', re.IGNORECASE)
	return regex

def get_ftp(ftp_url):
	URL = urlparse(ftp_url)
	url_host	= URL.netloc
	url_path 	= URL.path
	print(url_path)
	print(url_host)
	work_dir 	= os.getcwd()
	os.chdir(work_dir+"/data")
	try:
		ftp 	= FTP(url_host)
		ftp.login()
		ftp.cwd(url_path)
		ftp.retrlines('LIST')
		filenames = ftp.nlst()
		print(filenames)
		filematch = "*.*" #Any files

		for filename in ftp.nlst(filematch):
		    fhandle = open(filename, 'wb')
		    print('Getting ' + filename)
		    ftp.retrbinary('RETR ' + filename, fhandle.write)
		    fhandle.close()
		print("FTP Download done.")
	except Exception:
		pass
	os.chdir(work_dir)
	return URL

