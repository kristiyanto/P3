##########################################################################
## P3: Portable Proteomics Pipeline
## DANIEL KRISTIYANTO (daniel.kristiyanto@pnnl.gov)
## APR 27, 2016
##########################################################################


import gzip
import sys
import os
import glob
import os.path
import subprocess 
import re
import csv
from ftplib import FTP
import configparser
from datetime import datetime
import time
from urllib.parse import urlparse

# Set Env
working_dir = "/root/data"
#working_dir = "/Users/Daniel/Desktop/LABELLED"
config_file = "p3.config"
maxwait = 3 * 60 * 60 # Wait for 3 hours
os.chdir(working_dir) if os.path.isdir(working_dir) else sys.exit("{} not found. Please make sure it's properly mounted.".format(working_dir))

def main():

	start_time = time.time()
	print("Starting at: {}".format(str(datetime.now())))
	
	# Scan for config file
	if os.path.isfile(config_file) is False:
		print("p3.config is not found.")
		write_blank_p3(config_file)

	# Gather the files
	options = check_config(config_file)
	get_files(options)
	
		
	# Unzip the files 
	unzip_files()

	# Re-List the files
	fasta = scan_fasta()
	sp_files_ext = (".mzml",".mgf",".mzxml",".ms2",".pkl", ".mzXML")
	spectrum = scan_files(sp_files_ext)
	
	if (len(spectrum)!=0):
		print("{} Spectrum files found: {}".format(len(spectrum), spectrum))
	else:
		print("No spectrum files found. Accepted format:{} ".format(sp_files_ext))
	
	Q_METHOD = options.get("QUANTIFICATION","METHOD").upper()
	RUN_MSGF = (options.get("MSGF", "RUN_MSGF").upper() == "YES")
	

	# Run MSGF 
	if fasta and spectrum and RUN_MSGF:
		print("Running MSGF") 

		try:
			msgf_opts = options.get("MSGF", "MSGF_OPTIONS").upper()
			if (len_msgf_opts !=0): msgf(spectrum, fasta, options, msgf_opts)
		except:
			msgf(spectrum, fasta, options)
	else:
		print("Skipping MSGF identification...")

	

	if Q_METHOD == "SPECTRUM_COUNT":
		keep_waiting = True
		wait = 0
		lock = "SC_RESULT.txt.tmp"
		if os.path.isfile(lock): sys.exit("Spectrum Count Quantification is already run by other container. Quitting.")
		while keep_waiting:
			mzid_set = scan_files((".mzid"))
			if len(mzid_set) == 0:
				sys.exit("No .mzid file found")
			elif len(scan_files('.mzid.tmp')) != 0: 
				print("SCAN MZID TMP:" + scan_files('.mzid.tmp'))
				print("MSGF identification is being run by other containers. Waiting...")
				if wait > maxwait:
					keep_waiting = False
					sys.exit("Had waited for 2 hours. Quitting now. \n Perhaps remove .tmp files and rerun containers?")
				else: wait = wait + 5 
				time.sleep(5)
			else:
				touch(lock)
				scquant(mzid_set, options)
				os.remove(lock)
				keep_waiting = False

	if Q_METHOD == "ITRAQ4":          
		mzid = scan_files(".mzid")
		for m in mzid: itraq(m, options)
		lock_files = scan_files(".rda.tmp")
		wait = 0
		while len(lock_files) != 0:
			print("Identification is being run by other containers. Waiting...")
			if wait > maxwait:
				keep_waiting = False
				sys.exit("Had waited for 2 hours. Quitting now. \n Perhaps remove .tmp files and rerun containers?")
			else: 
				wait = wait + 10
				lock_files = scan_files(".rda.tmp")
				time.sleep(10)
		
		itraq_folding(options)

	stop_time = time.time()
	print("Finished at: {}".format(str(datetime.now())))
		
	
################################################ FUNCTIONS  ################################################

def msgf(spectrum, fasta, options, *opt):
	is_itraq = options.get("QUANTIFICATION","METHOD").upper() == "ITRAQ4"
	for file in spectrum:
		out = os.path.splitext(file)[0] + ".mzid"
		lock = out + ".tmp"
		print("(MSGF) Working on: {}".format(file))

		if os.path.isfile(out):
			print("{} already exists.".format(out))
		elif os.path.isfile(lock):
			print("Other container is working on {}. Next...".format(lock))
		else:
			try:
				cmd = ['java', '-Xmx3500M', '-jar', '/root/MSGFPlus.jar', '-s', file, '-d', fasta, '-o', out]
				print(cmd)
				if len(opt) != 0: cmd.extend([str(x) for x in opt])
				touch(lock)
				subprocess.call(cmd)
				print(str(cmd))
				os.remove(lock)
			except:
				print("MSGF Error: {}".format(file))
			if is_itraq: itraq(out, options)


def itraq(mzid, options):
	ext_set = (".mzXML", ".mzML", ".MZXML", ".mzml", ".MZML", ".mzml")
	opts_set = ("EVALUE_TRESHOLD", "pNA", "QUANTIFICATION_METHOD")
	sp_files = scan_files(ext_set)
	lock = mzid[:-5] + ".rda.tmp"
	out = mzid[:-5]+".rda"
	if os.path.isfile(out) or os.path.isfile(lock):
		print("{} ... next.".format(lock))
	else:
		cmd = ['Rscript', "/root/itraq.R"]
		mzml = ""
		for ext in ext_set:
			m = mzid[:-5] + ext
			if m in sp_files:
				mzml = m
				cmd.append(mzml)
				cmd.append(mzid)
				opts = get_itraq_opts(options)
				for k in opts_set:
					cmd.append(str(opts[k]))
				break
		cmd.append(lock[:-4])
		print("(ITRAQ) Working on {}".format(ext_set))
		if len(mzml)!=0: 
			print("Quantifying:" + mzml)
			try:
				touch(lock)
				print(cmd)
				subprocess.call(cmd)
				os.remove(lock)
			except:
				print("iTRAQ4 quantification failed.")
   

def itraq_folding(options):
	try: c_by = options.get("ITRAQ4","COMBINE_BY").lower()
	except: c_by = mean
	finally:
		cmd = ['Rscript', "fold.R", c_by]

def get_pride(prideID):
	cmd = ['Rscript', "fold.R", c_by]
	
def get_itraq_opts(options):
	R_OPTS = dict(EVALUE_TRESHOLD= 1, pNA= 0, QUANTIFICATION_METHOD="count")
	for k,v in R_OPTS.items():
		try:
			x = options.get("ITRAQ4", k)
			if len(x) != 0 : R_OPTS[k] = x
		except:
			continue
	print("ITRAQ. {} : {}".format(k, R_OPTS[k]))
	return R_OPTS             

			
def scquant(mzid_set, options):
	R_OPTS = dict(SCORE_TRESHOLD= 7.0, ERROR_TRESHOLD= 20, FDR=0.001, ITERATION=5000)
	opts_set = ("SCORE_TRESHOLD", "ERROR_TRESHOLD", "FDR", "ITERATION")
	for k,v in R_OPTS.items():
		try: 
			x = options.get("SPECTRUM_COUNT", k)
			if len(x) != 0 : R_OPTS[k] = x
		except:
			continue
	cmd = ['Rscript', "/root/scquant.R"]
	for k in opts_set:
		cmd.append(str(R_OPTS[k]))
	print(cmd)
	subprocess.call(cmd)
	
	
   
			  
################################################ ADMINISTRIVIA  ################################################
	
def check_config(cfg):
	options = configparser.ConfigParser()
	try:
		options.read(cfg)
		global_o = options.get("SOURCE", "REPO").upper() in ("LOCAL", "FTP", "PRIDE") 
		msgf_o = options.get("MSGF", "RUN_MSGF").upper() in ("YES", "NO")
		quant_o = options.get("QUANTIFICATION", "METHOD").upper() in ("SPECTRUM_COUNT", "ITRAQ4", "NONE")  
		if  global_o and msgf_o and quant_o:
			return options
		else:
			sys.exit()
	except:
		if not global_o: sys.exit("Configuration error: SOURCE. Check p3.config.")
		if not msgf_o: sys.exit("Configuration error: MSGF. Check p3.config.")
		if not quant_o: sys.exit("Configuration error: QUANTIFICATION. Check p3.config.")
				
def get_files(options):
	try: 
		src = options.get("SOURCE", "REPO").upper()
	except:
		sys.exit("Configuration file error. Check p3.config.")

	global_o = ["REPO", "FTP_1", "FTP_2", "PRIDE"]
	if src == "LOCAL":
		print("Reading local files.")

	elif src == "FTP" and options.get("SOURCE", "FTP_1") == "":
		sys.exit("Reading from FTP. \n ERROR FTP Source is not defined. Check p3.config")        
	elif src == "FTP" and options.get("SOURCE", "FTP_1") != "":
		ftp1 = options.get("SOURCE", "FTP_1")
		ftp2 = options.get("SOURCE", "FTP_2")
		if len(ftp2) != 0:
			print("Reading from {} and {}.".format(ftp1, ftp2)) 
			fetch_ftp(ftp1)
			fetch_ftp(ftp2)

		else:
			print("Reading from {}.".format(ftp1))
			fetch_ftp(ftp1)

		
	elif src == "PRIDE":
		prideID = options.get("SOURCE", "PRIDEID").upper()
		if len(prideID) != 0:
			print("Reading from PRIDE REPOSITORY. Pride ID: {}".format(prideID))
			get_pride(prideID)

		else:
			sys.exit("Reading from PRIDE REPOSITORY. \n Error: Pride Repository is not defined. Check p3.config.")
	else: 
		sys.exit("Error reading p3.config.")

def fetch_ftp(ftp_url):
	print("Downloading from {}".format(ftp_url))
	URL = urlparse(ftp_url)
	url_host = URL.netloc
	url_path = URL.path
	ftp = FTP(url_host)
	ftp.login()
	ftp.cwd(url_path)
	ftp.retrlines('LIST')
	filenames = ftp.nlst()
	
	filematch = "*.*"
	print("Downloading...")
	for filename in ftp.nlst(filematch)[0:3]:
		if os.path.isfile(filename):
			print("{} already exists. Next..".format(filename))
		else:
			try:
				fhandle = open(filename, 'wb')
				ftp.retrbinary('RETR ' + filename, fhandle.write)
				fhandle.close()
				print("{} downloaded.".format(filename))
			except:
				print("{} failed to download".format(filename))
				
def unzip_files():
	gzfiles = glob.glob("*.gz")
	if len(gzfiles) != 0:
		for file in gzfiles:
			base = os.path.basename(file)
			dest = base[:-3]
			if not os.path.isfile(dest):
				with gzip.open(base, 'rb') as infile:
					with open(dest, 'wb') as outfile:
						for line in infile:
							outfile.write(line)

def scan_files(ext):
	files = []
	for file in os.listdir(working_dir):
		if file.endswith(ext):
			files.append(file)
	return files



def scan_fasta():
	fasta = []
	ext = (".fasta")
	for file in os.listdir(working_dir):
		if file.lower().endswith(ext):
			fasta.append(file)
	if (len(fasta) > 1): 
		print("Multiple .fasta files found. \n Only {} will be used".format(fasta[0]))
		return fasta[0]
	elif (len(fasta) == 0): 
		print("No Fasta file found. MSGF will not be run.")
	else:
		print("Fasta file: {}".format(fasta[0]))
		return fasta[0]
	
def write_blank_p3(config_file):
	content = """
	# This is configuration file for P3. http://github.com/kristiyanto/p3
	[SOURCE]
	# REPO = LOCAL/FTP/PRIDEID
	REPO = LOCAL
	# If REPO = FTP, at least FTP_1 is required.
	# e.g: 
	# FTP_1 = ftp://massive.ucsd.edu/MSV000079527/sequence/ 
	# FTP_2 = ftp://massive.ucsd.edu/MSV000079527/peak/mgf/
	FTP_1 = ftp://massive.ucsd.edu/MSV000079527/sequence/
	FTP_2 = ftp://massive.ucsd.edu/MSV000079527/peak/mgf/
	# if REPO = PRIDEID, PRIDEID must be defined
	# e.g 
	# PRIDEID = PXD000001
	PRIDEID = 
	
	[MSGF]
	# RUN_MSGF = YES / NO ; if NO, MSGF identification will be skipped.
	# e.g:
	# RUN_MSGF = YES
	RUN_MSGF = YES
	# MSGF_OPTIONS: OPTIONAL. Additional options for identification other than -s -d and -o
	# Check https://omics.pnl.gov/software/ms-gf for more detail
	# e.g
	# MSGF_OPTIONS = -t 2.5Da
	
	[QUANTIFICATION]
	# METHOD = SPECTRUM_COUNT / ITRAQ4 / NONE 
	# If "NONE", quantification will not be run.
	# e.g
	METHOD = SPECTRUM_COUNT
	
	[SPECTRUM_COUNT]
	# Ignored if QUANTIFICATION METHOD is not SPECTRUM_COUNT
	# If left blank will be replaced with default
	# Please check documentation for more detail
	SCORE_TRESHOLD = 7.0
	ERROR_TRESHOLD = 20
	FDR = 0.01
	ITERATION = 5000
	
	[ITRAQ4]
	# Ignored if QUANTIFICATION METHOD is not ITRAQ4
	# If left blank will be replaced with default
	# Please check documentation for more detail
	EVALUE_TRESHOLD = 75
	pNA = 0
	# QUANTIFICATION_METHOD: / trapezoidation / max / sum / SI / SIgi / SIn / SAF / NSAF / count 
	# *CASE SENSITIVE* 
	QUANTIFICATION_METHOD = trap
	COMBINE_BY = mean
	"""
	try:
		with open(config_file, "w") as p3_config:
			p3_config.write(content)
			p3_config.close()
	except:
		sys.exit("Attempt to create a blank p3.config failed. Quitting.")
	sys.exit("Config Template: p3.config was written. Please evaluate configuration file and re-run the container.")

def touch(file, times=None):
	with open(file, 'a'):
		os.utime(file, times)

if __name__ == "__main__": main()

