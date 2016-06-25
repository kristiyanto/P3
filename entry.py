##########################################################################
## P3: Portable Proteomics Pipeline
## DANIEL KRISTIYANTO (daniel.kristiyanto@pnnl.gov)
## JUNE 25, 2016
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
from datetime import timedelta


# Set Env
working_dir = "/root/data"
config_file = "p3.config"
os.chdir(working_dir) if os.path.isdir(working_dir) else sys.exit("{} not found. Please make sure it's properly mounted.".format(working_dir))

def main():
	maxwait = 10800 
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
	sp_files_ext = (".mzXML", ".mzML", ".MZXML", ".mzml", ".MZML", ".mzml",".mzxml")

	spectrum = scan_files(sp_files_ext)
	
	if (len(spectrum)!=0):
		print("{} Spectrum files found: {}".format(len(spectrum), spectrum))
	else:
		print("No spectrum files found. Accepted format:{} ".format(sp_files_ext))
	
	Q_METHOD = options.get("QUANTIFICATION","METHOD").upper()
	RUN_MSGF = (options.get("MSGF", "RUN_MSGF").upper() == "YES")
	fasta = scan_fasta()

	sys.stdout.flush()
	# Run MSGF 
	if fasta and spectrum and RUN_MSGF:
		print("Running MSGF") 
		msgf(spectrum, fasta, options)
	else:
		print("Skipping MSGF identification...")

	try:
		maxwait = options.get("SOURCE", "MAXWAIT")
		print("Max wait updated to {} seconds".format(maxwait))
	except:
		pass

	mzid = scan_files(".mzid")
	for m in mzid: 
		out2 = os.path.splitext(m)[0] + ".txt"	
		if not os.path.isfile(out2): itraq(m, options)

	stop_time = time.time()
	elapsed_time = time.time() - start_time
	print("Started at: {} Finished at: {} Elapsed: {}".format(str(start_time), str(datetime.now()), str(timedelta(seconds=elapsed_time))))
		
	
################################################ FUNCTIONS  ################################################

def msgf(spectrum, fasta, options, *opt):
	#is_itraq = options.get("QUANTIFICATION","METHOD").upper() == "ITRAQ4"
	#is_count = options.get("QUANTIFICATION","METHOD").upper() == "SPECTRUM_COUNT"
	tags = get_msgf_opts(options)
	for file in spectrum:
		out = os.path.splitext(file)[0] + ".mzid"
		lock = out + ".tmp"
		out2 = os.path.splitext(file)[0] + ".txt"
		#print("(MSGF) Working on: {}".format(file))

		if os.path.isfile(out):
			print("{} already exists.".format(out))
		elif os.path.isfile(lock):
			print("Other container is working on {}. Next...".format(lock))
		else:
			try:
				cmd = ['java', '-Xmx3500M', '-jar', '/root/MSGFPlus.jar', '-s', file, '-d', fasta, '-o', out]
				if len(tags) != 0: cmd.extend(tags)
				touch(lock)
				subprocess.call(cmd)
				os.remove(lock)
			except:
				print("MSGF Error: {}".format(file))
			if not os.path.isfile(out2): itraq(out, options)
		sys.stdout.flush()


def itraq(mzid, options):
	ext_set = (".mzXML", ".mzML", ".MZXML", ".mzml", ".MZML", ".mzml")
	opts_set = ("SPEC_EVALUE_TRESHOLD", "pNA", "QUANTIFICATION_METHOD", "COMBINE_BY")
	sp_files = scan_files(ext_set)
	lock = mzid[:-5] + ".rda.tmp"
	out = mzid[:-5]+".rda"
	if os.path.isfile(out) or os.path.isfile(lock):
		print("{} already exists... next.".format(lock))
	else:
		cmd = ['Rscript', "/root/itraq.R"]
		mzml = ""
		for ext in ext_set:
			m = mzid[:-5] + ext
			if m in sp_files:
				mzml = m
				cmd.append(mzml)
				cmd.append(mzid)
				if(options.get("QUANTIFICATION","METHOD").upper() == "SPECTRUM_COUNT"):
					opts = get_count_opts(options)	
				else:
					opts = get_itraq_opts(options)
				for k in opts_set:
					cmd.append(str(opts[k]))
				break
		cmd.append(lock[:-4])
		#print(" Working on {}".format(ext_set))
		if len(mzml)!=0: 
			print("Quantifying:" + mzml)
			try:
				touch(lock)
				subprocess.call(cmd)
				os.remove(lock)
			except:
				print("Quantification failed.")
	sys.stdout.flush()
   

def itraq_folding(options):
	try: c_by = options.get("ITRAQ4","COMBINE_BY").lower()
	except: c_by = mean
	finally:
		print("")
		#cmd = ['Rscript', "/root/fold.R", c_by]

def get_pride(prideID):
	cmd = ['Rscript', "/root/pride.R", prideID]
	subprocess.call(cmd)

	
def get_itraq_opts(options):
	R_OPTS = dict(SPEC_EVALUE_TRESHOLD= 0, pNA= 4, QUANTIFICATION_METHOD="sum", COMBINE_BY="SKIP")
	for k,v in R_OPTS.items():
		try:
			x = options.get("ITRAQ4", k)
			if len(x) != 0 : R_OPTS[k] = x
		except:
			continue
	print("ITRAQ. {} : {}".format(k, R_OPTS[k]))
	sys.stdout.flush()
	return R_OPTS             

def get_msgf_opts(options):
	MSGF_OPTS = []
	tmp = ("-t", "-m", "-inst", "-e", "-ti", "-ntt", "-tda","-minLength","-maxLength","-minCharge","-maxCharge","-n","-thread","-mod","-minNumPeaks","-addFeatures")
	for k in tmp:
		try:
			v = options.get("MSGF", k)
		except:
			continue
		if k == "-mod" and len (v) >0 and not os.path.isfile(v):
			sys.exit("MSGF Option: {} could not be found.".format(v)) 
		if len(v) > 0: 
			MSGF_OPTS.append(k)
			MSGF_OPTS.append(v)
	#print(MSGF_OPTS)
	sys.stdout.flush()
	return MSGF_OPTS

def get_count_opts(options):
	R_OPTS = dict(SPEC_EVALUE_TRESHOLD= 0, pNA= 4, QUANTIFICATION_METHOD="count", COMBINE_BY="SKIP")
	for k,v in R_OPTS.items():
		try:
			x = options.get("SPECTRUM_COUNT", k)
			if len(x) != 0 : R_OPTS[k] = x
		except:
			continue
	print("ITRAQ. {} : {}".format(k, R_OPTS[k]))
	sys.stdout.flush()
	return R_OPTS             

			
def scquant(mzid_set, options):
	R_OPTS = dict(SPEC_EVALUE_TRESHOLD= 0)
	opts_set = ("SPEC_EVALUE_TRESHOLD")
	for k,v in R_OPTS.items():
		try: 
			x = options.get("SPECTRUM_COUNT", k)
			if len(x) != 0 : R_OPTS[k] = x
		except:
			continue
	cmd = ['Rscript', "/root/scquant.R"]
	for k in opts_set:
		cmd.append(str(R_OPTS[k]))
	if os.path.isfile("SpectrumCount.txt"):
		sys.exit("SpectrumCount.txt already exsists. Quitting.")
	else:
		print(cmd)
		subprocess.call(cmd)
	sys.stdout.flush()
			  
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
			if not global_o: sys.exit("Configuration error: SOURCE. Check p3.config.")
			if not msgf_o: sys.exit("Configuration error: MSGF. Check p3.config.")
			if not quant_o: sys.exit("Configuration error: QUANTIFICATION. Check p3.config.")
	except:
		sys.exit("Error reading p3.config.")
				
def get_files(options):
	try: 
		src = options.get("SOURCE", "REPO").upper()
		if src == "FTP":
			try:
				ftp1 = options.get("SOURCE", "FTP_1")
			except:
				sys.exit("Failed to read FTP_1 option in config file.")
			try:
				ftp2 = options.get("SOURCE", "FTP_1")
			except:
				pass
	except:
		sys.exit("Configuration file error. Check p3.config.")

	if src == "LOCAL":
		print("Reading local files.")

	elif src == "FTP" and options.get("SOURCE", "FTP_1") == "":
		sys.exit("Reading from FTP. \n ERROR FTP Source is not defined. Check p3.config")        
	elif src == "FTP" and options.get("SOURCE", "FTP_1") != "":
		if len(ftp2) != 0:
			print("Reading from {} and {}.".format(ftp1, ftp2))
			sys.stdout.flush() 
			fetch_ftp(ftp1)
			fetch_ftp(ftp2)

		else:
			print("Reading from {}.".format(ftp1))
			fetch_ftp(ftp1)

			
	elif src == "PRIDE":
		prideID = options.get("SOURCE", "PRIDEID").upper()
		print(prideID)
		if len(prideID) != 0:
			if not os.path.isfile("pride_url.txt"):
				print("Reading from PRIDE REPOSITORY. Pride ID: {}".format(prideID))
				get_pride(prideID)
			with open("pride_url.txt", "r") as f:
				ftp1 = f.readline().rstrip()
			fetch_ftp(ftp1)
		else:
			sys.exit("Reading from PRIDE REPOSITORY. \n Error: Pride Repository is not defined. Check p3.config.")
		try:
			os.remove("pride_url.txt")
		except:
			pass

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
	#ftp.retrlines('LIST')
	filenames = ftp.nlst()
	
	for filematch in ("*.fasta","*.mzXML", "*.mzML", "*.MZXML", "*.mzml", "*.MZML", "*.mzml","*.mzxml","*.gz*"):
		print("Downloading...")
		for filename in ftp.nlst(filematch):
			print(filename + " ...")
			sys.stdout.flush()

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
		sys.stdout.flush()	

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
	# REPO = LOCAL/FTP/PRIDE
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
	# -t = 10ppm 
	# -m = 0 
	# -inst = 1 
	# -e = 1 
	# -ti = -1,2 
	# -ntt = 2 
	# -tda = 1 
	# -minLength = 6 
	# -maxLength = 50 
	# -minCharge = 2 
	# -maxCharge = 5 
	# -n = 1 
	# -thread = 7 
	# -mod = MSGFDB_Mods.txt 
	# -minNumPeaks = 5 
	# -addFeatures = 1
	
	[QUANTIFICATION]
	# METHOD = SPECTRUM_COUNT / ITRAQ4 / NONE 
	# If "NONE", quantification will not be run.
	# e.g
	# METHOD = SPECTRUM_COUNT
	METHOD = SPECTRUM_COUNT
	
	[SPECTRUM_COUNT]
	# This section is ignored if QUANTIFICATION METHOD is not SPECTRUM_COUNT
	# If left blank filtering will not be performed
	SPEC_EVALUE_TRESHOLD = 1e-10

	# COMBINE_BY will perform feature folding with the specified function. If "SKIP" is given, this task will not be performed.
	# COMBINE_BY: SKIP / mean / median / weighted.mean / sum / medpolish
	# For spectrum count use "sum" or "SKIP" (no counting).
	COMBINE_BY = sum

	[ITRAQ4]
	# This section is ignored if QUANTIFICATION METHOD is not ITRAQ4
	# If left blank filtering will not be performed
	SPEC_EVALUE_TRESHOLD = 1e-10
	
	# pNA = Should scans with NA's in any of the reporter ions be retained?
	# 4: Scans with NA in any of the reporter ions will be removed
	# 0: All scans will be retained
	pNA = 4

	# QUANTIFICATION_METHOD: / trapezoidation / max / sum / SI / SIgi / SIn / SAF / NSAF  
	# *CASE SENSITIVE* 
	QUANTIFICATION_METHOD = max

	# COMBINE_BY will perform feature folding with the specified function. If "SKIP" is given, this task will not be performed.
	# COMBINE_BY: SKIP / mean / median / weighted.mean / sum / medpolish
	COMBINE_BY = mean


	# OPTIONAL
	# MAXWAIT = 10800 # 	Optional. How long the container will wait for other containers 
	# 						to finish when multiple containers are running concurrently. 


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

