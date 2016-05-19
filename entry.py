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
from datetime import timedelta


# Set Env
working_dir = "/root/data"
#working_dir = "/Users/Daniel/Desktop/LABELLED"
config_file = "p3.config"
os.chdir(working_dir) if os.path.isdir(working_dir) else sys.exit("{} not found. Please make sure it's properly mounted.".format(working_dir))

def main():
	maxwait = 10800 # Wait for 3 hours
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
	sp_files_ext = (".mzXML", ".mzML", ".MZXML", ".mzml", ".MZML", ".mzml", ".mgf",".mzxml",".ms2",".pkl")

	spectrum = scan_files(sp_files_ext)
	
	if (len(spectrum)!=0):
		print("{} Spectrum files found: {}".format(len(spectrum), spectrum))
	else:
		print("No spectrum files found. Accepted format:{} ".format(sp_files_ext))
	
	Q_METHOD = options.get("QUANTIFICATION","METHOD").upper()
	RUN_MSGF = (options.get("MSGF", "RUN_MSGF").upper() == "YES")
	fasta = scan_fasta()


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

	# if Q_METHOD == "SPECTRUM_COUNT":
	# 	keep_waiting = True
	# 	wait = 0
	# 	lock = "SC_RESULT.txt.tmp"
	# 	if os.path.isfile(lock): sys.exit("Spectrum Count Quantification is already run by other container. Quitting.")
	# 	while keep_waiting:
	# 		mzid_set = scan_files((".mzid"))
	# 		if len(mzid_set) == 0:
	# 			sys.exit("No .mzid file found")
	# 		elif len(scan_files('.mzid.tmp')) != 0: 
	# 			#print("SCAN MZID TMP:" + scan_files('.mzid.tmp'))
	# 			print("MSGF identification is being run by other containers. Waiting...")
	# 			if wait > maxwait:
	# 				keep_waiting = False
	# 				sys.exit("Had waited for {}s. Quitting now. \n Perhaps remove .tmp files and rerun containers?".format(maxwait))
	# 			else: wait = wait + 5 
	# 			time.sleep(5)
	# 		else:
	# 			touch(lock)
	# 			scquant(mzid_set, options)
	# 			os.remove(lock)
	# 			keep_waiting = False

	# if Q_METHOD == "ITRAQ4":          
	mzid = scan_files(".mzid")
	for m in mzid: 
		out2 = os.path.splitext(m)[0] + ".txt"	
		if not os.path.isfile(out2): itraq(m, options)
	# 	lock_files = scan_files(".rda.tmp")
	# 	wait = 0
	# 	while len(lock_files) != 0:
	# 		print("Identification is being run by other containers. Waiting...")
	# 		if wait > maxwait:
	# 			keep_waiting = False
	# 			sys.exit("Had waited for {} hours. Quitting now. \n Perhaps remove .tmp files and rerun containers?".format(maxwait))
	# 		else: 
	# 			wait = wait + 10
	# 			lock_files = scan_files(".rda.tmp")
	# 			time.sleep(10)
		
	# 	itraq_folding(options)

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
				#print("This is the {}".format(cmd))
				subprocess.call(cmd)
				os.remove(lock)
			except:
				print("MSGF Error: {}".format(file))
			if not os.path.isfile(out2): itraq(out, options)


def itraq(mzid, options):
	ext_set = (".mzXML", ".mzML", ".MZXML", ".mzml", ".MZML", ".mzml")
	opts_set = ("SPECEVALUE_TRESHOLD", "pNA", "QUANTIFICATION_METHOD", "COMBINE_BY")
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
				#print(cmd)
				subprocess.call(cmd)
				os.remove(lock)
			except:
				print("Quantification failed.")
   

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
	R_OPTS = dict(SPECEVALUE_TRESHOLD= 0.01, pNA= 4, QUANTIFICATION_METHOD="sum", COMBINE_BY="mean")
	for k,v in R_OPTS.items():
		try:
			x = options.get("ITRAQ4", k)
			if len(x) != 0 : R_OPTS[k] = x
		except:
			continue
	print("ITRAQ. {} : {}".format(k, R_OPTS[k]))
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
	return MSGF_OPTS

def get_count_opts(options):
	R_OPTS = dict(SPECEVALUE_TRESHOLD= 0.01, pNA= 4, QUANTIFICATION_METHOD="count", COMBINE_BY="mean")
	for k,v in R_OPTS.items():
		try:
			x = options.get("ITRAQ4", k)
			if len(x) != 0 : R_OPTS[k] = x
		except:
			continue
	print("ITRAQ. {} : {}".format(k, R_OPTS[k]))
	return R_OPTS             

			
def scquant(mzid_set, options):
	R_OPTS = dict(SPEC_EVALUE_TRESHOLD= 1e-10)
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
	#
	# MAXWAIT = 10800 # 	Optional. How long the container will wait for other containers 
	# 						to finish when multiple containers are running concurrently. 
	# 
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
	METHOD = SPECTRUM_COUNT
	
	[SPECTRUM_COUNT]
	# Ignored if QUANTIFICATION METHOD is not SPECTRUM_COUNT
	# If left blank will be replaced with default
	# Please check documentation for more detail
	SPEC_EVALUE_TRESHOLD = 1e-10
	# COMBINE_BY will perform feature folding with the specified function. If "SKIP" is given, this task will not be performed.
	# COMBINE_BY: SKIP / mean / median / weighted.mean / sum / medpolish
	COMBINE_BY = mean

	[ITRAQ4]
	# Ignored if QUANTIFICATION METHOD is not ITRAQ4
	# If left blank will be replaced with default
	# Please check documentation for more detail
	SPEC_EVALUE_TRESHOLD = 1e-10
	pNA = 4
	# QUANTIFICATION_METHOD: / trapezoidation / max / sum / SI / SIgi / SIn / SAF / NSAF  
	# *CASE SENSITIVE* 
	QUANTIFICATION_METHOD = max
	# COMBINE_BY will perform feature folding with the specified function. If "SKIP" is given, this task will not be performed.
	# COMBINE_BY: SKIP / mean / median / weighted.mean / sum / medpolish
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

