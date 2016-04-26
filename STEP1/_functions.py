# DANIEL.KRISTIYANTO@PNNL.GOV
# MSGF

import os
import sys
import subprocess
import re
from os.path import isfile, join
from urllib.parse import urlparse
from ftplib import FTP
######################## CHECK OPTIONS ###########################
def get_msgf_opts(cfg_file):
	r_options = []
	global_options = SafeConfigParser()
	try: 
		global_options.read(cfg_file)
		sc_options = ["method","msgf_options","source"]
		for o in sc_options:
			print(o +": "+global_options.get("step1",o))
			r_options.append(global_options.get("step1",o))
		return global_options
	except:
		print("Configuration File Error.")

def get_sc_opts(cfg_file):
	r_options = []
	parser = SafeConfigParser()
	try: 
		parser.read(cfg_file)
		sc_options = ["score_treshold","error_treshold", "fdr", "iteration"]
		for o in sc_options:
			print(o +": "+parser.get("spectrum_count",o))
			r_options.append(parser.get("spectrum_count",o))
		return r_options
	except:
		print("Configuration File Error.")


def get_itraq_opts(cfg_file):
	r_options = []
	parser = SafeConfigParser()
	try: 
		parser.read(cfg_file)
		sc_options = ["evalue_treshold","pNA", "quant_method", "combine_by"]
		for o in sc_options:
			print(o +": "+parser.get("itraq4",o))
			r_options.append(parser.get("itraq4",o))
		return r_options
	except:
		print("Configuration File Error.")


######################## SUB PROCESSES ###########################

def write_config(working_dir):
	try:
		with open(os.path.join(working_dir,"p3.config"),"w") as f:
			f.write("[step1]\n")
			f.write("#ftp_src1=\n")
			f.write("#ftp_src2=\n")
			f.write("method=spectrum_count\n")
			f.write("msgf_options=\n")
			f.write("[spectrum_count]\n")
			f.write("score_treshold = 7.0\n")
			f.write("error_treshold = 20\n")
			f.write("fdr = 0.01\n")
			f.write("iteration = 5000\n")
			f.write("#[itraq4]\n")
			f.write("#evalue_treshold = 75\n")
			f.write("#pNA = 0\n")
			f.write("#quant_method = trap\n")
			f.write("#combine_by = mean\n")
			f.close()
		print("Missing p3.config. Default config file written.")
	except:
		print("Missing p3.config, and failed to create one.")
		raise SystemExit




def touch(fname, times=None):
	with open(fname, 'a'):
		os.utime(fname, times)

######################## SUB PROCESSES ###########################

def msgf(spectrum, db, out):
	#print(spectrum, db, out)
	lock = out+".lock"
	if (not os.path.isfile(out)) and (not os.path.isfile(lock)) :
		try:
			touch(lock)
			subprocess.call(['java', '-Xmx3500M', '-jar', 'MSGFPlus.jar', '-s', spectrum, '-d', db, '-o', out])			
			os.remove(lock)
		except:
			print("MZID conversion failed.")

def itraq(r_options):
	cmd = ['Rscript','itraq.R']
	cmd.extend([str(x) for x in r_options])
	runr = subprocess.call(cmd)

def scquant(r_options):
	cmd = ['Rscript','filter.R']
	cmd.extend([str(x) for x in r_options])
	runr = subprocess.call(cmd)


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
	ftp 	= FTP(url_host)
	ftp.login()
	ftp.cwd(url_path)
	ftp.retrlines('LIST')
	filenames = ftp.nlst()
	#print(filenames)
	filematch = "*.*" #Any files
	for filename in ftp.nlst(filematch):
		if os.path.isfile(filename):
			print("File {} already exists. Skipping..".format(filename))
		else:
			try:
				fhandle = open(filename, 'wb')
				print("Downloading {}".format(filename))
				ftp.retrbinary('RETR ' + filename, fhandle.write)
				fhandle.close()
			except:
				print("{} failed.".format(filename))
	print("FTP Download done.")
	os.chdir(work_dir)
	return URL

def download_ftp(input_csv):
	with open(input_csv) as f:
		reader = csv.reader(f)
		header = next(reader)
		for row in reader:
			for src in row:
				print(src)
				if (check_url(src)):
					get_ftp(src)

def unzip_gz(working_dir):
	spectrum_tmp = []
	for src_name in glob.glob(os.path.join(working_dir, '*.gz')):
		base = os.path.basename(src_name)
		dest_name = os.path.join(working_dir, base[:-3])
		if not os.path.isfile(dest_name):
			spectrum_tmp.append(dest_name)
			with gzip.open(src_name, 'rb') as infile:
				with open(dest_name, 'wb') as outfile:
					for line in infile:
						outfile.write(line)

