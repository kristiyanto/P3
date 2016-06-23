
# P3: Portable Proteomics Pipeline

P3 is a docker container for mass-spectometry data pre-processing pipelines. The pipelines consists of protein identification (using [MSGF+ tool developed by PNNL](https://omics.pnl.gov/software/ms-gf)), and quantification (Bioconductor / [MSnbase](http://bioconductor.org/packages/release/bioc/html/MSnbase.html). 


#### Input
The container takes mass-spectometry raw files (```*.mzml```, ```*.mzml.gz```, ```*.mgf```, ```*.mzxml```, ```*.ms2```, ```*.pkl```) and a peptide sequences file (```*.fasta```). Folder must be mounted to the container's "/root/data" (e.g. using ```-v``` tag). [Click here](http://container-solutions.com/understanding-volumes-docker/) for more information about Docker Volumes. 

The files can be retrieved from a mounted local storage ![Linking volumes between host OS and VM (Docker Engine)](media/vmvolume.png)
[Click here](https://www.virtualbox.org/manual/ch04.htmlftp) for more information about Folder Sharing from VM and host machine.

Alternatively, the container can also retrieve the files from the internet by providing either FTP address or PrideID in the ```p3.config``` file.

#### Configuration Files
If not provided, ```p3.config``` template will be created. The configuration file navigate how the pipeline should be run.

### Running the container
```
docker pull kristiyanto/p3
docker run --rm -v /path/to/files:/root/data kristiyanto/p3
```

#### Output
Once process is done, following files will be created:
- *.txt 	: a tab delimited file of the result (spec-evalue, identified peptides, quantification results, etc.)
- *.rda 	: R object of the results. Ensure MSnBase package is installed prior importing.
- *.mzid 	: results from MSGF+


### Requirements
To run P3, Docker engine must be installed. [Click here](https://docs.docker.com/engine/installation/) for a detailed information to install Docker engine on various operating system including Windows and MacOS.

Protein identification and quantification is a computationally intensive process. Depending on the size of the data, at least 4Gb available memory on the Docker Machine is required. Click here for more information on increasing the memory allocation for Docker engine on VirtualBox machine for MacOS and Windows Users.

![Adjusting RAM allocation for Docker Machine](media/ram.png)


### Contact Information

You are invited to contribute for new features, updates, fixes by sending pull requests.
