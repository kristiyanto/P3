
# P3: Portable Proteomics Pipeline
P3 is an open source pipelines for Mass Spectrometry Data Pre-Processing and Quantification, wrapped in Docker containers. 

Identification is performed by utilizing [MSGF+ tool developed by PNNL](https://omics.pnl.gov/software/ms-gf), and Quantification is conducted by utilizing [MSnbase](http://bioconductor.org/packages/release/bioc/html/MSnbase.html), a bioconductor package by Laurent Gatto _et. all._

P3 uses [Bioconductor Proteomics](https://github.com/Bioconductor/bioc_docker) image as the base image.

### Requirements
To run P3, Docker engine must be installed. [Click here](https://docs.docker.com/engine/installation/) for a detailed information to install Docker engine on various operating system including Windows and MacOS.

Protein identification and quantification is a computationally intensive process. Depending on the size of the data, at least 4Gb available memory on the Docker Machine is required. Click here for more information on increasing the memory allocation for Docker engine on VirtualBox machine for MacOS and Windows Users.

![Adjusting RAM allocation for Docker Machine](media/ram.png)

### Input 
The container takes Mass Spectrometry  (MS2) data (\*.mzml, \*.mgf, \*.mzxml, \*.ms2, \*.pkl) and peptide sequences files (\*.fasta). It is advised to put only the files to be computed in the folder as a separate sessions (e.g. different project for each folder).

Alternatively, if the files are accessible on FTP server, a CSV file with the information about the files can as the input. The container reads the information, download the files and run the computation.

Sample of format for the input csv file:

| Spectrum_Files | Database_File |
| --- | --- |
| ftp://link/to/the/folder | ftp://link/to/file.fasta |

In either case, the folder must be mounted to the container's "/root/data" (using ```-v``` tag). [Click here](http://container-solutions.com/understanding-volumes-docker/) for more information about Docker Volumes. 

### Output
For each MS2 file provided, a MZID file containing the protein identification is generated. When the pipelines completed, tab-delimited txt files ```LabelledQuant.txt``` or ```LabelFreeQuant.txt```, and ```evalue.txt``` are also generated. 

```LabelledQuant.txt``` or ```LabelFreeQuant.txt``` is the quantified result (output) for either Labelled or LabelFree respectively. ```evalue.txt``` contains more detailed information including Spectrum number, Pep Sequence, and e-value prior to the agggregation proccess that may be useful for further analysis.

### Running the container
Make sure to put the files to compute on a separate folder and mounted to the Docker Engine. For Windows and MacOS ```/Users``` or ```/c/Users``` is the default. 

* For Windows users, it usually starts with /c/Users/
* For MacOS users, it usually starts with /Users/
* For Linux users, there is no restriction.

![Linking volumes between host OS and VM (Docker Engine)](media/vmvolume.png)

[Click here](https://www.virtualbox.org/manual/ch04.htmlftp) for more information about Folder Sharing from VM and host machine.

To run the container (LabelFree/Spectrum Count):

```
docker pull kristiyanto/msgf:spectrumcount
docker run --rm -v /c/Users/path/to/your/file:/root/data kristiyanto/msgf:spectrumcount
```

To run the container (Labelled):

```
docker pull kristiyanto/msgf:itraq
docker run --rm -v /c/Users/path/to/your/file:/root/data kristiyanto/msgf:itraq
```

### Contact Information

You are invited to contribute for new features, updates, fixes by sending pull requests.
