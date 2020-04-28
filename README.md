# Flux-simulator pipeline
### Date: 04.28.2020

### Table of Contents
- [Install](#Install)
- Pipeline

### Install
- Download [flux-simulator](http://confluence.sammeth.net/display/SIM/2+-+Download)
- Install [JDK](http://www.oracle.com/technetwork/java/javase/downloads/jdk8-downloads-2133151.html)
Or, use pre-installed java.
- For user-defined expression profiles, see [this](http://confluence.sammeth.net/display/SIM/flux+simulator+.pro+file)




### Table of Contents
- [What it does?](#what-it-does)
- [Installation](#installation)
- [Testing](#testing)
- [Configuration](#configuration)
- [Usage](#usage)
- [Contact](#contact)

### What it does?
Snakemake is a workflow management package written in python3. It automatically manages and submits jobs
for you, so you can be just like [this](https://s-media-cache-ak0.pinimg.com/originals/1e/b9/92/1eb992ab12cf9376eb762c63c3a51911.jpg).

![Snakemake](doc/snakemake.png)

To get ready, you only need a configuration file that specifies your file names and comparisons, and a set of fastq files.
The pipeline will generate STAR alignment, kallisto expression quantification, rMATS splicing anlaysis, and a BigWiggle file
for visualizing the RNA-seq data.

### Installation
#### Install Snakemake
First, install the [Snakemake](http://snakemake.readthedocs.io/en/latest/index.html) python package, follow instructions [here](http://snakemake.readthedocs.io/en/latest/getting_started/installation.html).

If you are using python2 environment previously, I strongly encourage you to build a `py36` environment using 
[Anaconda](https://www.anaconda.com/download/). Anaconda will automatically manage different versions of python and also their dependencies.

Follow the instructions [here](https://conda.io/docs/user-guide/tasks/manage-python.html#installing-a-different-version-of-python) for
how to build a python3 environment.

#### Install software
Next, put the following software in your environment variables. If you don't have the software, click the link to download.
  - [STAR](#)
  - [kallisto](#)
  - [rMATS](#)
  - [Darts](#) *dev*
  - [sleuth](#) *dev*
  
If you don't know how to put things in environment variables, see [here](https://www.cyberciti.biz/faq/set-environment-variable-unix/).

### Testing
To test if you have successfully installed `Snakemake`, type following commands in your shell:
```
source activate py36 ## only if your default python is python2
PROJECT='test_snakemake' snakemake -n
```

You should see a list of task being printed on the screen.

### Configuration
You can open the template configuration file in "test_snakemake/config/config.yaml". Now let's see what's in there, step by step.

#### The genome index
First there is a chunk of genomic indices; you can either use mine or change to your own paths.
```
genome_index:
    hg19:
        star_idx: /u/nobackup/yxing/NOBACKUP/frankwoe/hg19/star_idx_gencode_v19
        kallisto_idx: /u/nobackup/yxing/NOBACKUP/frankwoe/hg19/kallisto_pctx_idx_gencodeV19/hg19_pc_tx_kal.idx
        gtf: /u/nobackup/yxing/NOBACKUP/frankwoe/hg19/gencode.v19.annotation.gtf
    mm10:
        star_idx: /u/nobackup/yxing/NOBACKUP/frankwoe/mm10/star_idx_gencodeM13
        gtf: /u/nobackup/yxing/NOBACKUP/frankwoe/mm10/gencode.vM13.annotation.gtf
```
Note the `genome` you specified for your pipeline run has to be found in the `genome_index`. For example, if you want to analyze your data in `hg38`,
you have to add the corresponding index under a new entry named `hg38` here.

#### The sample names
The next part is for specifying 'sample name' to 'fastq name' relationships. 
```
sample_dict:
    s1: ['s1_lane1', 's1_lane2']
    s2: ['s2_lane1', 's2_lane2']
```
Human-readable description for above: we have two samples named 's1' and 's2'; for 's1', look for fastq files from lane1, i.e 's1_lane1_1.fastq.gz' and 's1_lane1_2.fastq.gz', as 
well as from lane 2, i.e. 's1_lane2_1.fastq.gz' and 's1_lane2_2.fastq.gz'.

Note: currently, the pipeline only recognizes pattern of '*_1.fastq' and '*_2.fastq' for paired-end reads; that is, naming your files like '*_R1.fastq' and '*_R2.fastq' will throw
a 'FileNotFound' error.

#### Speicifying the comparisons
Now we have all the sample names, it's time to specify which should be compared against which:
```
rmats:
    readLength: 101
    comparison:
        s1-s2:          # this is your comparison name; make it informative
            group1:
                ["s1"]  # for multiple replicates, just append in the list
            group2:
                ["s2"]
```
You can change the comparison name here, i.e. "s1-s2"; and also add more comparisons if you need to. All specified comparisons will be run.


### Contact
If you still have problems/questions, or would like to incorporate certain modules, let me know!

Zijun Zhang <zj.z@ucla.edu>
