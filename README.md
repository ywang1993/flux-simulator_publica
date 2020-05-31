# Flux-simulator pipeline
#### Date: 04.28.2020
#### update: 05.25.2020


## Table of Contents
- [Dependencies](#Dependencies)
- [Pipeline](#Pipeline)
    * [run Flux Simulator](#Flux-Simulator)
    * [mapping](#mapping)
    * [run rMATS](#run-rmats)
    * [calculate PSI](#calculate-PSI-from-expression-profile)
- [Contact](#Contact)



## Dependencies
1. Download and install [flux-simulator](http://confluence.sammeth.net/display/SIM/2+-+Download). Mannuals can be found here:
    - Official [confluence website](http://confluence.sammeth.net/display/SIM/Home)
    - [wiki for flux simulator](http://fluxcapacitor.wikidot.com/simulator)
    - For user-defined expression profiles, follow instructions [here](http://confluence.sammeth.net/display/SIM/flux+simulator+.pro+file).

2. Install [JDK](http://www.oracle.com/technetwork/java/javase/downloads/jdk8-downloads-2133151.html). Or, use pre-installed java.

## Pipeline
### Flux Simulator
- Reqired Input
    * reference file `reference/gencode.v31lift37.annotation.gtf`
    * genome fasta by chromosome `reference/GRCh37_by_chrom/`
    * parameter file (`.PAR`) in folder `paraFiles`
        - [example 1 .PAR](paraFiles/example_unmodified.PAR)
        - [example 2 .PAR](paraFiles/example_from_NBT.PAR) used in [Sailfish](https://www.nature.com/articles/nbt.2862) paper

#### Run as a whole process 
Using parameter file [example_unmodified.PAR](paraFiles/example_unmodified.PAR) to generate 30M paired-end reads in `fastq` folder
```
flux-simulator -p XX.par
```


#### Run each step seperately
*This enable us to modify expression profile using real data*, instead of the simulated expression profile.

1. generate expression profile [count_30M.PRO](./profiles/count_30M.PRO) by [expression.PAR](./paraFiles/expression.PAR) 
    ```
    flux-simulator -p XX.par -x
    ```
2. Optional: modify expression profile [count_30M.PRO](./profiles/count_30M.PRO) (column 6) to [count_modified.PRO](./profiles/count_modified.PRO).  **Column 6** represent the number of initial molecules. It can be replaced by counts from real data (e.g. TPM, normalized by transcript length). 
   - **my example**: column 6 replaced by <img src="https://render.githubusercontent.com/render/math?math=TPM_{kallisto}*50"> from real data (293 T cell line) [abundance.tsv](./profile/abundance.tsv)

3. generate library [count_modified.LIB](./libraries/count_modified.LIB) by [lib.PAR](./paraFiles/lib.PAR) 
    ```
    flux-simulator -p XX.par -l
    ```
4. get sequences [50M.fastq](./fastq/50M.fastq) by [seq_50M.PAR](./paraFiles/seq_50M.PAR) 
    ```
    flux-simulator -p XX.par -s
    ```

### calculate PSI from expression profile
```
python src/calculatePSI_fromFlux.py --gtf ./reference/gencode.v31lift37.annotation.gtf --pro ./profiles/count_modified.PRO -od ./psi
```
```
$ python src/calculatePSI_fromFlux.py --help
usage: python calculatePSI_fromFlux.py --gtf path/to/gtf.txt --pro path/to/XX.PRO [options]

optional arguments:
  -h, --help    show this help message and exit
  --od OUTPATH  output folder, default: working folder

required named arguments:
  --gtf FN_GTF  file path to GTF file used for flux simulator
  --pro FN_PRO  file path of expression profile used in flux simulator
```
### Contact
Yuanyuan Wang <wyynju1993@gmail.com>


