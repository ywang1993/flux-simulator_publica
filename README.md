# Flux-simulator pipeline
### Date: 04.28.2020


## Table of Contents
- [Dependencies](#Dependencies)
- [Pipeline](#Pipeline)
    * [run Flux Simulator](#Flux-Simulator)
    * [mapping](#mapping)
    * [run rMATS](#run-rmats)
    * [calculate PSI](#calculate-PSI-from-expression-profile)
- [Contact](#Contact)



## Dependencies
#### flux simulator
1. Download and install [flux-simulator](http://confluence.sammeth.net/display/SIM/2+-+Download)
    Note: Mannuals for flux simulator can be found here:
    - Official [confluence website](http://confluence.sammeth.net/display/SIM/Home)
    - [wiki for flux simulator](http://fluxcapacitor.wikidot.com/simulator)
    - For user-defined expression profiles, follow instructions [here](http://confluence.sammeth.net/display/SIM/flux+simulator+.pro+file).

2. Install [JDK](http://www.oracle.com/technetwork/java/javase/downloads/jdk8-downloads-2133151.html). Or, use pre-installed java.

#### rmats
#### mapping softwares (e.g. STAR)

## Pipeline
### Flux Simulator
#### Required input
- parameter file `.PAR`
    [Example parameter file](paraFiles/example_unmodified.PAR)
- reference file `.GTF`
- genome fasta by chromosome


#### Run as a whole process 
Example parameter file [example_unmodified.PAR](paraFiles/example_unmodified.PAR)
```
flux-simulator -p XX.par
```


#### Run each step seperately
*This enable us to modify expression profile using real data*, instead of using simulated expression profile.

1. generate expression profile by [expression.PAR](./paraFiles/expression.PAR)
    ```
    flux-simulator -p XX.par -x
    ```
2. Optional: modify expression profile (column 6).
    Expression profile column 6 represent the number of initial molecules. It can be replaced by counts from real data, for example, kallisto output `abundance.tsv` (column 4).
3. generate library by [lib.PAR](./paraFiles/lib.PAR)
    ```
    flux-simulator -p XX.par -l
    ```
4. get sequences by [seq_50M.PAR](./paraFiles/seq_50M.PAR)
    ```
    flux-simulator -p XX.par -s
    ```
### mapping
### run rmats
### calculate PSI from expression profile
```
python ../src/calculatePSI.py --gtf ../reference/gencode.v31lift37.annotation.gtf --pro ../profiles/count_modified.PRO --rmats ../../rmats/post/fromGTF.MXE.txt,../../rmats/post/fromGTF.RI.txt,../../rmats/post/fromGTF.A3SS.txt,../../rmats/post/fromGTF.A5SS.txt,../../rmats/post/fromGTF.SE.txt
```
```
$ python ../src/calculatePSI.py --help
usage: python calculatePSI.py --gtf path/to/gtf.txt --pro path/to/XX.PRO --rmats path/to/from.GTF.SE.txt [options]

optional arguments:
  -h, --help        show this help message and exit
  --od OUTPATH      output folder, default: working folder

required named arguments:
  --gtf FN_GTF      file path to GTF file used for flux simulator
  --pro FN_PRO      file path of expression profile used in flux simulator
  --rmats FN_RMATS  file paths of from.GTF.AS.txt (generated by rmats);
                    seperated by comma
```
### Contact
Yuanyuan Wang <wyynju1993@gmail.com>


