# Flux-simulator pipeline
### Date: 04.28.2020


## Table of Contents
- [Install Softwares](#Install Softwares)
- [Pipeline](#Pipeline)
    - [Flux Simulator Pipeline](#Flux Simulator Pipeline)
- [Contact](#Contact)



## Install Softwares
#### flux simulator
1. Download and install [flux-simulator](http://confluence.sammeth.net/display/SIM/2+-+Download)
2. Install [JDK](http://www.oracle.com/technetwork/java/javase/downloads/jdk8-downloads-2133151.html). Or, use pre-installed java.

    Note: Mannuals for flux simulator can be found here:
    - Official [confluence website](http://confluence.sammeth.net/display/SIM/Home)
    - [wiki for flux simulator](http://fluxcapacitor.wikidot.com/simulator)
    - For user-defined expression profiles, follow instructions [here](http://confluence.sammeth.net/display/SIM/flux+simulator+.pro+file).


#### rmats

## Pipeline
### Flux Simulator Pipeline
#### Required input
- parameter file `.PAR`
  [Example parameter file](paraFiles/example_unmodified.PAR)
- reference file `.GTF`
- genome fasta by chromosome


#### Run as a whole process 
```
flux-simulator -p XX.par
```


#### Run each step seperately
This enable us to modify expression profile using real data, instead of using simulated expression profile.

1. generate expression profile
    ```
    flux-simulator -p XX.par -x
    ```
2. modify expression profile (column 6)
3. generate library
    ```
    flux-simulator -p XX.par -l
    ```
4. get sequences
    ```
    flux-simulator -p XX.par -s
    ```
### RUN rmats
### calculate PSI from expression profile
    ```
    python2 src/calculatePSI.py [fn_gtf] [fn_pro] [fn_rmats]
    ```

### Contact
Yuanyuan Wang <wyynju1993@gmail.com>


