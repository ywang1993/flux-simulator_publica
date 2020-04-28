# Flux-simulator pipeline
### Date: 04.28.2020


### Table of Contents
- [Mannual](#Mannual)
- [Installation](#Installation)
- [Contact](#Contact)


### Mannual
- Official [confluence website](http://confluence.sammeth.net/display/SIM/Home)
- [wiki for flux simulator](http://fluxcapacitor.wikidot.com/simulator)
- For user-defined expression profiles, follow instructions [here](http://confluence.sammeth.net/display/SIM/flux+simulator+.pro+file).


### Installation
- Download and install [flux-simulator](http://confluence.sammeth.net/display/SIM/2+-+Download)
- Install [JDK](http://www.oracle.com/technetwork/java/javase/downloads/jdk8-downloads-2133151.html).

Or, use pre-installed java.



### Pipeline
#### Required input
- parameter file `.PAR`
  [Example parameter file](paraFiles/example_unmodified.PAR)
- reference file `.GTF`
- genome fasta by chromosome


#### Run as a whole process 
flux-simulator -p XX.par


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


### Contact
Yuanyuan Wang <wyynju1993@gmail.com>

