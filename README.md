# Flux-simulator pipeline
### Date: 04.28.2020

### Table of Contents
- [Mannual](#Mannual)
- [Installation](#Installation)
- [Contact](#Contact)


### Mannual
- Official [confluence website](http://confluence.sammeth.net/display/SIM/Home)
- This is the [wiki for flux simulator](http://fluxcapacitor.wikidot.com/simulator)
- For user-defined expression profiles, follow instructions [here](http://confluence.sammeth.net/display/SIM/flux+simulator+.pro+file).


### Installation
- Download and install [flux-simulator](http://confluence.sammeth.net/display/SIM/2+-+Download)
- Install [JDK](http://www.oracle.com/technetwork/java/javase/downloads/jdk8-downloads-2133151.html).

Or, use pre-installed java.

### Pipeline
1. generate expression profile
```
flux-simulator -p XX.par -x
```
1. modify expression profile (column 6)
1. generate library
```
flux-simulator -p XX.par -l
```
1. get sequences
```
flux-simulator -p XX.par -s
```
### Contact
Yuanyuan Wang <wyynju1993@gmail.com>

