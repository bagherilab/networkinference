# Network inference performance complexity: a consequence of topological, experimental, and algorithmic determinants

<sub>LAST UPDATED: 2018-05-18</sub>

## Content overview

- __`_banjo`__ : [BANJO v2.2.0](https://users.cs.duke.edu/~amink/software/banjo/) (_accessed February 28, 2017_)
- __`_genie3`__ : [GENIE3](http://www.montefiore.ulg.ac.be/~huynh-thu/software.html) (_accessed October 3, 2017_)
- __`_mider`__ : [MIDER v2_JAN2015](http://nautilus.iim.csic.es/~gingproc/mider.html) (_accessed October 3, 2017_)
- __`_tigress`__ : [TIGRESS v2.1](http://projets.cbio.mines-paristech.fr/~ahaury/svn/dream5/html/index.html) (_accessed October 3, 2017_)
- __`scripts`__
    + `.pbs` scripts for submitting jobs on [QUEST](http://www.it.northwestern.edu/research/user-services/quest/index.html)
- __`lib`__ : Library functions
- __`pipeline`__ : Main code for simulation, nference, and analysis

## Algorithm edits

Some minor edits were made to the downloaded algorithms. For all cases, extraneous files associated with each of the algorithms (e.g. sample data, documentation) were removed for clarity and simplicity. Note that any associated license files are included.

- __`GENIE3`__
    + commented out `fprintf` calls in `GENIE3.m`
    + commented out `tic` and `toc` in `GENIE3.m`
- __`MIDER`__
    + commented out `fprintf` calls in `mider.m`
    + output of `mider.m` changed from `[Output]` to `con_array` to fit pipeline
    + value of `pb` in `mider.m` changed from `5` to `2`
- __`TIGRESS`__
    + commented out `fprintf`/`disp` calls in `score_edges.m`
    + commented out `tic` and `toc` in `tigress.m`
    + added multiple variable selection check in `lars.m`
    + added empty add variable check in `lars.m`
    + added zero padding for early exit in `stability_selection.m`

## Running on QUEST

Pipeline was written specifically for running on Northwestern's high performance computing core [QUEST](http://www.it.northwestern.edu/research/user-services/quest/index.html). QUEST currently has 20-28 cores per node so scripts are written requesting no more than 20 nodes at a time.

#### File structure

The working directory on QUEST contains the following:

`~/Matlab/`
> `_genie3` 
> `_mider`
> `_tigress`
>  `_banjo`
> `logs/` (standard out and error logs)  
> `results/` (results stored here)
> `*.m` (all Matlab files)  
> `*.pbs` (all submission scripts)

## Pipeline

#### Simulate data

Simulation is very fast and can be run on the login node. Load Matlab using:

```Shell
    >> module load matlab/r2016a
    >> matlab -nosplash -nodisplay -singleCompThread
```

Run the simulations using:

```Matlab
    for i = 1:36
        LM_CONTROLLER(1, i); % generate in silico data
    end

    LM_CONTROLLER(2); % compiles results into single array
    LM_CONTROLLER(3); % generates null models
```

#### Infer networks

Network inference using selection of network inference methods. There are 108 different motif/logic gate/stimulus combinations. 

```Shell
    >> msub -t corr[1-108] run_corr.pbs
    >> msub -t genie3[1-108] run_genie3.pbs
    >> msub -t mider[1-108] run_mider.pbs
    >> msub -t tigress[1-108] run_tigress.pbs
    >> msub -t banjo[1-108] run_banjo.pbs
```

#### Run nulls

Run inference algorithms on the null networks. There are 255 different stimulus/noise/parameter A combinations:

```Shell
    >> msub -t genie3[1-255] run_genie3_nulls.pbs
    >> msub -t corr[1-255] run_corr_nulls.pbs
    >> msub -t mider[1-255] run_mider_nulls.pbs
    >> msub -t tigress[1-255] run_tigress_nulls.pbs
    >> msub -t banjo[1-255] run_banjo_nulls.pbs
```

Concatenate results from nulls into a single matrix.

```Matlab
    for i = 1:15
        LM_CONTROLLER(6, i, ALGORITHM);
    end
```

The script `LM_summarize.py` can be used to check which network inference runs are missing. Run using:

```Shell
    >> ls -LRh results/ > summary.txt
    >> python3 LM_summarize.py
```

#### Calculate metrics

Calculation is relatively fast and can be run in using:

```Matlab
    for i = 1:108 
        LM_CONTROLLER(7, i, ALGORITHM);
    end
```
where `ALGORITHM` is a string denoting which algorithm to process (`GENIE3`, `CORR`, `TIGRESS`, `MIDER`, or `BANJO`).

#### Compile results

Final results are compiled into `.csv` files for the data browser.

```Matlab
    LM_CONTROLLER(8); % save simulation data
    LM_CONTROLLER(9, 0, ALGORITHM); % save inference data
    LM_CONTROLLER(10, 0, ALGORITHM); % save summary data
```