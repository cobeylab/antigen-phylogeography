## Introduction

This is a modification of Trevor Bedford's [antigen](https://github.com/trvrb/antigen) code that allows ecological parameters to vary between populations.

-------------------------------------------

## Building

The program requires Python 2.7.x and Java 1.7+ to compile and run.

The program can be built using:

```{bash}
[path-to]/antigen-phylogeography/make.py
```	

This will create a `dist` (distribution) directory that can be used as the basis for
a simulation run or a parameter sweep. Inside `dist`:

* `run.py`: a Python script that calls `java` to actually run the simulation. It can be
modified to change, e.g., the amount of memory Java uses.
* `example_run`: a directory containing an example parameters file and a copy of a
set of analyses in Mathematica
* `example_sweep`: a directory containing tools to run a parameter sweep of the model on
a SLURM cluster
* `lib`: the `.jar` files containing the simulation code as well as third-party libraries
* `src`: a copy of the source code, for the sake of knowing what version of the code
generated your results.
	
## Parameters

Parameters are loaded from the file `parameters.json` into an object defined by the
`Parameters` class defined in `src/antigen/Parameters.java`. You must specify all
parameters for every simulation run; no defaults are provided in order to prevent
confusion about which parameters are being used. If a parameter is left out, the
simulation will throw an exception.

The one exception to this rule is the `randomSeed` parameter: if left out of the file,
one will be generated at runtime. As an extra verification, all loaded parameters will be
written into `parameters_out.json`, including the random seed actually used for the
simulation, so that individual runs can be replicated.

## Running a single run

To run a single simulation, execute `run.py` in a directory containing a `parameters.json`
file. This will generate all output in that directory.

## Running a parameter sweep

To run a parameter sweep, do the following:
1. Make a copy of the generated `dist` directory to the cluster.
2. Modify `constant_parameters.json` to include the parameter values shared among all runs.
3. Modify the `run_sweep.py` script to perform the desired sweep. Specifically, you need
to modify the `generateJobs()` function to yield a SLURM job name, directory name, and
dictionary containing all the parameter values that vary from run to run.
4. Modify `job.sbatch` to match the memory requirements for the job.
5. To start, the variable `DRY` is set to `True`, allowing you to do a "dry run" of the
sweep where output directories and parameter files are generated, but jobs are not
actually submitted. Try it by running `run.py`.
6. If the dry-run directory hierarchy and parameter files look good, you can submit the
jobs by removing the test `results` directory, setting `DRY` to `False`, and running
`run.py` again.

## Output

The simulation will output a timeseries of region-specific prevalence and incidence to
[`out.timeseries`](example/out.timeseries).  It will also sample viruses periodically and output their 
geographic and antigenic locations to [`out.tips`](example/out.tips) and a tree connecting these samples 
to [`out.branches`](example/out.branches).  This file contains pairs of viruses, child and parent, 
representing nodes in a genealogy.  Average values are output to [`out.summary`](example/out.summary).

The vaccine strains will be periodically written to `vaccine.timeseries`. They are also written
to a table in `samples.sqlite`, which may later be expanded to include all simulation output.

If you have Mathematica, [`antigen-analysis-noFig.m`](example/antigen-analysis-noFig.m) will generate a many summary statistics for the simulation, including the proportion of the trunk in each deme and the antigenic lead of each deme.

## Assembling the output of a parameter sweep

An example script `make_database.py` is provided to show how to gather together results
from output runs into a single SQLite database. The script `summarize.py` adds an extra
database table containing some summary statistics for each combination of parameters.

The script `write_params_to_sqlite.py` will record all the parameters in parameters_out.json from every simulation into a single SQLite database, in order to mitigate the need to keep a large number of files.

## Memory

The `-Xmx22G` is required, because as an individual-based model, the memory requirements are
typically quite large. Each host requires a minimum of 40 bytes of memory, plus 8 bytes per
Phenotype recorded in its immune history.  If the yearly attack rate is 15% and the host life span
is 30 years, at equilibrium the average size of the immune history will be 4.5 references.  This
gives memory usage of: population size x 76 bytes.  With 45 million hosts (used in the default
parameters), the equals 3.42GB. Large epidemics, either stochastic or due to seasonal forcing, require even more memory.

In addition to hosts and immune histories, the simulation tracks the virus genealogy through
VirusTree.  This is harder to profile, and will continually grow in memory usage throughout the
simulation.  With the default parameters, VirusTree takes 5.5 MB at the end of a simulated year and
may up to 220 MB at the end of the default 40 simulated years.

Memory can be easily profiled by calling `jmap -histo <PID>`.

-------------------------------------------

Copyright Trevor Bedford, Sarah Cobey, Frank Wen, and Ed Baskerville, 2010-2016. Distributed under the GPL v3.