ctmixtures
==========

Cultural transmission models with heterogeneous mixtures of social learning modes

**BUILD STATUS**:  [![Build Status](https://travis-ci.org/mmadsen/ctmixtures.svg?branch=master)](https://travis-ci.org/mmadsen/ctmixtures)


## Purpose ##

The goal of the model is to simulate the population-level effects of realistic mixtures of social learning rules. I 
conjecture that at the population level, mixtures of social learning rules will often converge to appear unbiased, 
and thus the result of "neutral" or "random copying" processes -- in other words, will meet the quantitative expectations
developed for pure diffusion processes, or the Wright-Fisher and Moran models from population genetics.  

This is not a completely new conjecture, of course.  We know mathematically that many different stochastic processes converge 
to a standard Fokker-Planck diffusion model, for example.  Within the archaeological literature, Mesoudi and Lycett (2009:42-43)
note that "perhaps some mix of conformity, anti-conformity and innovation combine to produce aggregate, population-level data 
that are indistinguishable from random copying.  However, to our knowledge, this claim has not yet been tested explicitly."

In addition to testing that conjecture explicitly, my goal is to map the "basin of attraction" -- i.e., what specific 
circumstances DO converge to apparent neutrality, and which retain statistical evidence of the original social 
learning rules.  

The simulation model contained here is written in Python, but relies on a number of modules which employ C code for 
performance.  Models are non-graphical, and are meant to run in batch mode with little or no console output, because a single 
run means little given the variability of these stochastic processes.  Data are logged to a MongoDB database,
and the result of experiments emerges from statistical analysis of the results.  (That was by way of explaining that
there aren't any GUI components, graphs, or little charts like one often sees in agent-based simulations).  

That said, here's how to run the simulation model yourself.  


## Getting Started ##

Download this repository to your system:
```Shell
git clone https://github.com/mmadsen/ctmixtures.git
```

The major dependencies are:

1.  MongoDB database server
1.  Python 2.7
1.  A list of Python modules

**WINDOWS USERS** should note that I have not tested this software on Microsoft Windows.  Please see my [notes on Windows](http://notebook.madsenlab.org/software.html), 
but there is no reason this will not work on Windows given appropriate tools installed.  


### MongoDB ###

MongoDB is a "no-SQL" database server that allows one to use very flexible database schemas, but is still searchable
and can have excellent performance when indexed for the type of queries being performed.  I use this, rather than 
logging raw data to disk, because it eliminates the need to post-process very large text files, and allows me to add
variables and derived statistics to the data set without (again) post-processing large text files.  The data format is JSON,
which is standard and well documented, so it is easy to operate on MongoDB from virtually every programming 
environment out there.  

You can install it easily on Linux via package managers.  For example, on Ubuntu:

```Shell
sudo apt-key adv --keyserver hkp://keyserver.ubuntu.com:80 --recv 7F0CEB10
echo 'deb http://downloads-distro.mongodb.org/repo/ubuntu-upstart dist 10gen' | sudo tee /etc/apt/sources.list.d/mongodb.list
sudo apt-get update
sudo apt-get install mongodb-org
```

On OS X, both the Homebrew and MacPorts package managers will install it complete with instructions to start it at boot, and there are step-by-step manual instructions at [mongodb.org](http://docs.mongodb.org/manual/tutorial/install-mongodb-on-os-x/).

Windows users should consult the instructions at the [MongoDB website](http://docs.mongodb.org/manual/tutorial/install-mongodb-on-windows/).


### Python 2.7 ###

You may already have Python 2.7 installed on your system.  This is especially be true on modern versions of Linux.  

On OS X, especially in contemporary versions (10.8, 10.9), you may have either 2.6 or 2.7, but I would strongly 
advise you to leave the system Python **completely alone** and use a local distribution instead.  The system itself relies strongly
on the system Python, so upgrading it or changing versions of packages may have unforeseen effects elsewhere on your system.  

Instead, the easiest solution is to install [Anaconda Python](https://store.continuum.io/cshop/anaconda/), a free Python distribution engineered to (a) run from a local directory on your system, giving you the ability to install any modules you need without conflict, and (b) pre-packaged with most of the key high-performance numerical and statistical packages for Python.  Some of the latter (e.g., SciPy) can be an extremely lengthy build process if you do it yourself.  Let Anaconda do it for you instead. To install Anaconda on Ubuntu:

```Shell
wget http://09c8d0b2229f813c1b93-c95ac804525aac4b6dba79b00b39d1d3.r79.cf1.rackcdn.com/Anaconda-2.0.1-Linux-x86_64.sh
bash Anaconda-2.0.1-Linux-x86_64.sh
```


### Other Required Modules or Dependencies ###

Once you have Anaconda or a similar Python distribution installed, there are a number of other modules that need 
to be installed.  The easiest way to do this is to use `pip` and let it run all the modules in the `requirements.txt` file in this repository. Do this:

```Shell
cd ctmixtures
pip install -r requirements.txt
```

You will also need the SWIG code generation system, which I used in the 
[slatkin-exact-tools](https://github.com/mmadsen/slatkin-exact-tools)
project to interface between python and C code. You can download these tools with `git clone https://github.com/mmadsen/slatkin-exact-tools.git`. The easiest way to install these is via the Linux package managers or Homebrew/MacPorts on OS X.  

After SWIG is installed, use the supplied script `install-slatkin-tools.sh` to compile and install the 
Slatkin Exact Tools. 


### Install CTMixtures Itself ###

You don't have to install this software anywhere, you can easily just adjust your PYTHONPATH to include "." (i.e., the 
current directory).  But you can also install the `ctmixture` modules in your local python distribution and use them from any 
directory by running `python setup.py install` in this directory.  


## Running a Simulation ##

The simplest driver program for running a simulation is in the `simulations` directory, and is called `sim-ctmixture-notimeaveraging.py`.  
This script runs a single simulation model, and is given both a configuration file (in JSON format), and several
command-line switches.  This simulation model runs a mixture of CT transmission rules for some period of time without 
sampling.  Given an interval to measure the survival of traits (ala Kandler and Steele), a sample is taken to start 
trait tracking, and then a sample is taken to complete the analysis.  In addition to completing the survival analysis, the
final sample taken also calculates and records a large number of sample- and population-specific statistics about 
trait frequencies, diversity, etc.  

Before running any simulation, ensure that MongoDB is running on your local machine (although it can be configured to 
run on another server), and that no username and password is required to access it (which is the default configuration).  
You may also want a good GUI tool for inspecting the database, since ALL output is sent there except for some
debugging tracing (if specified).  I strongly recommend using [Robomongo](http://robomongo.org/) on OS X and Linux.  

Configuration files specify a number of basic model elements, and are also used to construct large batches of 
simulation runs.  Examples occur in the `conf` directory, and a typical configuration looks like:

```
{
    "REPLICATIONS_PER_PARAM_SET" : 10,
    "POPULATION_SIZES_STUDIED" : [400],
    "INNOVATION_RATE" : [0.5, 1.0, 2.0, 5.0],
    "NUMBER_OF_DIMENSIONS_OR_FEATURES" : [1,2,4],
    "NUMBER_OF_TRAITS_PER_DIMENSION" : [4],
    "SAMPLE_SIZES_STUDIED" : [20,50,100],
    "SIMULATION_CUTOFF_TIME" : 3000000,
    "CONFORMISM_STRENGTH" : [0.0, 0.1, 0.2, 0.05, 0.3, 0.4, 0.5],
    "ANTICONFORMISM_STRENGTH" : [0.0, 0.1, 0.2, 0.05, 0.3, 0.4, 0.5],
    "INTERACTION_RULE_CLASS" : {"ctmixtures.rules.ConformistCopyingRule": 0.5, "ctmixtures.rules.AntiConformistCopyingRule": 0.5},
    "POPULATION_STRUCTURE_CLASS" : "ctmixtures.population.FixedTraitStructurePopulation",
    "INNOVATION_RULE_CLASS" : "ctmixtures.rules.InfiniteAllelesMutationRule",
    "NETWORK_FACTORY_CLASS" : "ctmixtures.population.SquareLatticeFactory",
    "TRAIT_FACTORY_CLASS" : "ctmixtures.traits.LocusAlleleTraitFactory"
}
```

This configuration specifies the names of several Python classes (from this repository) which are used to construct 
the simulation model.  This makes it easy to swap in a different class to construct different models (e.g., to 
put agents on a different spatial structure, for example).  

To examine mixtures of cultural transmission rules, the important configuration line here is the `INTERACTION_RULE_CLASS`.  
This is a JSON map or dictionary, with keys which are the Python classname of a class implementing a cultural 
transmission rule, and a proportion of the population which this rule should make up.  In the example above, we are
configuring a population with 50% conformist social learners, and 50% anti-conformists.  The simulation library assigns agents
rules at random, according to these proportions.  

A simulation program also takes a number of command line parameters, which specify the details of a single 
simulation run.  Running the script without any parameters gives a usage message which lists the parameters.  Currently 
(as of 6/15/14), these are:

```
mark:ctmixtures/ (master) $ simulations/sim-ctmixture-notimeaveraging.py                                                                                                                                         [7:14:43]
usage: sim-ctmixture-single.py [-h] --experiment EXPERIMENT [--debug DEBUG]
                               [--dbhost DBHOST] [--dbport DBPORT]
                               --configuration CONFIGURATION --popsize POPSIZE
                               --numloci NUMLOCI --maxinittraits MAXINITTRAITS
                               --conformismstrength CONFORMISMSTRENGTH
                               --anticonformismstrength ANTICONFORMISMSTRENGTH
                               --innovationrate INNOVATIONRATE --periodic
                               {1,0} [--kandlerinterval KANDLERINTERVAL]
                               [--simulationendtime SIMULATIONENDTIME]
```

A typical simulation command might then be:

```
mark:ctmixtures/ (master*) $  simulations/sim-ctmixture-notimeaveraging.py --experiment foo --configuration conf/conformism-mixture.json  --simulationendtime 200000  --kandlerinterval 100 --numloci 2 --maxinittraits 5  --periodic 0  --conformismstrength 0.1 --anticonformismstrength 0.1  --innovationrate 1.0 --debug 1 --popsize 100 

```

Since debugging output is turned on (`--debug 1`), the output will resemble this:

```
2014-06-16 18:06:09,707 DEBUG: experiment name: foo
2014-06-16 18:06:09,708 DEBUG: configured theta = 1.0, using numloci 2 * per-locus mutation rate 0.00502512562814 = all-loci innovation rate: 0.0100502512563
2014-06-16 18:06:09,708 DEBUG: Taking a Kandler trait survival sample of 10000 timesteps, beginning at tick 190000
Couldn't import dot_parser, loading of dot files will not be possible.
2014-06-16 18:06:10,202 DEBUG: Configuring CT Mixture Model with structure class: ctmixtures.population.FixedTraitStructurePopulation graph factory: ctmixtures.population.SquareLatticeFactory interaction rule: {u'ctmixtures.rules.ConformistCopyingRule': 0.5, u'ctmixtures.rules.AntiConformistCopyingRule': 0.5}
2014-06-16 18:06:10,213 DEBUG: Lattice model:  popsize 100, lattice will be 10.0 by 10.0
2014-06-16 18:06:10,215 DEBUG: Conformist rule operating at strength: 0.1
2014-06-16 18:06:10,215 DEBUG: Anticonformist rule operating at strength: 0.1
2014-06-16 18:06:10,215 DEBUG: creating 50 obj for rule ctmixtures.rules.ConformistCopyingRule
2014-06-16 18:06:10,215 DEBUG: creating 50 obj for rule ctmixtures.rules.AntiConformistCopyingRule
2014-06-16 18:06:10,216 INFO: Starting urn:uuid:cc021fce-53bc-4aea-b66c-0258e1436792
2014-06-16 18:06:11,287 DEBUG: time: 100000  copies by locus: [49947, 50053]  innovations: 1047 innov by locus: [510, 537]
2014-06-16 18:06:12,263 DEBUG: Starting Kandler remaining trait tracking at 190000
2014-06-16 18:06:12,377 DEBUG: time: 200000  copies by locus: [100002, 99998]  innovations: 2044 innov by locus: [1042, 1002]
2014-06-16 18:06:12,377 DEBUG: Stopping Kandler remaining trait tracking at 200000
2014-06-16 18:06:12,384 DEBUG: locus: 0 snapshot one: set([720, 952, 994, 680, 983])  snapshot two: set([680, 1005, 1006, 1039, 1043, 1013, 1046, 1047, 952]) 
2014-06-16 18:06:12,384 DEBUG: locus: 0 intersection: set([952, 680])
2014-06-16 18:06:12,384 DEBUG: locus: 1 snapshot one: set([736, 937, 941, 943, 944, 913, 914])  snapshot two: set([736, 992, 1003, 1007, 913, 981]) 
2014-06-16 18:06:12,384 DEBUG: locus: 1 intersection: set([736, 913])
2014-06-16 18:06:13,118 INFO: Completed: urn:uuid:cc021fce-53bc-4aea-b66c-0258e1436792  Elapsed: 3.41031599045
```

The final line, which is printed even without debugging output, gives a guaranteed globally unique identifier for this simulation run,
across any number of systems in a cluster, etc.  All records in the database are identified with this "UUID" and thus we can 
run simulations across many different systems, merge their databases, and still have unique records of experiments.  

If you open Robomongo, and select a connection with `127.0.0.1` (your local machine), and drill down in the left panel, you will see
a database automatically created for the experiment (named "foo"), which contains the raw output of the simulation (`foo_samples_raw`).  
Expand this database, and you will see two "collections" (i.e., tables).  One contains timing data for simulation runs, which aids in 
predicting how long large batches will take and improving simulation performance.  

The other table contains the simulation output, and is named `mixture_model_stats`.  Expand this, and you will see a row of data.  In the 
default mode ("tree mode"), you can expand variables which have an arrow to their left, and you can see that many statistics are 
captured for the simulation run.  These are a mixture of the original configuration parameters (which are ALWAYS kept with the data), 
and observable values (e.g., `num_trait_configurations` records the number of combinations of the two dimensions (i.e., trait combinations from two 
different dimensions or loci) present in the population when sampled).  Some variables record arrays of values, or dictionaries with keys and values.  

I'll write more about what each statistic means, but this should get you running with the model, if you have a row of data in MongoDB with about 37
main entries (at the time this document was updated).  

## Simulating CT Mixtures with Time Averaging ##

TBD.



## Running Batches of Models ##

TBD.  




