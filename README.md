# RAMP starting kit on Primary Vertex reconstruction

#### Set up

Open a terminal and

1. install the `ramp-workflow` library (if not already done)
  ```
  $ pip install git+https://github.com/paris-saclay-cds/ramp-workflow.git
  ```
  
2. Follow the ramp-kits instructions from the [wiki](https://github.com/paris-saclay-cds/ramp-workflow/wiki/Getting-started-with-a-ramp-kit) (i.e. set-up the environment with conda or pip)


### Downloading the data

One can get the data with
```
python download.py
```
It will then get download as compressed files and unpacked to data/train data/test (might take a while). By default, this will extract 5k train and 1k test events.

### Building the baseline solution
The baseline solution was adapted from the Primary Vertex reconstruction used in LHCb. It can be compiled in the following way:
```
mkdir build
cd build
cmake ..
make
cd ..
```
You will then need to copy the resulting library files (.so, .dylib, ... depending on your architecture) from the build/PatPV and build/baseline directories into the submissions/starting_kit/ directory.

### Testing the kit
After setting up the environment,
```
ramp_test_submission --quick-test
```
should run the baseline PV finding and print out the scores, running on a subset of the total data.
Just doing 
```
ramp_test_submission
```
will use all locally available data.

#### Local notebook

Get started on this RAMP with the [dedicated notebook](vertex_finding_starting_kit.ipynb).



#### Help
Go to the `ramp-workflow` [wiki](https://github.com/paris-saclay-cds/ramp-workflow/wiki) for more help on the [RAMP](http:www.ramp.studio) ecosystem.




