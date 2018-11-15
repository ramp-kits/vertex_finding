# RAMP starting kit on Primary Vertex reconstruction

[![Build Status](https://travis-ci.org/ramp-kits/vertex_finding.svg?branch=master)](https://travis-ci.org/ramp-kits/vertex_finding)

### Set up

1. clone this repository
  ```
  git clone https://github.com/ramp-kits/vertex_finding.git
  cd vertex_finding
  ```

2. install the dependancies
  - with [conda](https://conda.io/miniconda.html)
  ```
  conda update conda                   # make sure conda is up-to-date
  conda env create -f environment.yml  # use environment.yml to create the 'vertex_finding' env
  source activate vertex_finding       # activates the virtual env
  ```
  - without `conda` (best to use a **virtual environment**)
  ```
  python -m pip install -r requirements.txt
  ```

3. download the data
  ```
  python download_data.py
  ```
  After download, the data will be unpacked to `data/train` and `data/test` (might take a while). By default, this will extract 5k train and 1k test events.

### Building the baseline solution (using C++)

The baseline solution was adapted from the Primary Vertex reconstruction used in LHCb.   
It requires `CMake`, the `Boost` library (both installed with `conda`) and a decent C++ compiler to be built.

- manual build
    ```bash
    mkdir build && cd build
    cmake .. && make
    cd ..
    ```
- automated build, when calling the baseline submission for the first time
    ```bash
    ramp_test_submission --quick-test --submission baseline
    ``` 

### Testing the kit

After setting up the environment, run the starting kit (random values) and the baseline solutions.
```bash
ramp_test_submission --quick-test
ramp_test_submission --quick-test --submission baseline
```
They should run on a subset of the data and print out the scores.

To process all of the locally available data, remove the `--quick-test` flag.
```
ramp_test_submission --submission baseline
```

### Local notebook

Get started on this RAMP with the [dedicated notebook](vertex_finding_starting_kit.ipynb).

### Help
Go to the `ramp-workflow` [wiki](https://github.com/paris-saclay-cds/ramp-workflow/wiki) for more help on the [RAMP](http:www.ramp.studio) ecosystem.
