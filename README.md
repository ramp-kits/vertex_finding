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
It will then get download as compressed files and unpacked to data/train data/test (might take a while)
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




