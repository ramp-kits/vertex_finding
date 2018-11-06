# RAMP starting kit on Primary Vertex reconstruction

#### Set up

Open a terminal and

1. install the `ramp-workflow` library (if not already done)
  ```
  $ pip install git+https://github.com/paris-saclay-cds/ramp-workflow.git
  ```
  
2. Follow the ramp-kits instructions from the [wiki](https://github.com/paris-saclay-cds/ramp-workflow/wiki/Getting-started-with-a-ramp-kit) (i.e. set-up the environment with conda or pip)


### Test the kit
After setting up the environment,
```
ramp_test_submission --quick-test
```
should run the baseline PV finding and print out the scores.
#### Local notebook

Get started on this RAMP with the [dedicated notebook](pv_finding.ipynb).



#### Help
Go to the `ramp-workflow` [wiki](https://github.com/paris-saclay-cds/ramp-workflow/wiki) for more help on the [RAMP](http:www.ramp.studio) ecosystem.




