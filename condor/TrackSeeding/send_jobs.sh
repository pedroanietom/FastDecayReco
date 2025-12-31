#!/bin/bash

for run in $(head -n 42 42_tracking_list.txt); do
  condor_submit -append "RunNumber=$run" condor.job
done

