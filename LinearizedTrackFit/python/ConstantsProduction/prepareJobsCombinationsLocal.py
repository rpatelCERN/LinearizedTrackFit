from prepareJobCombinationsLocalPreEstimates import *
from prepareCombinationsLocalSecondEstimates import *
from jobMonitoring import *

__author__ = 'demattia'


prepare_pre_estimates()

submit_and_monitor(jobs_queue, maximum_parallel_jobs=4, monitor_interval=10)
