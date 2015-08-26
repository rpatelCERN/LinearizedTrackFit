from prepareJobCombinationsLocalPreEstimates import *
from prepareCombinationsLocalSecondEstimates import *
from prepareCombinationsLocalLTF import *
from jobMonitoring import *

__author__ = 'demattia'


# prepare_pre_estimates()
# prepare_second_estimates()
prepare_linearized_track_fit()
submit_and_monitor(jobs_queue, maximum_parallel_jobs=4, monitor_interval=10)
