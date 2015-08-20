import os
import subprocess
import time

__author__ = 'demattia'


def check_pid(pid):
    """ Check For the existence of a unix pid. """
    try:
        # Kill with value 0 does not kill the process
        os.kill(pid, 0)
        return True
    except OSError:
        return False


def submit_and_monitor(jobs_queue, maximum_parallel_jobs, monitor_interval):
    """
    Now submit the jobs to keep the total number of jobs running in parallel < N.
    """
    print jobs_queue
    # maximum_parallel_jobs = 4

    failed_processes = []
    running_processes = []

    while len(jobs_queue) > 0 or len(running_processes) > 0:
        print "jobs left", len(jobs_queue)

        # Get the list of all running process ids
        failed_processes.extend([p for p in running_processes if p[1].poll() is not None and p[1].poll() != 0])
        running_processes = [p for p in running_processes if p[1].poll() is None]
        print "running jobs:", [p[0] for p in running_processes]

        for i in range(min(maximum_parallel_jobs - len(running_processes), len(jobs_queue))):
            job_command = jobs_queue.pop()
            print job_command
            # job_command = ""
            process = subprocess.Popen(job_command,
                                       stdout=subprocess.PIPE,
                                       stderr=subprocess.STDOUT)
            running_processes.append([job_command.lstrip(os.getcwd()).rstrip("jobFile.sh"), process])

        if len(jobs_queue) == 0 and len(running_processes) == 0:
            break
        # Every 10 seconds
        time.sleep(monitor_interval)

    print "All jobs finished"
    if len(failed_processes) == 0:
        print "No failed jobs"
    else:
        print "Failed jobs:"
        for p in failed_processes:
            print p[0]
