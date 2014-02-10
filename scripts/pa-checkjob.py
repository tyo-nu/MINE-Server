import argparse
import sys
import os
import time
import traceback
from biokbase.probabilistic_annotation.Helpers import job_info_dict
from biokbase.probabilistic_annotation.Client import _read_inifile
from biokbase.userandjobstate.client import UserAndJobState, ServerError as JobStateServerError

desc1 = '''
NAME
      pa-checkjob -- check status of probabilistic annotation jobs

SYNOPSIS      
'''

desc2 = '''
DESCRIPTION
      Check status of probabilistic annotation jobs submitted by the user.  For
      each job, information about the job is displayed.  A job that has completed
      is then deleted from the system.

      The --jobID optional argument is the identifier of a specific job to
      check.

      The ujs-url optional argument specifies an alternate URL for the user and
      job state service.
'''

desc3 = '''
EXAMPLES
      Check all jobs:
      > pa-checkjob
      
      Check a specific job:
      > pa-checkjob --jobID 52b317cbe4b0ef8357331c59

SEE ALSO
      pa-annotate

AUTHORS
      Matt Benedict, Mike Mundy 
'''

def print_job(info):
    # Check if the job had an error.
    if info['error']:
        print "Job '%s' (%s) ended with error '%s'." %(info['id'], info['description'], info['status'])
        print 'Error details:'
        print ujsClient.get_detailed_error(info['id'])
        ujsClient.delete_job(info['id'])

    # Check if the job is complete.
    elif info['complete']:
        print "Job '%s' (%s) completed successfully." %(info['id'], info['description'])
        ujsClient.delete_job(info['id'])
    
    # Job is still running.
    else:
        print "Job '%s' (%s) has status '%s' and is working on task %s of %s.  Check again later." \
            %(info['id'], info['description'], info['status'], info['total_progress'], info['max_progress'])
            
    return
    
if __name__ == "__main__":
    # Parse options.
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, prog='pa_checkjob', epilog=desc3)
    parser.add_argument('-j', '--jobID', help='job ID', action='store', dest='jobID', default=None)
    parser.add_argument('--ujs-url', help='url for user and job state service', action='store', dest='ujsURL', default='https://kbase.us/services/userandjobstate/')
    usage = parser.format_usage()
    parser.description = desc1 + '      ' + usage + desc2
    parser.usage = argparse.SUPPRESS
    args = parser.parse_args()
    
    # Get the list of jobs for the user.
    auth = _read_inifile()
    ujsClient = UserAndJobState(args.ujsURL)
    try:
        jobList = ujsClient.list_jobs([ auth['user_id'] ], 'RCE')
    except JobStateServerError as e:
        print e.message
        exit(1)
    
    # See if the user has any jobs in the list.
    if len(jobList) == 0:
        print 'There are no jobs for you.'
        exit(1)

    # Print info about the specific job if requested.
    if args.jobID is not None:
        for job in jobList:
            info = job_info_dict(job)
            if args.jobID == info['id']:
                print_job(info)
                exit(0)
        print "Job '%s' was not found." %(args.jobID)
        exit(1)
        
    # Print all of the jobs in the list.
    for job in jobList:
        print_job(job_info_dict(job))
            
    exit(0)