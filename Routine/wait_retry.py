#!/usr/bin/python
import sys
import httplib2
import urllib2
import json
import time
from subprocess import *

### OBSOLETE (Du to the SGE 2 slots queue limit for plugins, Routine can't be launched more than 1 time simultaneously.)
### OBSOLETE (The purpose of this script is to check if no other Routine plugin is running, and then re-launch the Routine plugin)

## NEW 20/01/2017 : On S5 sequencer, SGE slot limit is now 1. Routine plugin cannot be launched whithout blocking plugins analysis pipeline...
## Routine is now performed by run_routine.py
## This script is still usefull to force one run_routine.py at a time, in order to let a run to finish all plugins before beginning another run.

### USE : python wait_retry.py PK 
pk = sys.argv[1]

# log de ce script
log = open('/results/plugins/Routine/log.txt', 'a+')
now = time.strftime('%c')
#print "Current date & time " + time.strftime('%c') + "\n"
log.write("Current date & time " + time.strftime('%c') + "\n")

h = httplib2.Http()
h.add_credentials('ionadmin', 'ionadmin')
headers = {'Content-type': 'application/json','Accept': 'application/json'}
url = 'http://129.21.10.5' + '/rundb/api/v1/results/' + str(pk) + '/plugin/'

# wait : 30 min and retry
maxtry = 30
started = True
while started:
        #time.sleep(120)
	time.sleep(1800)
	cmd = Popen(['pgrep','-f','run_routine.py'], stdout=PIPE)
	pid_found = cmd.communicate()[0]
	if pid_found:
		started = True
		#print time.strftime('%H:%M:%S') + " : run_routine.py still for another run...\n"
                log.write(time.strftime('%H:%M:%S') + " : run_routine.py still for another run...\n")
	else:
		started = False
	maxtry = maxtry - 1
	if maxtry == 0: # arret de ce script apres x essais.
		exit(1)

if not started:
        log.write("...run_routine.py finished!\n")
	#print "...run_routine.py finished!\n"
	pluginName = 'Routine'
	pluginUpdate  = {'plugin': [pluginName]}
        log.write("-> Re-Launching Routine...\n")
	#print "-> Re-Launching Routine...\n"
	resp, content = h.request(url, 'POST', body=json.dumps(pluginUpdate),headers=headers )
	log.write("-> Routine launched at " + time.strftime('%H:%M:%S') + "\n")
        #print "-> Routine launched at " + time.strftime('%H:%M:%S') + "\n"
log.close()	
exit(0)
