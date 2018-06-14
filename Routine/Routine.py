#!/usr/bin/env python
import httplib2
import urllib2
import os
import json
from ion.plugin import *
from subprocess import *
import time

class Routine(IonPlugin):
	""" Plugin object to run all default plugins for routine analysis."""
	version = '0.1'
	json_dat = {}
			
	def launch(self):
                print "START Routine plugin"
                
                try:
                        with open('startplugin.json', 'r') as j:
                                self.json_dat = json.load(j)
                except:
                        print 'Error reading plugin json.'

                pk = self.json_dat['runinfo']['pk']

                print "run_id  %s" % pk

                ## TN report FIX 09/02/2017                                                                 
                resultsName = self.json_dat['expmeta']['results_name']
                if resultsName.endswith("_tn"): # c'est le rapport TN, inutile de lancer Routine dessus     
                        print "--- TN REPORT ---"
                        print " [%s] Waiting for full report to be completed..." % (time.strftime("%H:%M:%S"))
                        goodPK = int(pk)-1
                        fullRunCompleted = False
                        while not fullRunCompleted: # check when full report is completed, then re-launch   
                                time.sleep(60) #600                                                         
                                api_url = self.json_dat['runinfo']['api_url'] + '/v1/results/' + str(goodPK)
                                f = urllib2.urlopen(api_url)
                                d = json.loads(f.read())
                                state = d['status']
                                if state == 'Completed':
                                        fullRunCompleted = True
                        time.sleep(600) # wait 10min                                                        
                        print " [%s] %s FULL REPORT is now ready" % (time.strftime("%H:%M:%S"),d['resultsName'])
                        
                        h = httplib2.Http()
                        h.add_credentials('ionadmin', 'ionadmin')
                        headers = {'Content-type': 'application/json','Accept': 'application/json'}
                        url = self.json_dat['runinfo']['net_location'] + '/rundb/api/v1/results/' + str(goodPK) + '/plugin/'
                        pluginName = 'Routine'
                        pluginUpdate  = {'plugin': [pluginName]}
                        resp, content = h.request(url, 'POST', body=json.dumps(pluginUpdate),headers=headers)
                        print resp
                        print content
                        print " [%s] Routine plugin launched on full report, exiting." % (time.strftime("%H:%M:%S"))

                else:
                        htmlOut = open('plugin_log.html', 'w')
                        htmlOut.write('''                                                                   
                        <div><object width="1280" height="800" data="plugin_log.txt"></object></div>        
                        <div><object width="1280" height="800" data="error_log.txt"></object></div>         
                        ''')
                        ### PLUGIN SLOT FIX : CHECK if no other Routine plugin is running (SGE plugins slots limit) #################                                                                                  
                        cmd = Popen(['pgrep','-f','run_routine.py'], stdout=PIPE)
                        pid_found = cmd.communicate()[0]
                        if pid_found:
                                print " [%s] WARNING : ROUTINE plugin is already Started on another Run" % (time.strftime("%H:%M:%S"))
                                print " [%s] launching wait-retry.py..." % (time.strftime("%H:%M:%S"))
                                # lancer script wait_retry.py FIX lancement dans un shell separe sinon quand ce script se termine, wait_retry.py se termine aussi                                                      
                                cmd = Popen(['python','/results/plugins/Routine/wait_retry.py',str(pk)],stdout=open('plugin_log.txt','w'),stderr=open('error_log.txt','w'),preexec_fn=os.setpgrp)
                                print " [%s] ...done." % (time.strftime("%H:%M:%S"))
                                print " [%s] This plugin will now terminate. Then Routine will be automatically re-launched when no other Routine is Started" % (time.strftime("%H:%M:%S"))
                                return False
                        else:
                                print "No other instance of Routine found"
                                print " [%s] Launching run_routine.py ..." % (time.strftime("%H:%M:%S"))
                                # lancer script dans un shell separe                                        
                                cmd = Popen(['python','/results/plugins/Routine/run_routine.py',self.json_dat['runinfo']['results_dir']],stdout=open('plugin_log.txt','w'),stderr=open('error_log.txt','w'),preexec_fn=os.setpgrp)
                                print " [%s] ...done." % (time.strftime("%H:%M:%S"))

                return True
		
if __name__ == '__main__':
	PluginCLI(Routine())
