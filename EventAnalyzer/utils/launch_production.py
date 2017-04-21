#!/usr/bin/python
from sys import argv, exit
import LaunchOnFarm
import commands
from datetime import datetime

def createFilesListFromEOSPath(path, files):
    if path[-1]!='/':
        path += '/'
    files = files.split('\n')
    out_filename = 'tmp_filelist_'+datetime.now().strftime('%y%m%d-%H%M')
    out_content = ''
    out_file = open(out_filename, 'w')
    for f in files:
        out_file.write("'"+path+f+"',\n")
    return out_filename

def main(argv):
    command_out = commands.getstatusoutput('eos ls '+argv[1]+' | grep root')
    if command_out[0]!=0:
        print 'Directory \''+argv[1]+'\' not found on /eos/cms!'
        return -1

    files_list = createFilesListFromEOSPath(argv[1], command_out[1])
    print files_list, 'created and populated!'

    LaunchOnFarm.useLSF = True
    #LaunchOnFarm.Jobs_Queue = '1nd'
    LaunchOnFarm.Jobs_Queue = '8nh'
    LaunchOnFarm.SendCluster_Create(argv[3]+'/', argv[3])
    num_jobs = LaunchOnFarm.SendCluster_LoadInputFiles(files_list, 500)
    for i in range(num_jobs):
        LaunchOnFarm.SendCluster_Push(['CMSSW', argv[2]])
    LaunchOnFarm.SendCluster_Submit()

    return 0

if __name__=='__main__':
    if len(argv)<4:
        print 'Usage:', argv[0], ' [list of input files] [python configuration] [request name]'
        exit(0)
    main(argv)
