#!/usr/bin/python

import os

#das_client.py --query "run dataset=/DoubleMuon/Run2016C-PromptReco-v2/AOD | grep run.run_number"

def getRunsInDataset(dataset):
    dbs = '''das_client.py --query "run dataset='''+dataset+''' | grep run.run_number" --limit 999'''
    dbsOut = os.popen(dbs)
    files = []
    for line in dbsOut:
        line = line.rstrip()
        try:
            files.append(int(line))
        except ValueError:
            continue
    return files

def getFillFromRun(run):
    dbs = '''das_client.py --query "run run=%i | grep run.lhcFill"''' % (run)
    dbsOut = os.popen(dbs)
    fill = 0
    for line in dbsOut:
        try:
            fill = int(line)
        except ValueError:
            continue
    return fill
    

if __name__=='__main__':
    dataset = '/DoubleEG/Run2016B-23Sep2016-v3/AOD'
    #dataset = '/DoubleEG/Run2016C-23Sep2016-v1/AOD'
    #dataset = '/DoubleEG/Run2016G-23Sep2016-v1/AOD'
    runs = getRunsInDataset(dataset)
    print 'Dataset', dataset, 'contains runs:', runs
    out = open(dataset.replace('/', '_'), 'w')
    for run in runs:
        print '%i:%i' % (run, getFillFromRun(run))
        out.write('%i:%i\n' % (run, getFillFromRun(run)))
        
