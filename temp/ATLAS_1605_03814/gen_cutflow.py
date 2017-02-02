#!/usr/bin/env python

import sys, os, copy, subprocess
from math import sqrt
from atom_validation import *
from tex_formatting import *

###########################################################
if len(sys.argv) == 1:
    print('[runname] or [runname] [atom.yaml]')

runname = sys.argv[1]
card_path = runname + '.card'
try:
    atomFile = sys.argv[2]
except:
    atomFile = runname + '.yaml'
    
if not os.path.exists(card_path):
    print( card_path + ' does not exist.')
    exit()

if not os.path.exists(atomFile):
    print( atomFile + ' does not exist.')
    exit()

###########################################################

ananame, mode = '', ''
Nexp, Nev0 = '', '';
discript = {}
cutList = []
for line in open(card_path):

    if 'Analysis:' in line: ananame = line.split()[1]
    if 'Run name:' in line: 
        runname2 = line.split()[2]
        if runname2 != runname:
            print runname2, runname
            print 'inconsistent runname!!  Check the Run name in the card'
            exit()
    if line[0] != '#' and 'Num. Events Exp:' in line: Nexp = float(line.split('Num. Events Exp:')[1].split()[0] )
    if line[0] != '#' and 'Process:' in line: discript['Process:'] = line.split('Process:')[1]
    if line[0] != '#' and 'Generator:' in line: discript['Generator:'] = line.split('Generator:')[1]
    if line[0] != '#' and 'Parameters:' in line: discript['Parameters:'] = line.split('Parameters:')[1]
    if line[0] != '#' and 'Num. NoCut Events Exp:' in line: Nev0 = float(line.split('Num. NoCut Events Exp:')[1].split()[0] )
 

    if 'start of cutflow table' in line:
        mode = 'cutflow'        
        continue
    if 'Backup Lists' in line: break

    if mode == 'cutflow':
        dict = {}
        try:
            dict['Title'] = line.split("'")[1]
        except: continue
        dict['Atom Name'] = line.split("'")[3]
        result = line.split("'")[4]
        if len(line.split("'")) == 7: 
            parent_name = line.split("'")[5]
            dict['Parent Atom Name'] = parent_name 
        if len(line.split("'")) > 7:
            print 'ERROR in the following line:'
            print line
            print "too many ' "
            exit()
        try:
            dict['Cumulative'] = float(result)
        except: 
            if len( result.split('(') ) > 1:
                dict['Cumulative'] = float(result.split('(')[0])
                dict['Cum. Error'] = float( result.split('(')[1].split(')')[0] )
            else:
                pass
        cutList.append(dict)            

if ananame == '':
    print('"Analysis:" was not found in the card.')
    exit()

if runname == '':
    print('"Run name:" was not found in the card.')
    exit()

###########################################################

for line in open(atomFile):
    if "Processed Events:" in line:
        discript['Number of Atom MC events:'] = line.split('Processed Events:')[1]

disc_lines = []
disc_lines.append('\\begin{itemize}')
disc_lines.append('\\item  Process: ' + discript['Process:'])
disc_lines.append('\\item  Parameters: ' + discript['Parameters:'])
disc_lines.append('\\item  Number of Atom MC events: ' + str(discript['Number of Atom MC events:']) ) 
disc_lines.append('\\item  Event Generator: { \\tt ' + discript['Generator:'] + ' }' ) 
disc_lines.append('\\end{itemize}')

###########################################################

perc = 100.
prev, prev_err = 1, 0
cutList0 = copy.deepcopy(cutList)
cutList = []
Nnorm = 1.

if Nev0 != '':
    print 'Event Mode: Normalized by Num. NoCut Events =', Nev0 
    Nnorm = Nev0
for cut0 in cutList0:
    cut = copy.deepcopy(cut0)
    cut['ProcIds'] = [0]    
    eff = cut['Cumulative'] / Nnorm    
    cut['Cumulative'] = eff * perc
    if 'Cum. Error' in cut.keys():
        err = cut['Cum. Error'] / Nnorm
    else:
        err = sqrt(eff)/sqrt(Nexp)    
    cut['Cum. Error'] = err * perc
    if 'Parent Atom Name' in cut.keys():
        found = False
        for i in range(len(cutList)):
            dm = cutList[i]
            if dm['Atom Name'] == cut['Parent Atom Name']:
                prev = dm['Cumulative'] / perc
                prev_err = dm['Cum. Error']
                cut['Parent ID'] = i
                found = True
        if not found:
            print 'Error at', cut['Atom Name']
            print 'Parent Atom Name = ["' + cut['Parent Atom Name'] + '"] does not exist.'
            print 'Abort'
            exit()
    if prev > 0: 
        cut['Individual'] = eff/prev * perc
    else:
        cut['Individual'] = 0.        
    cut['Indiv. Error'] = (sqrt(err) + sqrt(prev_err))**2 * perc  
    cutList.append(cut)
    ##
    prev = eff   
    prev_err = err

#for cut in cutList:
#    for key, val in cut.items():        
#        print key, '    ', val

###########################################################

cutflow = {}
cutflow['Caption'] = ''
cutflow['Comments'] = ''
cutflow['CutFile'] = atomFile
cutflow['Cuts'] = cutList

val = Validation()

texoutput, screenoutput = val.processCutFlow(ananame, runname, cutflow) # read Atom output (runname.yaml)
print screenoutput

###########################################################
###########################################################

tex = tex_format()

pwd = os.path.dirname(os.path.realpath(__file__))
texDir = os.path.join(pwd, 'tex')
texFile = os.path.join(texDir, runname)
fout = open(texFile+'.tex', 'w')

fout.write( tex.header )

fout.write( '\\begin{document}' + '\n' )
secname = ''
for x in runname.split('_'): secname += x + ' ' 
fout.write( '\\section{' + secname + '}' + '\n')

fout.write('\n')

for line in disc_lines: fout.write(line + '\n')
fout.write('\n')

for line in texoutput:
    if 'newline' in line: continue
    fout.write( line + '\n')
fout.write( tex.end_document )

###################################

exit()







