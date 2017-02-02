#!/usr/bin/env python
import sys, os, subprocess

try:
    card_path = sys.argv[1]
except:
    print('[card]')
    exit()

if not os.path.exists(card_path):
    print( card_path + ' does not exist.')
    exit()

ananame = ''
runname = ''
option_mode = False
options = []
for line in open(card_path):
    try:
        if line.split()[0] == 'Analysis:':
            ananame = line.split()[1]
    except: pass

    try:
        if line.split()[0] == 'Run' and line.split()[1] == 'name:':
            runname = line.split()[2]
    except: pass

    try:
        if line.split()[0] == 'Eventfile:':
            eventFile = line.split()[1]
            if not os.path.exists(eventFile):
                print(eventFile + ' does not exist')
                exit()            
    except: pass

    if 'Options start' in line: 
        option_mode = True
        continue
    if 'Options end' in line: break
    if option_mode: options.append(line)

if ananame == '':
    print('"Analysis:" was not found in the card.')
    exit()

if runname == '':
    print('"Run name:" was not found in the card.')
    exit()

if eventFile == '':
    print('"EventFile:" was not found in the card.')
    exit()

script = runname + '.script'
scriptFile = open( script, 'w')
scriptFile.write( 'add Analysis ' + ananame + '\n')
scriptFile.write( 'add InputFile ' + eventFile  + '\n')
scriptFile.write( 'add OutputFile ' + runname + '.yaml'  + '\n')
for op in options: scriptFile.write( op )
scriptFile.write( 'launch'  + '\n')
scriptFile.close()

subprocess.call(['cat',script])
subprocess.call(['atom','-b',script])













