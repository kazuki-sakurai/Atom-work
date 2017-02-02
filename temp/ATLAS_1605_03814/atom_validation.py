import yaml, atom, tex_formatting, out_formatting
import os, subprocess, shutil, re
import math


class Validation:
    def __init__(self):
        self.download_script = 'download.sh'
        self.running_script = 'job.sh'
        self.work_path = './work/'
        self.output_path = './output/'
        self.options_path = ''
        self.save_output = True
        self.tformatter = tex_formatting.tex_format()
        self.oformatter = out_formatting.output_format()
        self.atom_numevents = 0

    def getAnalysesList(self, paths, analyses_list = []):
        analist = {}
        skipcheck=True
        tmp_analyses_names = [ x.split(':')[0] for x in analyses_list ]
        tmp_analyses_runs = {}
        [ tmp_analyses_runs.update({y[0]:y[1]}) if len(y) > 1 else tmp_analyses_runs.update({y[0] : ''}) for y in [ x.split(':') for x in analyses_list ] ]
        if len(analyses_list) > 0:
            skipcheck = False
        for item in paths:
            for root, subFolders, files in os.walk(item):
                for item in files:
                    ext = item[-5:]
                    if ext.lower() != '.yaml':
                        continue
                    infile = open(os.path.join(root,item), 'r')
                    header = infile.readline()
                    infile.close()
                    if header.startswith('#') and header.find('validation card') != -1:
                        ana_name = item[:-5] #strip .yaml extension
                        if skipcheck or ana_name in tmp_analyses_names:
                            analist.update({ ana_name : [ os.path.join(root,item), tmp_analyses_runs[ana_name] ] })
        return analist

    def writeMainLaTeX(self, analysis_list):
        with open(self.output_path + '/validation.tex', 'w') as fout:
            fout.write(self.tformatter.begin_document)
            fout.write('\n')
            for item in analysis_list:
                secname = item
                underscore_positions = [n for n in xrange(len(item)) if item.find('_', n) == n]
                if len(underscore_positions) == 2:
                    #it's EXP_xxxx_yyyy -> 'EXP xxxx.yyyy'
                    #or EXP_hepex_zzzzz -> EXP hep-ex/zzzz'
                    if secname[underscore_positions[0]+1].upper() == 'H':
                        secname = secname[:underscore_positions[0]] + ' ' + \
                            secname[underscore_positions[0]+1:underscore_positions[1]].lower() + '/' + \
                            secname[underscore_positions[1]+1:]
                    else:
                        tmplist = list(secname)
                        tmplist[underscore_positions[0]] = ' '
                        tmplist[underscore_positions[1]] = '.'
                        secname = ''.join(tmplist)
                elif len(underscore_positions) >= 3:
                    #it's EXP_CONF_xxxx_yyy -> 'EXP-CONF-xxxx-yyy'
                    #it's EXP_PAS_SUS_xx_yyy -> 'EXP-PAS-SUS-xx-yyy'
                    tmplist = list(secname)
                    for pos in underscore_positions:
                        tmplist[pos] = '-'
                    secname = ''.join(tmplist)
                fout.write('\\section{' + secname + '}\n')
                fout.write('\\input{' + item + '}\n')
                fout.write('\n')
            fout.write('\n')
            fout.write(self.tformatter.end_document)

    def processYaml(self, file_in, analysis_name, run_name):
        outputTeX = []
        outputScreen = ''
        source = yaml.load(file_in)
        repo = ''
        print source.keys()
        if 'Repository' in source.keys():
            repo = source['Repository']
        if 'Runs' in source.keys():
            for idx, item in enumerate(source['Runs']):
                try:
                    name = item['Run Name']
                    evfiles = item['Event Files']
                    optfile = item['Atom Options']
                    descr = item['Description']
                except:
                    print "Error in '" + file_in.name + "': mandatory fields not found."
                    continue
                if len(run_name) > 0 and (name != run_name) and ('"'+name+'"' != run_name) and ("'"+name+"'" != run_name):
                    continue
                run_path = 'run_' + ('00'+str(idx+1))[-2:]
                outfiles = []
                if 'Plot' in item.keys():
                    if 'PlotFile' in item['Plot'].keys():
                        outfiles.append(item['Plot']['PlotFile'])
                if 'CutFlow' in item.keys():
                    [ outfiles.append(x['CutFile']) for x in item['CutFlow'] if 'CutFile' in x.keys() ]
                if len(outfiles) == 0:
                    print "Error in '" + file_in.name + "': no output files found."
                cwdpath = self.initRun(analysis_name, run_path, evfiles, optfile, list(set(outfiles)), repo)
                print 'cwdpath = ', cwdpath
                #self.runAtom(cwdpath)
                tmptex, tmpscreen = self.initSection(name,descr)
                outputTeX.extend(tmptex)
                outputScreen += tmpscreen
                for kind, info in item.iteritems():
                    if kind == 'Plot':
                        if 'PlotFile' in info.keys():
                            pname = info['PlotFile']
                            if 'Plots' in info.keys():
                                for pindex, plotinfo in enumerate(info['Plots']):
                                    tmptex, tmpscreen = self.processCmpPlot(cwdpath, analysis_name, name, idx, pname, plotinfo, pindex)
                                    outputTeX.extend(tmptex)
                                    outputScreen += tmpscreen
                            else:
                                print "Error! Cannot find list of plots"
                        else:
                            print "Error! Cannot find plot name"
                    elif kind == 'CutFlow':
                        for cutflow in info:
                            #print cutflow
                            print cutflow.keys()
                            print 'CutFile = ', cutflow['CutFile']
                            for dic in cutflow['Cuts']:
                                for key, val in dic.items():
                                    print key, val
                            return
                            tmptex, tmpscreen = self.processCutFlow(cwdpath, analysis_name, name, cutflow)
                            outputTeX.extend(tmptex)
                            outputScreen += tmpscreen
                    elif kind == 'Event Files' or kind == 'Atom Options' or kind == 'Run Name' or kind == 'Description':
                        continue
                    else:
                        print "Warning: section '%s' not implemented" % kind
                # outputTeX.extend(tmptex)
                # outputScreen += tmpscreen
                self.finishRun(cwdpath, list(set(outfiles)))
        if len(outputTeX) > 0:
            self.appendLaTeX(analysis_name,outputTeX)
        return outputScreen


    def initRun(self, analysis_name, run_dir, event_file_names, option_file_name, output_files, repository_url):
        #create temp dir
        temp_path = os.path.realpath(os.path.normpath(self.work_path + '/' + analysis_name + '/' + run_dir))
        if not os.path.isdir(temp_path):
            os.makedirs(temp_path)
        if not os.path.isdir(temp_path+'/out'):
            os.makedirs(temp_path+'/out')
        #get the event file
        if len(repository_url) > 0:
            conn_type = ''
            if repository_url.upper().startswith('HTTP'):
                conn_type = 'HTTP'
            elif repository_url.upper().startswith('FILE'):
                conn_type = 'FILE'
            elif repository_url.upper().startswith('SSH'):
                conn_type = 'SSH'
            elif repository_url.upper().startswith('FTP'):
                conn_type = 'FTP'
            else:
                print "Unkown protocol for URL '" + repository_url + "'. Defaulting to ssh."
                conn_type = 'SSH'
            pos = repository_url.find("://")
            if pos >= 0:
                pos = pos+3
            else:
                pos = 0
            urlpath =  repository_url[pos:]
            for elem in event_file_names:        
                res = subprocess.Popen([self.download_script,conn_type,urlpath,elem], cwd=temp_path)
                res.wait()
        # build the option file
        optpath = option_file_name 
        if option_file_name.find('/') < 0:
            optpath = self.options_path + option_file_name
        shutil.copy(optpath, temp_path + '/options.dat')
        with open(temp_path+ '/options.dat', 'a') as fopt:
            fopt.write('add Analysis ' + analysis_name + '\n')
            if len(repository_url) > 0:
                for item in event_file_names:
                    fopt.write('add InputFile ./' + item + '\n')            
            else:
                for item in event_file_names:
                    fopt.write('add InputFile ' + item + '\n')            
            for item in output_files:
                fopt.write('add OutputFile ./' + item + '\n')
            fopt.write('launch\n')
        return temp_path

    def finishRun(self,working_path,output_files):
        # move the plots out of the way and clean up
        listOfFiles = os.listdir(working_path + '/out')
        if len(listOfFiles) > 0:
            if not os.path.isdir(self.output_path + '/figures'):
                os.makedirs(self.output_path + '/figures')
            for f in listOfFiles:
                outtmp = self.output_path + '/figures/' + os.path.basename(f)
                if os.path.isfile(outtmp):
                    os.remove(outtmp)
                shutil.move(working_path + '/out/' + f, self.output_path + '/figures')
        if self.save_output:
            if not os.path.isdir(self.output_path + '/saved_files'):
                os.makedirs(self.output_path + '/saved_files')
            for f in output_files:
                outtmp = self.output_path + '/saved_files/' + os.path.basename(f)
                if os.path.isfile(outtmp):
                    os.remove(outtmp)
                shutil.move(working_path + '/' + f, self.output_path + '/saved_files')
        # shutil.rmtree(working_path)

    def tex_num_format(self, val):
        return "$ "+("{:.2e}".format(float(val))).replace('e+0','\\cdot 10^').replace(".00"," ").replace("\\cdot 10^0","").strip()+" $"

    def initSection(self,run_name,run_description):
        evfinder = re.compile(r'Number\sof\sAtom\sMC\sevents:\s\$\s*(.*)\s*\$\.')
        res = evfinder.findall(run_description)
        tmpdesc = run_description
        if res[-1] == '0' and self.atom_numevents != 0:
            tmpdesc = evfinder.sub('Number of Atom MC events: '+self.tex_num_format(self.atom_numevents)+'.',run_description)
        texoutput = []
        screenoutput = '\n'
        screenoutput += '-' * 158 
        screenoutput += '\n'
        screenoutput += '-' * 158 
        screenoutput += '\n'
        screenoutput += 'Validation of ' + run_name + '\n'
        screenoutput += '-' * 158 
        screenoutput += '\n'
        texoutput.append('\n')
        texoutput.append('\\subsection{' + run_name + '} \n')
        texoutput.append('\n')
        tmpdesc = tmpdesc.replace('\item','\n\item').replace('\end','\n\end')
        texoutput.append(tmpdesc)
        texoutput.append('\n')
        return texoutput, screenoutput

    def appendLaTeX(self, analysis_name, output):
        with open(self.output_path + '/' + analysis_name + '.tex', 'w') as fout:
            for item in output:
                fout.write(item + '\n')


    def runAtom(self, working_path):
        res = subprocess.Popen([self.running_script],cwd=working_path,stdout=subprocess.PIPE)
        lines_iterator = iter(res.stdout.readline, b"")
        for line in lines_iterator:
            print line, # yield line
        out, err = res.communicate()
        evtsearch = re.findall(r'Processed\s(\d+)\sevents',out,re.MULTILINE)
        try:
            self.atom_numevents = int(evtsearch[-1])
        except:
            pass
        return res.returncode
        # bypasss for debugging
        # self.atom_numevents = 10000
        # return 0

    def processCmpPlot(self, working_path, analysis_name, run_name, run_idx, plot_file_name, data, index):
        cmphistos = {}
        runhistos = {}
        atom.atom_plots.getAtomCmpData(cmphistos)
        texoutput = []
        screenoutput = ''
        refhisto = None
        runhisto = None
        pname = '/' + analysis_name + '/' + data['PlotName']
        refname = ''
        pfile = working_path + '/' + plot_file_name
        if 'Reference' in data.keys():
            refname = data['Reference']
        else:
            refname = '/CMP' + pname + '/0'
        if refname[5:] in cmphistos.keys(): #need to strip the /CMP/ part
            refhisto = cmphistos[refname[5:]]
        else:
            print 'Error! Cannot find reference plot for comparison.'
            return texoutput, screenoutput
        runhistos = atom.atom_plots.getHistos([pfile])[pfile]
        if pname in runhistos.keys():
            runhisto = runhistos[pname]
        else:
            print 'Error! Cannot find plot in plot file.'
            return texoutput, screenoutput
        plot_caption = "Comparison plot for %s. Blue: ATOM, red: original." % '********'
        if 'PlotCaption' in data.keys():
            plot_caption = data['PlotCaption']
        table_caption = "List of most discrepant bins for %s." % '********'
        if 'TableCaption' in data.keys():
            table_caption = data['TableCaption']
        plot_comments = ''
        if 'Comments' in data.keys():
            plot_comments = data['Comments']
        #produce comparison table
        badbins = atom.atom_plots.findDeviations(runhisto,refhisto)
        #produce comparison plots
        pparser = atom.plotinfo.mkStdPlotParser([],'')
        pheaders = pparser.getHeaders(pname)
        poptions = pparser.getHistogramOptions(pname)
        pspecial = pparser.getSpecial(pname)
        outbasename = os.path.splitext(plot_file_name)[0]+'_P'+str(index+1)
        plotfilenames = self.producePlots(working_path,analysis_name,runhisto,refhisto,pheaders,poptions,refname,pspecial,outbasename)
        #format results
        texoutput = self.tformatter.make_plot_section(analysis_name,analysis_name + '_R' + str(run_idx+1) + '_P' +str(index+1),plot_caption,table_caption, plotfilenames, badbins, plot_comments)
        screenoutput = self.oformatter.show_plot(analysis_name,run_name,badbins)
        return texoutput, screenoutput

    def producePlots(self, working_path, ana_name, mchisto,refhisto,headers,options,refname,special,file_basename):
        exp_name = ana_name.split('_')[0]
        plotlin = atom.atom_plots.Plot()
        plotlog = atom.atom_plots.Plot()
        plotlin['Legend'] = '1'
        plotlin['LogY'] = '0'
        plotlog['Legend'] = '1'
        plotlog['LogY'] = '1'
        for key, val in headers.iteritems():
            plotlin[key] = val
            plotlog[key] = val
        refhisto.setAnnotation('ErrorBars', '1')
        refhisto.setAnnotation('PolyMarker', '*')
        refhisto.setAnnotation('ConnectBins', '0')
        refhisto.setAnnotation('Title', exp_name)
        mchisto.setAnnotation('ErrorBars', '1')
        for key, val in options.iteritems():
            mchisto.setAnnotation(key, val)
        mchisto.setAnnotation('Title', 'Atom')
        atom.atom_plots.setStyle(mchisto, 0)
        atom.atom_plots.setStyle(mchisto, 1)
        plotlin['RatioPlot'] = '1'
        plotlin['RatioPlotReference'] = refname
        plotlog['RatioPlot'] = '1'
        plotlog['RatioPlotReference'] = refname
        outputlin = ''
        outputlog = ''
        outputlin += str(plotlin)
        outputlog += str(plotlog)

        if special:
            tmpspec = ''
            tmpspec += "\n"
            tmpspec += "# BEGIN SPECIAL %s\n" % h
            tmpspec += special
            tmpspec += "# END SPECIAL\n\n"
            outputlin += tmpspec
            outputlog += tmpspec

        tmpplot = atom.atom_plots.getFlatString([refhisto,mchisto])
        outputlin += tmpplot
        outputlog += tmpplot
        flin = file_basename+'_lin'
        flog = file_basename+'_log'
        atom.atom_plots.writeOutput(outputlin, flin, working_path)
        atom.atom_plots.writeOutput(outputlog, flog, working_path)
        atom.atom_plots.runMakePlots(flin+'.dat',working_path,'tex')
        atom.atom_plots.runMakePlots(flog+'.dat',working_path,'tex')
        self.runLaTeX(working_path,flin)
        self.runLaTeX(working_path,flog)
        flin = flin + '.pdf'
        flog = flog + '.pdf'
        shutil.move(working_path+'/'+flin,working_path+'/out')
        shutil.move(working_path+'/'+flog,working_path+'/out')
        return [ flin, flog ]

    def runLaTeX(self, working_path, file_basename):
        print 'Converting plots into PDF...'
        tmp = subprocess.Popen(['latex',file_basename+'.tex'], stdout=subprocess.PIPE, cwd=working_path)
        tmp.wait()
        tmp = subprocess.Popen(['dvipdf',file_basename+'.dvi'], cwd=working_path)
        tmp.wait()
        print 'Done.'

    def processCutFlow(self, analysis_name, run_name, data):
        texoutput = ''
        screenoutput = ''
        # get Atom efficiencies
        parents = {}
        atomeffs = {}
        atomerrs = {}
        atomders = {}
        atomdererrs = {}
        cfile = data['CutFile']
        print 'reading Atom output'
        parents, atomeffs, atomerrs, atomders, atomdererrs = self.getCutEfficiencies(cfile, analysis_name)
        print 'done'            
        # get Exp efficiencies
        list_cuts = []
        name_list = []        
        eff_exp = []
        err_exp = []
        Reff_exp = []
        Rerr_exp = []
        eff_atom = []
        err_atom = []
        Reff_atom = []
        Rerr_atom = []
        derv_atom = []
        dere_atom = []
        ratio_eff = []
        ratio_eff_sig = []
        ratio_R = []
        ratio_R_sig = []
        parent_index = []
        idx_dict = { '' : 0 }
        if 'Cuts' in data.keys():
            previous_name = ''
            for idx, item in enumerate(data["Cuts"]):
                if 'Title' in item.keys():
                    list_cuts.append(item['Title'])
                    name_list.append(item['Atom Name'])                    
                    subproclist = [ 0 ]
                    if 'ProcIds' in item.keys():
                        subproclist = item['ProcIds']
                    if 'Parent Cut' in item.keys():
                        parcut = item['Parent Cut']
                    else:
                        parcut = previous_name
                    name = item['Atom Name']
                    idx_dict.update({ name : idx+1 })
                    eff_exp.append(item['Cumulative'])
                    err_exp.append(item['Cum. Error'])
                    Reff_exp.append(item['Individual'])
                    Rerr_exp.append(item['Indiv. Error'])
                    val = 100.0
                    err = [0.0, 0.0]
                    if name in atomeffs.keys():
                        val = atomeffs[name]
                        err = atomerrs[name]
                    elif len(eff_atom) > 0:
                        val = eff_atom[-1]
                        err = err_atom[-1]
                    eff_atom.append(val)
                    err_atom.append(err)
                    if item['Cumulative'] > 0:
                        ratio_eff.append(val/item['Cumulative'])
                    else:
                        ratio_eff.append(0.)
                    valpull = val - item['Cumulative']
                    tmperr = []
                    if isinstance(item['Cum. Error'], (float,int,long)):
                        tmperr = [-item['Cum. Error'], item['Cum. Error']]
                    else:
                        tmperr = item['Cum. Error']
                    errpull = math.sqrt(err[0]**2+tmperr[1]**2) if valpull > 0.0 else math.sqrt(err[1]**2+tmperr[0]**2)
                    if errpull != 0.0:
                        ratio_eff_sig.append(valpull/errpull)
                    else:
                        ratio_eff_sig.append(0.0)
                    if name in atomders.keys():
                        dval, derr = self.cutDerivative(name,parcut,atomders,atomdererrs,parents)
                        derv_atom.append(dval) 
                        dere_atom.append(derr)
                    else:
                        derv_atom.append(0.) 
                        dere_atom.append([0.,0.])                        
                    # build atom ratios
                    if name in atomeffs.keys():
                        if 'Parent ID' in item.keys():
                            parentID = item['Parent ID']
                            eff_now = atomeffs[name]
                            err_now = atomerrs[name]
                            eff_prev = eff_atom[parentID]
                            err_prev = err_atom[parentID]
                            rval = -1
                            if eff_prev > 0: rval = eff_now / eff_prev * 100.                           
                            rerr1 = 100. * math.sqrt(err_now[0]**2 + err_prev[0]**2) 
                            rerr2 = 100. * math.sqrt(err_now[1]**2 + err_prev[1]**2) 
                            rerr = [rerr1, rerr2] 
                            parent_index.append(parentID)                                                        
                            # check
                            if item['Parent Atom Name'] != name_list[parentID]:
                                print 'ERROR!!  Derected a wrong parent cut'
                                print 'Abort'
                                exit() 
                        else:
                            rval, rerr = self.cutRatio(name,parcut,atomeffs,atomerrs,parents)
                            parent_index.append(idx_dict[parcut])                            
                        Reff_atom.append(rval)
                        Rerr_atom.append(rerr)
                        if item['Individual'] > 0:
                            ratio_R.append(rval/item['Individual'])
                        else:
                            ratio_R.append(0)
                        Rvalpull = rval -item['Individual']
                        if isinstance(item['Indiv. Error'], (float,int,long)):
                            tmperr = [-item['Indiv. Error'], item['Indiv. Error']]
                        else:
                            tmperr = item['Cum. Error']
                        Rerrpull = math.sqrt(rerr[0]**2+tmperr[1]**2) if Rvalpull > 0.0 else math.sqrt(rerr[1]**2+tmperr[0]**2)
                        if Rerrpull != 0.0:
                            ratio_R_sig.append(Rvalpull/Rerrpull)
                        else:
                            ratio_R_sig.append(0.0)
                    else:
                        Reff_atom.append(100.0)
                        Rerr_atom.append([0.0,0.0])
                        ratio_R.append(1.0)
                        ratio_R_sig.append(0.0)
                        parent_index.append('')
                    if name in atomeffs.keys():
                        previous_name = name
        else:
            print "Cuts not found in '" + cfile + "'. Exiting."
            return texoutput, screenoutput
        table_caption = ''
        if 'Caption' in data.keys():
            table_caption = data['Caption']
        table_comment = ''
        if 'Comments' in data.keys():
            table_comment = data['Comments']

        #format results
        texoutput = self.tformatter.make_cutflow_table(analysis_name,run_name,table_caption, list_cuts,
                    eff_exp, err_exp, eff_atom, err_atom, ratio_eff, ratio_eff_sig, parent_index, 
                     Reff_exp, Rerr_exp, Reff_atom, Rerr_atom, ratio_R, ratio_R_sig,derv_atom,dere_atom,table_comment)
        screenoutput = self.oformatter.show_cutflow(analysis_name, run_name, name_list,
                    eff_exp, err_exp, eff_atom, err_atom, ratio_eff, ratio_eff_sig, parent_index,
                     Reff_exp, Rerr_exp, Reff_atom, Rerr_atom, ratio_R, ratio_R_sig)

        return texoutput, screenoutput


    def getCutEfficiencies(self, filename, analysis, subproclist = [0]):
        stream = open(filename, 'r') 
        data = yaml.load(stream)
        if data is not None:
            if "Analyses" in data.keys():
                analyses = data["Analyses"]
                if len(analyses) == 0:
                    print 'ERROR!!'
                    print 'There is no analyses in the output yaml file!'
                    print 'Abort'
                if analysis in analyses.keys():
                    anadata = analyses[analysis]
                    effi = {}
                    err = {}
                    dereffi = {}
                    dererr = {}
                    parent = {}                
                    if "Cuts" in anadata.keys():
                        cuts = anadata["Cuts"]
                        subprocfracs, subprocerrs = self.getSubProcessFractions(filename,subproclist)
                        for cutItem in cuts:
                            name = cutItem["Cut Name"]  
                            __parent = cutItem["Parent Cut"]
                            if len(cutItem["Data"]) > 0:
                                parent[name] = __parent
                                effi[name] = 0.
                                err[name] = [ 0., 0. ]
                                dereffi[name] = 0.
                                dererr[name] = [0., 0.]
                                for subproc in cutItem["Data"]:
                                    subid = subproc["Sub-process ID"]
                                    if not subid in subproclist: continue
                                    subval = subprocfracs[subproclist.index(subid)]
                                    effi[name] += 100.*float(subproc["Cut Value"])*subval
                                    err[name] = [ err[name][0] + float(subproc["Cut Stat Error"][0])**2 * subval**2, err[name][1] + float(subproc["Cut Stat Error"][1])**2 * subval**2]  # add error on subprocs?
                                    dereffi[name] += float(subproc["Cut Log. Derivative"])*float(subproc["Cut Value"])*subval
                                    dererr[name] = [ dererr[name][0] + float(subproc["Cut Log. Der. Stat Error"][0])**2 * (float(subproc["Cut Value"])*subval)**2, \
                                        dererr[name][1] + float(subproc["Cut Log. Der. Stat Error"][1])**2 * (float(subproc["Cut Value"])*subval)**2 ]
                                err[name] = [ -100.*math.sqrt(err[name][0]), 100.*math.sqrt(err[name][1]) ]
                                if effi[name] == 0.0:
                                    if dereffi[name] != 0.0:
                                        dereffi[name] = 10000.
                                        dererr[name] = [ -10000., 10000. ]
                                else:
                                    dereffi[name] = dereffi[name]/effi[name]*100.
                                    dererr[name] = [ -math.sqrt(dererr[name][0])/effi[name]*100., math.sqrt(dererr[name][1])/effi[name]*100. ]
                    return parent, effi, err, dereffi, dererr
            else:
                print "Error: could not find analyses in '%s'" % filename
                return None, None, None, None, None
        else:
            print "Error: could not find data in '%s'" % filename
            return None, None, None, None, None

    def cutDerivative(self, name, parent_name, der_effs, der_errs, parents):
        der_tmp = 0.0
        der_tmp_err_1 = 0.0
        der_tmp_err_2 = 0.0
        curname = name
        missing = False
        while parents[curname][0] != parent_name:
            der_tmp += der_effs[curname]**2
            der_tmp_err_1 +=  (der_effs[curname]*der_errs[curname][0])**2
            der_tmp_err_2 +=  (der_effs[curname]*der_errs[curname][1])**2
            curname = parents[curname][0]
            if not (curname in parents.keys()):
                missing = True
                break 
        if not missing:
            der_tmp = math.sqrt(der_tmp)
            if der_tmp != 0.0:
                der_tmp_err_1 = math.sqrt(der_tmp_err_1)/der_tmp
            else:
                der_tmp_err_1 = 0.0
            if der_tmp != 0.0:
                der_tmp_err_2 = math.sqrt(der_tmp_err_2)/der_tmp
            else:
                der_tmp_err_2 = 0.0
            return der_tmp, [der_tmp_err_1, der_tmp_err_2]
        else:
            return 0.0, [0.0, 0.0]

    def cutRatio(self, name, parent_name, efficiencies, errors, parents):
        curname = name
        missing = False
        while parents[curname][0] != parent_name:
            curname = parents[curname][0]
            if not (curname in parents.keys()):
                missing = True
                break 
        if not missing:
            effparent = 100.0
            errparent = [0.0, 0.0]
            if parent_name != '':
                if(parents[curname][1]):
                    effparent = efficiencies[parent_name]
                else:
                    effparent = 1. - efficiencies[parent_name]
                errparent = errors[parent_name]
            effdaughter = efficiencies[name]
            errdaughter = errors[name]
            ratio = 0.0
            errratio = [ 0.0, 0.0 ]
            if effparent != 0.0:
                ratio = effdaughter/effparent*100.
                if ratio != 0.0:
                    errratio = [ -math.sqrt((errdaughter[0]/effdaughter)**2 + (errparent[1]/effparent)**2)*ratio,
                                math.sqrt((errdaughter[1]/effdaughter)**2 + (errparent[0]/effparent)**2)*ratio]
                else:
                    errratio = errdaughter # Is this right?
            else:
                ratio = 0.0
                errratio = errparent # Is this right?
            return ratio, errratio
        else:
            return -1.0, [0.0, 0.0]

    def getSubProcessFractions(self, filename, subproclist):
        stream = open(filename, 'r') 
        data = yaml.load(stream)
        if "Cross Sections" in data.keys():
            xsecs = data["Cross Sections"]
            den = 0.
            fracs = []
            errs = []
            for xsec in xsecs:
                if xsec["Process ID"] in subproclist:
                    val = float(xsec['Cross Section'])
                    err = xsec['Cross Section Error']
                    err = map(float, err) # <==!!
                    den += val
                    fracs.append(val)
                    if val != 0.0:
                        errs.append( err )
                    else:
                        errs.append([0.0,0.0])
            if den != 0.0:
                errsump = sum([ ( x[1] /den)**2 for x in errs ])
                errsumm = sum([ ( x[0] /den)**2 for x in errs ])
                fracs = [ x/den for x in fracs ]
                errs =  [ [ -math.sqrt( (1. - 2.*x/den)* y[0]**2/den**2 + x**2/den**2 * errsumm), math.sqrt( (1. - 2.*x/den)* y[1]**2/den**2 + x**2/den**2 * errsump) ] for x, y in zip(fracs, errs) ]
                return fracs, errs
            else:
                return [ 0.0] * len(subproclist), [ [ 0.0, 0.0] ] * len(subproclist)
        else:
            return [ 0.0] * len(subproclist), [ [ 0.0, 0.0] ] * len(subproclist)

    def isCutAncestor(self, name, parent_name, parents):
        curname = name
        missing = False
        while parents[curname][0] != parent_name:
            curname = parents[curname][0]
            if not (curname in parents.keys()):
                missing = True
                break 
        return not missing




