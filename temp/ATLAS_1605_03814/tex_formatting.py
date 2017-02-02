
class tex_format:

    header = '''
\\documentclass[12pt]{article}

\\setlength\\topmargin{0mm}
\\setlength\\headheight{0mm}
\\setlength\\headsep{0mm} 
\\setlength\\oddsidemargin{0mm}
\\setlength\\evensidemargin{0mm}
\\setlength\\textwidth{165mm}
\\setlength\\textheight{220mm}    

\\usepackage{fancyvrb}
\\usepackage{amssymb}
\\usepackage{amsmath}
\\usepackage{hyperref}        
\\usepackage{enumerate}
\\usepackage{slashed}
\\usepackage{graphicx}
\\usepackage{color}
\\usepackage{subcaption}
\\usepackage{float}
\\usepackage{ulem}
\\usepackage{url}
\\usepackage{colortbl}

\\setcounter{tocdepth}{3}

\\newcommand{\\met}{\\slashed E_T}
\\newcommand{\\mht}{\\slashed H_T}

'''

    begin_document_CF = '''
\\documentclass[12pt]{article}

\\setlength\\topmargin{0mm}
\\setlength\\headheight{0mm}
\\setlength\\headsep{0mm} 
\\setlength\\oddsidemargin{0mm}
\\setlength\\evensidemargin{0mm}
\\setlength\\textwidth{165mm}
\\setlength\\textheight{220mm}    

\\usepackage{fancyvrb}
\\usepackage{amssymb}
\\usepackage{amsmath}
\\usepackage{hyperref}        
\\usepackage{enumerate}
\\usepackage{slashed}
\\usepackage{graphicx}
\\usepackage{color}
\\usepackage{subcaption}
\\usepackage{float}
\\usepackage{ulem}
\\usepackage{url}
\\usepackage{colortbl}

\\setcounter{tocdepth}{3}

\\newcommand{\\met}{\\slashed E_T}
\\newcommand{\\mht}{\\slashed H_T}

\\title{Validation Cut-Flow Tables}
\\date{}

\\begin{document}

\\maketitle
\\tableofcontents
\\newpage
    '''


    begin_document = '''
\\documentclass[12pt]{article}

\\setlength\\topmargin{0mm}
\\setlength\\headheight{0mm}
\\setlength\\headsep{0mm} 
\\setlength\\oddsidemargin{0mm}
\\setlength\\evensidemargin{0mm}
\\setlength\\textwidth{165mm}
\\setlength\\textheight{220mm}    

\\usepackage{fancyvrb}
\\usepackage{amssymb}
\\usepackage{amsmath}
\\usepackage{enumerate}
\\usepackage{slashed}
\\usepackage{graphicx}
\\usepackage{color}
\\usepackage{subcaption}
\\usepackage{float}
\\usepackage{ulem}
\\usepackage{url}
\\usepackage{colortbl}

\\newcommand{\\met}{\\slashed E_T}
\\newcommand{\\mht}{\\slashed H_T}

\\begin{document}
    '''
    end_document = '''        
\\end{document}
    '''

    def rout(self, val, n):
        try:
            return '$ '+ str(round(val, n)) +' $'
        except:
            return val

    def rout_2err(self, val, errm, errp, n):
        try:
            return '$ '+ str(round(val, n)) +'^{+'+ str(round(errp, n)) +'}_{'+ str(round(errm, n)) +'} $'
        except:
            return val + ' +' + errp + ' ' + errm

    def color_warnings(self, val1, threshold1, val2, threshold2, with_warning=False):
        result =[]
        if abs(val1) > threshold1:
            if abs(val2) > threshold2:
                result.append('\\cellcolor{red}\\bf ')
                if with_warning:
                    result.append('\\cellcolor{magenta} ')
            else:
                result.append('\\cellcolor{red}\\bf ')
                if with_warning:
                    result.append('\\cellcolor{cyan} ')
        result = result + [ '', '' ]
        return result[:2]

    def make_plot_section(self, ananame, vname, plot_caption, table_caption, 
                   plot_files, bin_list, plot_comment):
        
        table_lists = []
        for elem in bin_list:
            line_list = []
            line_list.append(self.rout(elem[0],2))
            line_list.append('[ ' + self.rout(elem[1],2) + ', ' + self.rout(elem[2],2) + ']')
            line_list.append(self.rout_2err(elem[3],elem[4], elem[5],2))
            pullval = (elem[3]-1.)/elem[5] if elem[3] > 1.0 else (1.-elem[3])/elem[4]
            line_list.append(self.rout(pullval,2))
            table_lists.append(line_list)

        total_lines = []
        total_lines.extend(self.plot_writer(vname,plot_caption,plot_files))
        total_lines.extend(self.plot_table_writer(table_lists, vname, table_caption))
        return total_lines


    def plot_writer(self, vname, plot_caption, plot_files, plot_comment):
        plot_lines = []
        # plot_lines.append('\\newline')
        plot_lines.append(plot_comment)
        #plot_lines.append('\\newline')
        plot_lines.append('\\begin{figure}[ht!]')       
        plot_lines.append('\\centering')
        plot_lines.append('\\begin{subfigure}[b]{.45\\textwidth}')      
        plot_lines.append('\\centering')        
        plot_lines.append('\\includegraphics[width=\\textwidth]{figures/'+plot_files[0]+'}')        
        plot_lines.append('\\end{subfigure}')
        plot_lines.append('~')
        plot_lines.append('\\begin{subfigure}[b]{.45\\textwidth}')      
        plot_lines.append('\\centering')        
        plot_lines.append('\\includegraphics[width=\\textwidth]{figures/'+plot_files[1]+'}')        
        plot_lines.append('\\end{subfigure}')
        plot_lines.append('\\caption{\\footnotesize Linear scale (left) and log scale (right) comparisons for the '+plot_caption[0].lower()+plot_caption[1:]+'}')
        plot_lines.append('\\label{fig:'+ vname + '}')
        plot_lines.append('\\end{figure}')
        return plot_lines

    def plot_table_writer(self, table_lists, vname, table_caption):

        table_lines = []
        for line_list in table_lists:
            line = ' & '.join(line_list)
            table_lines.append(line)

        #######################################

        header = ['$x_{bin}$', 
                  '$[ x_{min}, x_{max} ]$', 
                  '$R\equiv\\frac{\\rm Atom}{\\rm Exp} \pm \\sigma_R$', 
                  '$\\frac{R-1}{\\rm \\sigma_R}$']

        headline = ' & '.join(header)
            
        #################################################

        dsla = '\\' + '\\'

        texlines = []

        texlines.append('\\renewcommand{\\arraystretch}{1.3}')
        texlines.append('\\begin{table}[ht!]')
        texlines.append('\\begin{center}')
        texlines.append('\\scalebox{1.0}[1.0]{ ')

        texlines.append('\\begin{tabular}{c|c||>{\columncolor{yellow}}c|c}')
        texlines.append('\\hline')
        texlines.append( headline + dsla )
        texlines.append('\\hline')
        for line in table_lines: texlines.append( line + dsla )
        texlines.append('\\hline')

        texlines.append('\\end{tabular}')
        texlines.append('}')
        texlines.append('\\caption{\\footnotesize ' + table_caption + '}') 
        texlines.append('\\label{tab:' + vname + '}')
        texlines.append('\\end{center}')
        texlines.append('\\end{table}')

        #for t in texlines: print t
        return texlines


    def make_cutflow_table(self, ananame, vname, table_caption, 
                   texname_list, eff_exp, err_exp, eff_atom, err_atom, ratio_eff, ratio_eff_sig, 
                   i_denom, Reff_exp, Rerr_exp, Reff_atom, Rerr_atom, ratio_R, ratio_R_sig, der_atom, dererr_atom, table_comment):
        
        # err_exp[0] = ' '
        # err_atom[0] = ' '
        table_lists = []

        for i in range(len(texname_list)):

            if i_denom[i] != i: table_lists.append(['\hline'])            
            pm = ' '         
            if i > 0: pm = ' $\\pm$ '
            line_list = []
            colors4 = self.color_warnings(1. - ratio_eff[i], 0.300000001, ratio_eff_sig[i], 4.00000001, False)
            colors10 = self.color_warnings(1. - ratio_R[i], 0.300000001, ratio_R_sig[i], 4.00000001, True)
            line_list.append( str(i+1) )    # 1
            line_list.append( colors10[1] + str(texname_list[i]) )  # 2 
            line_list.append( self.rout(eff_exp[i], 2) + pm + self.rout(err_exp[i], 2) )  # 3
            line_list.append( self.rout_2err(eff_atom[i],err_atom[i][0],err_atom[i][1], 2) )  # 4
            line_list.append( colors4[0] + self.rout(ratio_eff[i], 2) ) # 5
            line_list.append( self.rout(ratio_eff_sig[i], 2) ) # 6
            line_list.append( str(i_denom[i]) ) # 7
            line_list.append( self.rout(Reff_exp[i], 2) + pm + self.rout(Rerr_exp[i], 2) )  # 8
            line_list.append( self.rout_2err(Reff_atom[i], Rerr_atom[i][0], Rerr_atom[i][1], 2)  )  # 9
            line_list.append( colors10[0] + self.rout(ratio_R[i], 2) ) # 10
            line_list.append( self.rout(ratio_R_sig[i], 2) ) # 11
            line_list.append( self.rout_2err(der_atom[i], dererr_atom[i][0], dererr_atom[i][1], 2) ) # 12

            table_lists.append(line_list)

        return self.cutflow_table_writer(table_lists, vname, table_caption, table_comment)



    def cutflow_table_writer(self, table_lists, vname, table_caption, table_comment):

        table_lines = []
        for line_list in table_lists:
            line = ' & '.join(line_list)
            table_lines.append(line)

        #######################################

        header = ['\\#', 'cut name', '$\\epsilon_{\\rm Exp}$ (\\%)', '$\\epsilon_{\\rm Atom}$ (\\%)', 
                  '$\\frac{\\rm Atom}{\\rm Exp}$', 
                  '$\\frac{({\\rm Exp} - {\\rm Atom})}{\\rm Error}$',
                  '$\\#/?$',  
                  '$R_{\\rm Exp}$ (\\%)', '$R_{\\rm Atom}$ (\\%)', 
                  '$\\frac{\\rm Atom}{\\rm Exp}$', 
                  '$\\frac{({\\rm Exp} - {\\rm Atom})}{\\rm Error}$',
                  '$\\partial \\log \\epsilon_{\\rm Atom}/\\partial \\log x_{\\rm cut}$'
                  ]

        headline = ' & '.join(header)
            
            
        #################################################

        dsla = '\\' + '\\'

        texlines = []

        # texlines.append('\\\\')
        texlines.append(table_comment)
        texlines.append('\\newline')
        texlines.append('\\renewcommand{\\arraystretch}{1.3}')
        texlines.append('\\begin{table}[h!]')
        texlines.append('\\begin{center}')
        texlines.append('\\resizebox{1.1\linewidth}{!}{ ') 
        texlines.append('\\hskip-0.8in ')

        texlines.append('\\begin{tabular}{c|l||c|c|>{\columncolor{yellow}}c|c||c|c|c|>{\columncolor{yellow}}c|c|c}')
        texlines.append('\\hline')
        texlines.append( headline + dsla )
        texlines.append('\\hline')
        for line in table_lines: 
            if '\hline' in line:
                texlines.append( line )
            else:
                texlines.append( line + dsla )
        texlines.append('\\hline')

        texlines.append('\\end{tabular}')
        texlines.append('}')
        texlines.append('\\caption{\\footnotesize ' + table_caption + '}') 
        texlines.append('\\label{tab:cflow_' + vname + '}')
        texlines.append('\\end{center}')
        texlines.append('\\end{table}')

        #for t in texlines: print t
        return texlines


