class output_format:

    col = {
        'clear': '\033[0m',
        'black': '\033[30m',
        'red': '\033[31m',
        'green': '\033[32m',
        'yellow': '\033[33m',
        'blue': '\033[34m',
        'purple': '\033[35m',
        'cyan': '\033[36m',
        'white': '\033[37m'
    }

    def strout(self, val, n, sp, defcol='clear', warncol='clear', flag=False):

        if flag == False:
            return self.col[defcol] + str(round(val, n)).rjust(sp) + self.col['clear']
        else:
            return self.col[warncol] + str(round(val, n)).rjust(sp) + self.col['clear']

    def show_plot(self, ananame, vname, binlist):
        output = []
        sp = {}
        sp['xval'] = 6
        sp['xrange'] = 2*6+8
        sp['ratio'] = 14
        sp['ratioerr'] = 2*6+6
        sp['pull'] = 6

        output.append( ' ' * 158 )
        output.append( str(ananame) + ': ' + str(vname) )
        output.append( 'Most discrepant bins:' )
        output.append( '-' * 158 )

        expname = ananame.split('_')[0]

        head  = ('X').rjust(sp['xval'])   
        head += '[ Xmin, Xmax ]'.rjust(sp['xrange'])   
        head += ' | ' 
        head += ('R=Atom/' + expname).rjust(sp['ratio']) 
        head += ('[ ErrR+, ErrR- ]').rjust(sp['ratioerr'])  
        head += ' | ' 
        head += 'Pull'.rjust(sp['pull'])   

        output.append( head )

        output.append( '-' * 158 )

        for elem in binlist:

            warning = []
            pm = '    '
            xval = self.strout(elem[0],2,sp['xval'])
            xmin = self.strout(elem[1],2,sp['xval'])
            xmax = self.strout(elem[2],2,sp['xval'])
            pullval = (elem[3]-1.)/elem[5] if elem[3] > 1.0 else (1.-elem[3])/elem[4]
            colorbin = 'red' if abs(pullval) > 4. else ('blue' if abs(pullval) > 3. else 'clear' )
            ratio = self.strout(elem[3],2,sp['ratio'],colorbin)
            sigratio_low = self.strout(elem[4],2,6)
            sigratio_hi = self.strout(elem[5],2,6)
            pull = self.strout(pullval,2,sp['pull'])
            xerr = '[' + xmin + ',' + xmax + ' ]'
            rerr = '[' + sigratio_low + ',' + sigratio_hi + ' ]'
            line = xval
            line += xerr.rjust(sp['xrange']+2*2*4)
            line += ' | '
            line += ratio 
            line += rerr.rjust(sp['ratioerr']+2*2*4)
            line += ' | '
            line += pull

            output.append( line )

        return ('\n'.join(output) + '\n')


    def show_cutflow(self, ananame, vname, name_list, 
                     eff_exp, err_exp, eff_atom, err_atom, ratio_eff, ratio_eff_sig, 
                     i_denom, Reff_exp, Rerr_exp, Reff_atom, Rerr_atom, ratio_R, ratio_R_sig):

        expname = ananame.split('_')[0]
        sp = {}
        sp['num'] = 3
        sp['eff_exp'] = 6
        sp['err_exp'] = 5
        sp['eff_atom'] = 7
        sp['err_atom'] = 5
        sp['ratio_eff'] = 11
        sp['ratio_eff_sig'] = 9
        sp['i_denom'] = 2
        sp['Reff_exp'] = 7
        sp['Rerr_exp'] = 5
        sp['Reff_atom'] = 7
        sp['Rerr_atom'] = 5
        sp['ratio_R'] = 11
        sp['ratio_R_sig'] = 8
        sp['name'] = max([len(name) for name in name_list]) + 2

        output = []
        output.append( ' ' * 158 )

        output.append( str(ananame) + ': ' + str(vname) )

        output.append( '-' * 158 )

        head  = '#'.rjust(sp['num']) + ('cut name').rjust(sp['name']) 
        head += ' | ' 
        head += ('eff_' + expname + ' ').rjust(sp['eff_exp'] + 4 + sp['err_exp']) 
        head += ('eff_Atom ').rjust(sp['eff_atom'] + 4 + sp['err_atom']) 
        head += ('Atom/' + expname).rjust(sp['ratio_eff']) 
        head += ('Sig').rjust(sp['ratio_eff_sig'])  
        head += ' | '
        head += '#/?'
        head += ('R_' + expname + ' ').rjust(sp['Reff_exp'] + 4 + sp['Rerr_exp']) 
        head += ('R_Atom ').rjust(sp['Reff_atom'] + 4 + sp['Rerr_atom']) 
        head += ('Atom/' + expname).rjust(sp['ratio_R'] - 1) 
        head += ('Sig').rjust(sp['ratio_R_sig'])  

        output.append( head )

        output.append( '-' * 158 )

        for i in range(len(name_list)):

            warning = []
            if i_denom[i] != i: output.append( '---'  )


            pm = '  '
            xxeff_exp = self.strout(eff_exp[i], 2, sp['eff_exp'])               
            xxeff_atom = self.strout(eff_atom[i], 2, sp['eff_atom'])                        
            xxerr_exp = ' '.rjust(sp['err_exp'])
            xxerr_atom = ' '.rjust(sp['err_atom'])
            xxratio_eff = ' '.rjust(sp['ratio_eff'])
            xxratio_eff_sig = ' '.rjust(sp['ratio_eff_sig'])

            xxerr_exp = self.strout(err_exp[i], 2, sp['err_exp'])
            xxerr_atom = self.strout((err_atom[i][1]-err_atom[i][0])/2., 2, sp['err_atom'])         
            if abs(1. - ratio_eff[i]) > 0.3 and abs(ratio_eff_sig[i]) > 3.: 
                xxratio_eff = self.strout(ratio_eff[i], 2, sp['ratio_eff'], 'red')
            elif abs(1. - ratio_eff[i]) > 0.3:
                xxratio_eff = self.strout(ratio_eff[i], 2, sp['ratio_eff'], 'blue')
            else:
                xxratio_eff = self.strout(ratio_eff[i], 2, sp['ratio_eff'], 'green')
            xxratio_eff_sig = self.strout(ratio_eff_sig[i], 2, sp['ratio_eff_sig'])

            pm = ' +- '

            line  = str(i+1).rjust(sp['num']) + name_list[i].rjust(sp['name']) + ' | ' 
            line += xxeff_exp + pm + xxerr_exp 
            line += xxeff_atom + pm + xxerr_atom 
            line += xxratio_eff     
            line += xxratio_eff_sig

            #######################################################
            line += ' | '               

            xxReff_exp = ' '.rjust(sp['Reff_exp'])
            xxRerr_exp = ' '.rjust(sp['Rerr_exp'])
            xxReff_atom = ' '.rjust(sp['Reff_atom'])
            xxRerr_atom = ' '.rjust(sp['Rerr_atom'])
            xxratio_R = ' '.rjust(sp['ratio_R'])
            xxratio_R_sig = ' '.rjust(sp['ratio_R_sig'])
            if i > 0:
                xxReff_exp = self.strout(Reff_exp[i], 2, sp['Reff_exp'])
                xxRerr_exp = self.strout(Rerr_exp[i], 2, sp['Rerr_exp'])
                xxReff_atom = self.strout(Reff_atom[i], 2, sp['Reff_atom'])
                xxRerr_atom = self.strout((Rerr_atom[i][1]-Rerr_atom[i][0])/2., 2, sp['Rerr_atom'])
                if abs(1. - ratio_R[i]) > 0.3 and abs(ratio_R_sig[i]) > 5.:
                    xxratio_R = self.strout(ratio_R[i], 2, sp['ratio_R'], 'red')
                elif abs(1. - ratio_R[i]) > 0.3:
                    xxratio_R = self.strout(ratio_R[i], 2, sp['ratio_R'], 'blue')
                else:
                    xxratio_R = self.strout(ratio_R[i], 2, sp['ratio_R'], 'green')
                xxratio_R_sig = self.strout(ratio_R_sig[i], 2, sp['ratio_R_sig'])
            line += str(i_denom[i]).rjust(sp['i_denom'])
            line += xxReff_exp + pm + xxRerr_exp 
            line += xxReff_atom + pm + xxRerr_atom 
            line += xxratio_R 
            line += xxratio_R_sig

            output.append( line  )

        output.append( ' ' * 158 )

        return ('\n'.join(output) + '\n')
