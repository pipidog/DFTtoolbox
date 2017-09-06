# this code is to read xsf file and generate the input of LMTART
import numpy as np
import numpy.linalg as linalg
import sys, os

# Main =======================================================
class dftstr:        
    def grep(self,txtlines,kws,reg=False):
        # search for line numbers from a textlines where contains the keywords
        if reg:
            ln=[ n for n,txt in enumerate(txtlines) if (re.search(kws,txt)!=None)]
        else:
            ln=[ n for n,txt in enumerate(txtlines) if (txt.find(kws)!=-1)]
            
        return ln
        
    def cart2xred(self,sublat_cart,prim):
        # convert coordinate from cartisian to reduced coordinate
        # V: the vector, prim the primitive vectros, both numpy array
        sublat_xred=np.zeros((sublat_cart.shape[0],3))
        prim_inv=linalg.inv(prim.transpose())
        for n, sublat_cart_n in enumerate(sublat_cart):
            sublat_cart_n.shape=(3,1)
            sublat_xred[n,:]=prim_inv.dot(sublat_cart_n).transpose()

        return sublat_xred
        
    def xred2cart(self,sublat_red,prim):
        sublat=np.zeros(sublat_red.shape)
        for n, sublat_red_n in sublat_red:
            sublat[n,:]=np.dot(sublat_red_n,prim)
    
    def recip_vec(self,a_vec):
        # generate the reciprocal lattice row vectros 
        b_vec=np.zeros((3,3))
        I=np.eye(3,3)
        for n in range(0,3):
            bn=2*np.pi*linalg.inv(a_vec).dot(I[:,n])
            bn.shape=(1,3)
            b_vec[n,:]=bn
            
        return b_vec
        
    def kdiv(self,a_vec,kpath,kdense=20):
        b_vec=self.recip_vec(a_vec)
        kcart=np.zeros((kpath.shape[0],3))
        # convert to cartisian coordinates
        for n, kpath_n in enumerate(kpath):
            kcart[n,:]=kpath_n.dot(b_vec)
            
        # calculate # of k-point of each segment
        kdiv=[]
        for n, kcart_n in enumerate(kcart):
            if n!=kcart.shape[0]-1:
                kdiv.append(int(round(kdense*linalg.norm(kcart[n+1,:]-kcart[n,:]))))
                
        return kdiv
        
    def ptable(self):
        # this method generates the periodic table by reading periodic_table.dat
        ptable=\
        ['N/A','H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F',\
        'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar',\
        'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co',\
        'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr',\
        'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh',\
        'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe',\
        'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu',\
        'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf',\
        'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl',\
        'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th',\
        'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es',\
        'Fm', 'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs',\
        'Mt', 'Ds', 'Rg', 'Cn', 'Nh', 'Fl', 'Mc', 'Lv', 'Ts',\
        'Og']
        
        return ptable
        
    def atom_name(self,at_list):
        # this method convert atomic number to its name and vice versa
        ptable=self.ptable()
        at_name=[]
        for at_n in at_list:
            if type(at_n)==int:
                at_name+=[ptable[at_n]]
            elif type(at_n)==str:
                at_name+=[self.grep(ptable,at_n)]
            
        return at_name
        
    def getxsf(self,wkdir,prefix,coord='red'):
        if (wkdir[-1]!='/') & (wkdir[-1]!='\\'):
            wkdir+='/'
        # read prim vectors
        with open(wkdir+prefix+'.xsf') as file:
            flines=file.readlines()        
            
        primvec_ln=self.grep(flines,'PRIMVEC')[0]
        a_vec=np.zeros((3,3))
        for n, a in enumerate(flines[primvec_ln+1:primvec_ln+4]):
            a_vec[n,:]=[float(a_i) for a_i in a.split()]
                
        primcoord_ln=self.grep(flines,'PRIMCOORD')[0]
        tot_sublat=int(flines[primcoord_ln+1].split()[0])
        atom=[]
        sublat=np.zeros((tot_sublat,3))
        ptable=self.ptable()
        for n, pos in enumerate(flines[primcoord_ln+2:primcoord_ln+tot_sublat+2]):
            tmp=pos.split()
            # check atom name or atom number
            try:
                atom.append(int(tmp[0]))
            except ValueError:
                atom.append(self.grep(ptable,tmp[0])[0])
                
            sublat[n,:]=[float(pos_i) for pos_i in tmp[1:]]            
            
        # convert cartisian coordinate to reduced coordinates
        if coord=='red':
            sublat=self.cart2xred(sublat,a_vec)
        elif coord=='cart':
            pass
        else:
            print('Error: struct.dftstr.getxsf, coord={0} not defined!'.format(coord))
            sys.exit()

        return atom, a_vec, sublat
        
    def getkpf(self,wkdir,prefix):
        if (wkdir[-1]!='/') | (wkdir[-1]!='\\'):
            wkdir+='/'
        with open(wkdir+prefix+'.kpf') as file:
            flines=file.readlines()
            
        kline_start=self.grep(flines,'Real form of k-point coordinates')[0]+1
        kline_end=self.grep(flines,'END of FILE')[0]-4
        tot_kpt=kline_end-kline_start+1
        klabel=[]
        kpath=np.zeros((tot_kpt,3))
        for n, kpt in enumerate(flines[kline_start:kline_end+1]):
            tmp=kpt.split()
            klabel.append(tmp[-1])
            kpath[n,:]=[float(kpt_i) for kpt_i in tmp[0:3]]

        return klabel, kpath
        
    def printxsf(self,wkdir,atom,a_vec,sublat,prefix='printxsf'):
        with open(wkdir+prefix+'.xsf','w') as file:
            file.write('CRYSTAL\n')
            file.write('    PRIMVEC\n')
            [file.write('     {0[0]:12.8f}  {0[1]:12.8f}  {0[2]:12.8f}\n'.format(a_vec_n)) for a_vec_n in a_vec]
            file.write('    PRIMCOORD\n')
            file.write('      {0}  1\n'.format(len(atom)))
            [file.write('      {0:2s} {1[0]:12.8f} {1[1]:12.8f} {1[2]:12.8f}\n'.format(str(atom[n]),sublat[n]))\
            for n in range(0,len(atom))]
        
    def str2input(self,form,atom,a_vec,sublat):
        ptable=self.ptable()
        # convert atomic name format to atomic number format
        if type(atom[0])==str:
            atom=[ self.grep(ptable,atn)[0] for atn in atom]
        
        #spec=sorted(tuple(set(atom)))         
        spec=[]
        [spec.append(at_n) for at_n in atom if spec.count(at_n)==0]
        if form=='abinit':
            print('')
            print(' =========== Abinit Structure Inputs ===========')
            print('chkprim 0   # warnin on non-primitive cell')
            print('nsym 0      # auto symmetry finder')
            print('ntypat %i' % len(set(atom)))
            print('natom %i' % len(atom))
            print('znucl '+'%i '*len(spec) % tuple(spec))
            print('typat ')
            for n, atom_n in enumerate(atom):
                [print('%3i ' % (n+1),end='') for n, spec_n in enumerate(spec) if spec_n==atom_n]
                if (divmod(n,10)[1]==9) | (n+1==len(atom)):
                    print('')
            print('acell 3*1.0 angstrom')
            print('rprim')
            for n in range(0,3):
                print('%12.8f %12.8f %12.8f' % tuple(a_vec[n,:].tolist()))
            
            print('xred')    
            for n in range(0,len(atom)):
                print('%12.8f %12.8f %12.8f' % tuple(sublat[n,:].tolist()))
        elif form=='qe':
            print('')
            print(' =========== Quantum Espresso Structure Inputs ===========')
            print('ibrv = 0,')
            print('celldm(1) = 1.89,')
            print('ntyp = %i' % len(spec))
            print('nat = %i' % len(atom))
            print('/')
            print('CELL_PARAMETERS alat')
            for n in range(0,3):
                print('    %12.8f %12.8f %12.8f' % tuple(a_vec[n,:].tolist()))
            print('ATOMIC_SPECIES')
            for spec_n in spec:
                print('    %2s     1.000   PSP_NAME.upf' % ptable[spec_n]) 
            print('ATOMIC_POSITIONS crystal')
            for n in range(0,len(atom)):
                print('%2s  ' % ptable[atom[n]],end='')
                print('%12.8f %12.8f %12.8f 0 0 0' % tuple(sublat[n,:]))
        elif form=='elk':
            print('')
            print(' =========== Elk Structure Inputs ===========')
            print('primcell')
            print('  .true.')
            print('scale')
            print('  1.89')
            print('avec')
            [print('  %12.8f'*3 % (a_n[0],a_n[1],a_n[2])) for a_n in a_vec]
            print('atoms')
            print('%3i        :nspecies' % len(spec))
            for spec_n in spec:
                sublat_label=np.nonzero(np.array(atom)==spec_n)[0]
                print('%2s.in      :spfname' % ptable[spec_n])
                print('%2i         :natoms; atposl, bfcmt below' % len(sublat_label))
                for label in sublat_label:
                    print('%12.8f %12.8f %12.8f   0.0   0.0   0.0' % tuple(sublat[label,:].tolist())) 
        elif form=='lmtart':
            print('')
            print(' =========== LMTART Structure Inputs ===========')
            print(' ------------- ini -------------' )
            print(' <SECTION=MAIN>          ! MAIN ATOMIC DATA:')
            print(' Natom =%2i               ! # of atoms in the unit cell' % len(atom))
            print(' Nsort =%2i               ! # of sorts in the unit cell' % len(spec))
            print(' Nspin =2                ! # of spins')
            print(' Norbs =1                ! 1-without/2 -with spin orbit coupling')
            print(' Par0  = 1.88973         ! lattice parameter in a.u.')
            print(' Is(:) =',end='')
            for n, atom_n in enumerate(atom):
                [print('%2i ' % (n+1),end='') for n, spec_n in enumerate(spec) if spec_n==atom_n]
                if (n+1==len(atom)):
                    print('! atom-to-sort pointer array')
            for n, spec_n in enumerate(spec):
                print('<SECTION=SORT>      ! SORT DATA')
                print(' Name = %2s           ! ATOM LABEL'% ptable[spec_n])
                print(' Znuc = %3i          ! nuclear charge' % spec_n)
                print(' Split = 0.0         ! initial splitting')
            print('<SECTION=FFTS>      ! FFT GRIDS:')
            print(' nDiv(:) =6 6 6      ! Tetrahedron mesh')
            print(' ------------- str -------------' )
            print('<SECTION=CTRS>             ! CONTROL STRUCTURE:')
            print(' Natom =%2i                 ! # of atoms' % len(atom))
            print(' Icrd = 0                  ! basis coordinates')
            print(' BtoA = 1.0                !  b over a ratio')
            print(' CtoA = 1.0                !  c over a ratio')
            print('<SECTION=TRAN>           ! PRIMITIVE TRANSLATIONS:')
            for a_vec_n in a_vec:
                print('%12.8f, %12.8f, %12.8f'% tuple(a_vec_n))
            print('<SECTION=BASS>           ! BASIS ATOMS :')
            for sublat_n in sublat:
                print('%12.8f, %12.8f, %12.8f' % tuple(sublat_n))
        else:
            print('Error: structure.str_input, output format does not support')

    def kpath2input(self,form,a_vec,kpath,klabel=None,kdense=20):
        if klabel==None:
            klabel=[]
            for n in range(0,kpath.shape[0]):
                klabel.append('k'+str(n))
                
        kdiv=self.kdiv(a_vec,kpath,kdense)
        if form=='abinit':
            print('')
            print(' =========== Abinit k-path Inputs ===========')
            print('kptopt -%i' % len(kdiv))
            print('ndivk '+'%i '*len(kdiv) % tuple(kdiv))
            print('kptbounds')
            [print('%12.8f  '*3 % tuple(kpath_n)) for kpath_n in kpath]
        elif form=='qe':
            print('')
            print(' =========== Quantum Espresso k-path Inputs ===========')
            print('K_POINTS crystal_b')
            print('%i' % (len(kdiv)+1))
            kdiv.extend([1])
            for n, kpath_n in enumerate(kpath):
                print('%12.8f  '*3 % tuple(kpath_n),end='')
                print('%2i  ! %2s' % (kdiv[n],klabel[n]))
        elif form=='elk':
            print('')
            print(' =========== Elk k-path Inputs ===========')
            print('! band path: '+'%s-'*kpath.shape[0] % tuple(klabel))
            print('plot1d')
            print('    %i  %i : nvp1d, npp1d '% (kpath.shape[0],sum(kdiv)))
            for n, kpath_n in enumerate(kpath):
                if n==0:
                    print('  %12.8f'*3 % tuple(kpath_n), end='')
                    print('  : vlv1d')
                else:
                    print('  %12.8f'*3 % tuple(kpath_n))
        elif form=='lmtart':
            print('')
            print(' =========== LMTART k-path Inputs ===========')
            print('<SECTION=DIRS>           ! HIGH-SYMMETRY LINES:')
            print('%i ' % len(kdiv))
            for n, kpath_n in enumerate(kpath):
                if n!=kpath.shape[0]-1:
                    print('%s-%s' % (klabel[n],klabel[n+1]))
                    print('%12.8f, %12.8f, %12.8f,' % tuple(kpath[n,:]),end='')
                    print('%12.8f, %12.8f, %12.8f' % tuple(kpath[n+1,:]))
        else:
            print('Error: structure.kpath_input, output format does not support')
            
    def xsf2input(self,form,wkdir,prefix):
        atom, a_vec, sublat=self.getxsf(wkdir,prefix)
        self.str2input(form,atom,a_vec,sublat)
    
    def kpf2input(self,form,wkdir,prefix,kdense=20):
        atom, a_vec, sublat=self.getxsf(wkdir,prefix)
        klabel, kpath=self.getkpf(wkdir,prefix)
        self.kpath2input(form,a_vec,kpath,klabel,kdense)
            
            
            
# test =============
if __name__=='__main__':
	# parameters ==================================
	wkdir=os.path.dirname(os.path.realpath(__file__))+'/tests'
	prefix='Li2OsO3'
	form='lmtart'
	# main ========================================
	mystr=structure()
	mystr.xsf2input(form,wkdir,prefix)
	mystr.kpf2input(form,wkdir,prefix)

    
