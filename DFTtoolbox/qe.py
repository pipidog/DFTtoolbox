import numpy as np
import numpy.linalg as la
from DFTtoolbox.postproc import dftpp
from DFTtoolbox.struct import dftstr
import sys, os, re, shutil

class init(dftstr):
    def __init__(self,wkdir):
        if (wkdir[-1] is not '/') or (wkdir[-1] is not '\\'):
            wkdir+='/'
        self.wkdir=wkdir
    
    def ground(self,prefix,soc,mag,dftu,kdense=20):
        # collecing all input infomation ----------------
        ptable=self.ptable()
        #spec=sorted(tuple(set(atom)))         
        atom, a_vec, sublat=self.getxsf(self.wkdir,prefix)
        klabel, kpath=self.getkpf(self.wkdir,prefix) 
        kdiv=self.kdiv(a_vec,kpath,kdense)
        if type(atom[0])==str:
            atom=[ self.__grep__(ptable,atn)[0] for atn in atom]          
        # convert atomic name format to atomic number format
        spec=[]
        [spec.append(at_n) for at_n in atom if spec.count(at_n)==0]
        
        # generate pw.scf.in file
        file=open(self.wkdir+'pw.scf.in','w')
        file.write('&CONTROL\n')
        file.write('  title = \'{0}\',\n'.format(prefix))
        file.write('  prefix = \'pwscf\' ,\n')
        file.write('  calculation = \'scf\',\n')
        file.write('  restart_mode = \'from_scratch\' ,\n')
        file.write('  outdir = \'./outdir\',\n')
        file.write('  pseudo_dir = \'./\' ,                 ! pseudopotential folder\n')
        file.write('  verbosity = \'high\',                 ! must use \'high\' to print bands\n')
        file.write('  wf_collect = .true. ,\n')
        file.write('  max_seconds =  21000 ,\n')
        file.write('/\n')
        file.write('&SYSTEM\n')
        file.write('  ! < basic parameters > ===================\n')
        file.write('  ibrav = 0,\n')
        file.write('  celldm(1) = 1.89,\n')
        file.write('  ntyp = {0} ,\n'.format(len(spec)))
        file.write('  nat = {0} ,\n'.format(len(atom)))
        file.write('  ecutwfc = 33 ,                      ! > 30 for uspp/paw, >40 for ncpp\n')
        file.write('  ecutrho = 330 ,                     ! 4*ecutwfc for paw/ncpp, 10*ecutwfc for uspp\n')
        file.write('  occupations = \'smearing\' ,\n')
        file.write('  degauss = 0.02 ,\n')
        file.write('  ! < spin wavefunctions > =================\n')
        if (soc is 'on'):
            file.write('  nspin = 4,\n')
            file.write('  noncolin = .true. ,\n')
            file.write('  lspinorb = .true. ,\n')
            file.write('  starting_spin_angle = .false.       ! whether use soc-wf as initial wf\n')
        elif (soc is 'off') and (mag is 'on'):
            file.write('  nspin = 2,\n')
        elif (soc is 'off') and (mag is 'off'):
            file.write('  nspin = 1,\n')
        if (mag is 'on'):
            file.write('  ! < magnetic constrains > ================\n')
            file.write('  constrained_magnetization = \'none\'  ! constrain scheme\n')
            file.write('  lambda = 1.0                        ! constraint parameters, 0~1000\n')   
            file.write('  tot_magnetization=-1                ! major-s - minor-s (LSDA) > 0, if constrained_magnetization=\'total\'\n')
            file.write('  ! < initial moments > ====================\n')
            for n in range(0,len(spec)):
                file.write('  starting_magnetization({0}) = 0.0,    ! -1~+1, % of polar-valence-e\n'.format(n+1))
                if (soc is 'on'):
                    file.write('  angle1({0}) = 0.0,                    ! angle between m and z\n'.format(n+1))
                    file.write('  angle2({0}) = 0.0,                    ! angle between m on x-y plane and x\n'.format(n+1))
        if (dftu is 'on'):
            file.write('  ! < DFT+U parameters > ===================\n')
            file.write('  lda_plus_u = .true. \n')
            if (soc is 'off'):
                file.write('  lda_plus_u_kind = 0\n')
                [file.write('  Hubbard_U({0}) = 0.0                  ! (eV)\n'.format(n+1)) for n in range(0,len(spec))]
                [file.write('  Hubbard_J0({0}) = 0.0                 ! (eV)\n'.format(n+1)) for n in range(0,len(spec))]            
            elif (soc is 'on'):
                file.write('  lda_plus_u_kind = 1\n')            
                [file.write('  Hubbard_U({0}) = 0.0                  ! (eV)\n'.format(n+1)) for n in range(0,len(spec))]
                [file.write('  Hubbard_J(1,{0}) = 0.0                ! (eV)\n'.format(n+1)) for n in range(0,len(spec))]
        file.write('/\n')
        file.write('&ELECTRONS\n')
        file.write('  electron_maxstep = 100,\n')
        file.write('  conv_thr = 1.D-5 ,\n')
        file.write('  mixing_beta = 0.7 ,\n')
        file.write('  diagonalization = \'david\' ,\n')
        file.write('/\n')
        file.write('ATOMIC_SPECIES\n')
        for spec_n in spec:
            file.write('    %2s     1.000   PSP_NAME.upf       ! UPF file\n' % ptable[spec_n]) 
        file.write('CELL_PARAMETERS alat\n')
        for n in range(0,3):
            file.write('    %12.8f %12.8f %12.8f\n' % tuple(a_vec[n,:].tolist()))
        file.write('ATOMIC_POSITIONS crystal\n')
        for n in range(0,len(atom)):
            file.write('  %2s  ' % ptable[atom[n]])
            file.write('%12.8f %12.8f %12.8f 0 0 0\n' % tuple(sublat[n,:]))
        file.write('K_POINTS automatic\n  ')
        [file.write('{0} '.format(int(np.round(45/la.norm(a_n))))) for a_n in a_vec]
        file.write(' 1  1  1')
        file.close()
        
        # read pw.scf.in file
        with open(self.wkdir+'pw.scf.in','r') as file:
            flines=file.readlines()
        
        ln=self.grep(flines,'calculation = ')
        flines[ln[0]]='  calculation = \'bands\'\n'
        
        ln=self.grep(flines,'&SYSTEM')
        flines.insert(ln[0]+1,'  ! keep all varialble the same as pw.scf.in\n')
        
        ln=self.grep(flines,'K_POINTS')
        flines=flines[0:ln[0]]
        flines.append('K_POINTS crystal_b\n')
        flines.append('  {0}\n'.format(len(kdiv)+1))
        kdiv.extend([1])
        for n, kpath_n in enumerate(kpath):
            flines.append('{0[0]:12.8f} {0[1]:12.8f} {0[2]:12.8f} {1:2d} ! {2:2s}\n'.\
            format(kpath_n,kdiv[n],klabel[n]))
            
        with open(self.wkdir+'pw.bands.in','w') as file:
            file.writelines(flines)
            
        # generate projwfc.fat.in 
        with open(self.wkdir+'projwfc.fat.in','w') as file:
            file.write('&projwfc\n')
            file.write('  outdir=\'./outdir\'\n')
            file.write('  prefix=\'pwscf\'\n')
            file.write('  lsym=.false.\n')
            file.write('  filproj = \'fatband\'\n')
            file.write('/')
        # generate projwfc.dos.in 
        with open(self.wkdir+'projwfc.dos.in','w') as file:
            file.write('&projwfc\n')
            file.write('  outdir=\'./outdir\'\n')
            file.write('  prefix=\'pwscf\'\n')
            file.write('  lsym=.true.\n')
            file.write('  filproj = \'pdos\'\n')
            file.write('/')
        print('Note:')
        print('1. pw.scf.in, pw.bands.in, projwfc.dos.in, projwfc.fat.in are generated')
        print('2. If you want DOS, you must run it immediately after scf run')
        print('3. Recommend procedure for ground state calculations:')
        print()
        print(' $ > pw.x < pw.scf.in < pw.scf.out')
        print(' $ > projwfc.x < projwfc.dos.in < projwfc.dos.out')
        print(' $ > pw.x < pw.bands.in < pw.bands.out')
        print(' $ > projwfc.x < projwfc.fat.in < projwfc.fat.out')
        print()
        print(' * If you have the above file names, you will need to tell DFTtoolbox in postprocess')
        
    def relax(self,prefix):
        # collecing all input infomation ----------------
        ptable=self.ptable()
        #spec=sorted(tuple(set(atom)))         
        atom, a_vec, sublat=self.getxsf(self.wkdir,prefix)
        if type(atom[0])==str:
            atom=[ self.__grep__(ptable,atn)[0] for atn in atom]          
        # convert atomic name format to atomic number format
        spec=[]
        [spec.append(at_n) for at_n in atom if spec.count(at_n)==0]
        
        # generate pw.relax.in file
        file=open(self.wkdir+'pw.relax.in','w')
        file.write('&CONTROL\n')
        file.write('  title = \'{0}\',\n'.format(prefix))
        file.write('  calculation = \'vc-relax\', ! set to \'relax\' if want cell fixed\n')
        file.write('  restart_mode = \'from_scratch\' ,\n')
        file.write('  outdir = \'./outdir\',\n')
        file.write('  pseudo_dir = \'./\' ,\n')
        file.write('  prefix = \'pwscf\' ,\n')
        file.write('  verbosity = \'high\',\n')
        file.write('  etot_conv_thr = 1.0D-3 ,\n')
        file.write('  forc_conv_thr = 1.0D-3 ,\n')
        file.write('  tstress = .true. ,\n')
        file.write('  tprnfor = .true. ,\n')
        file.write('  nstep = 50\n')
        file.write('  wf_collect = .true. ,\n')
        file.write('  max_seconds =  21000 ,\n')
        file.write('/\n')
        file.write('&SYSTEM\n')
        file.write('  ibrav = 0,\n')
        file.write('  celldm(1) = 1.89,\n')
        file.write('  ntyp = {0} ,\n'.format(len(spec)))
        file.write('  nat = {0} ,\n'.format(len(atom)))
        file.write('  ecutwfc = 33 , ! > 30 for uspp/paw, >40 for ncpp\n')
        file.write('  ecutrho = 330 ,! 4*ecutwfc if ncpp/paw, 10*ecutwfc if uspp \n')
        file.write('  occupations = \'smearing\' ,\n')
        file.write('  degauss = 0.02 ,\n')
        file.write('  nspin = 1,\n')
        file.write('/\n')
        file.write('&ELECTRONS\n')
        file.write('  electron_maxstep = 100,\n')
        file.write('  conv_thr = 1.D-5 ,\n')
        file.write('  mixing_beta = 0.7 ,\n')
        file.write('  diagonalization = \'david\' ,\n')        
        file.write('/\n')
        file.write('&IONS\n')
        file.write('  ion_dynamics = \'bfgs\'\n')
        file.write('/\n')
        file.write('&CELL\n')
        file.write('  cell_dynamics = \'bfgs\' ,\n')
        file.write('  press_conv_thr = 0.5 ,\n')
        file.write('  cell_dofree = \'all\' , ! check if you want constraints\n')
        file.write('/\n')
        file.write('ATOMIC_SPECIES\n')
        for spec_n in spec:
            file.write('    %2s     1.000   PSP_NAME.upf\n' % ptable[spec_n]) 
        file.write('CELL_PARAMETERS alat\n')
        for n in range(0,3):
            file.write('    %12.8f %12.8f %12.8f\n' % tuple(a_vec[n,:].tolist()))
        file.write('ATOMIC_POSITIONS crystal\n')
        for n in range(0,len(atom)):
            file.write('  %2s  ' % ptable[atom[n]])
            file.write('%12.8f %12.8f %12.8f 0 0 0\n' % tuple(sublat[n,:]))
        file.write('K_POINTS automatic\n  ')
        [file.write('{0} '.format(int(np.round(45/la.norm(a_n))))) for a_n in a_vec]
        file.write(' 1  1  1')
        file.close()
        print('pw.relax.in generated. To run it, simply use:')
        print('  $ >  pw.x < pw.relax.in > pw.relax.out')
        
class postproc(dftpp):
    '''
    **** ABSTRACT ****
    This module defines the "espresso" class. This class contains
    several tools to prepare and postprocess quantum espresso data.
    **** INITLALIZE ****    
    INHERIT: 
    DFTtoolbox.postprocess.postprocess
    
    __init__(self,wkdir)
        => To initialize this class requires to input the working directory.
        All QE calculated data, such as band.dat, projwfc.dat, etc, for 
        postprocessing should put in this folder. Also any output file generated 
        by this class is also put here. 
    
    **** METHODS ****
    _kdiv_qe(self,kdiv)
        * convert kdiv from qe format to DFTtoolbox format
        => [kdiv]: kdiv used in QE input file
    band_read(self,Ef,spin=1,bandat='bands.dat')
        * read the output file of bands.x in QE and store it in DFTtoolbox
        standard format in band.npz
        => [Ef]: float, 
        Fermi energy in unit of eV (check QE scf.out file)
        => [spin]: int, 1 / 2, default=1
        if spin separate (i.e. spin-polarized LDA / GGA w/o SOC), spin=2
        else spin=1
        => [bandat]: str, defualt='bands.dat'
        * the file name of the output file of bands.x
    band_plot(self,kdiv='default',klabel='default',Ebound='default',lw=2,fontsize=18):
        * this method will automatically load band.npz and plot the band
        structure. Therefore, run band_read is required before calling this 
        method.
        => [kdiv]: list, int
        the kdiv used in QE input file, default =[0,tot_k], i.e. no div
        => [klabel]: list, str
        the name of each high symmetry points, e.g.  ['$\Gamma','X','M']
        default=['k1','k2',..,'kn']
        => [Ebound]: list, two float values, e.g. [-3.5, 3.5]
        upper and lower energy bounds, defult=[min(Ek),max(Ek)]
        => [lw]: int, default=2
        the line width 
        => [fontsize]: int, default=18
        the fontsize of the plot
    fatband_read(self,Ef,bandat='bands.dat',projdat='projwfc.dat',projout='projwfc.out')
        * this method read the fatband data of QE output. Note that, one needs to run 
        bands.x first to obtain bands.dat. Then run projwfc.x to obtain the output
        projwfc.dat and its log file projwfc.out (where state information is stored). 
        A typical bands.in file should look like:
        
        &bands
          prefix='Li2OsO3'
          outdir='./outdir'
          filband='bands.dat'
        /
        
        Also a typical projwfc.in file should look like:
         
         &projwfc
            outdir='./outdir'
            prefix='Li2OsO3'
            lsym=.false.
            filproj = 'projwfc.dat'
         /
        
        when running projwfc.in by projwfc.x < projwfc.in > projwfc.out 
       
        the projwfc.out should also kept for this code to read. Once the method is
        performed, the code will automatically call band_read to read bands.dat and
        store the data in band.npz. As for the data read from projwfc will be stored
        in fatband.npz.
        
        =>[Ef]: float
        the Fermi energy found in pw.scf.out
        => [bandat]: string, default=bands.dat
        the file name of the output data by bands.x
        =>[projdat]: string, default=projwfc.dat
        the file name of the output data by projwfc.x
        =>[projout]: string, default=projwfc.out
        the file name of the log file when running projwfc.x    
    
    fatband_plot(self,grp_type,state_grp,kdiv='default',klabel='default',\
        Ebound='default',ini_fig_num=1,marker_size=30,colorcode='b',fontsize=18):
        * this method plots the fatband sturctures. Note that it will load
        band.npz and fatband.npz to fetch all fatband data. Therefore, one must
        perform fatband_read (will automatically call band_read) before using 
        this method. 
        
        => [grp_tyep]: string, 'atom' / 'state'
        tells the code the meaning of state_grp. if 'atom', the code will interpret
        all values in state_grp as the atom label. if 'state', it means state label.
        => [state_grp]: list, e.g:[[1,2,3],[7,8,10],[12,17]]
        a list which shows the states to be grouped into a fatband plot. 
        therefore, the above example will output three fatband plots where
        the first plot shows the wieght of sum over state 1 2 and 3. 
        Note that the the meaning of state_grp depends on grp_tyep
        => [kdiv], [klabel], [Ebound]
        the same as band_plot
        => [ini_fig_num]: int, default=1
        in case you plot many plots in your code, you don't want it conflict with 
        others. so you can set the ini_fig_num to tell the code the figure number
        of each fatband plot. e.g. if ini_fig_num=3 and you have 3 state_grp, then
        your fatband plot will be fig.3 ~ fig.5
        => [marker_size]: int, default=30
        the size of your fatband markers
        => [colorcode]: str, 'r' / 'g' / 'b' , default: 'b'
        the default colormap for fatband plot. 'r': red, 'g': green, 'b': blue
        => [fontsize]: int, default=18
        the fontsize.
        
    **** Version ****
    02/01/2016: first built
    **** Comment ****
    1. Run /test/postprocess_test.py to test this module 
    '''
    def __init__(self,wkdir):
        if (wkdir[-1]!='/') & (wkdir[-1]!='\\') :
            wkdir=wkdir+'/'
        self.wkdir=wkdir
        
    def _kdiv_qe(self,kdiv):
        # convert kdiv from QE format to postprocess format
        if kdiv=='defaut':
            kdiv_eq='default'
        else:
            kdiv_qe=[0]
            [kdiv_qe.append(val+kdiv_qe[n]) for n, val in enumerate(kdiv)]           
            kdiv_qe.pop(-1)            
        
        return kdiv_qe
        
    def band_read(self,Ef,bandfile='pw.bands.out'):
        # read band eigenvalues fromm pw.bands.out, note to set verbosity='high'
        print('band_read start ...')
        with open(self.wkdir+bandfile,'r') as file:
            bandlines=file.readlines()
            
        # check band parameters
        if self.grep(bandlines,'SPIN')!=[]:
            spin=2
        else:
            spin=1
        k_ln=self.grep(bandlines,'bands (ev)')    
        tot_k= int(len(k_ln)/spin)
        
        # check how many lines with band eigenvalues
        for n, bandlines_n in enumerate(bandlines[k_ln[0]+2:k_ln[1]]):
            if bandlines_n.split()==[]:
                tot_banvalline=n
                break
                    
        # define function to read band states @ k-th point
        def kband_read(kn):
            kpoint=bandlines[k_ln[kn]].replace('-','  -').split()[2:5]
            kpoint=[ float(kpoint_n) for kpoint_n in kpoint]
            band_val=[]
            for bandlines_n in bandlines[k_ln[kn]+2:k_ln[kn]+3+tot_banvalline]:
                band_val.extend(bandlines_n.replace('-','  -').split())
            
            return kpoint, band_val
        kpoint, band_val=kband_read(0)
        tot_ban=int(len(band_val)*spin) 
        print(' => total band (with spin) = %d, total k_point = %d' % (tot_ban, tot_k))
        
        # read eigenvalues
        Ek=np.zeros((tot_k,tot_ban))
        k_point=np.zeros((tot_k,3))
        for kn, k_ln_n in enumerate(k_ln):
            kpoint, band_val=kband_read(kn)
            if kn <=tot_k-1:
                Ek[kn,0:int(tot_ban/spin)]=band_val
            else:
                Ek[kn-tot_k,int(tot_ban/spin):]=band_val
            
        Ek=Ek-Ef
        np.savez(self.wkdir+'band.npz',bandfile=bandfile,Ek=Ek,k_point=k_point,spin=spin)
        
    def band_plot(self,kdiv='default',klabel='default',Ebound='default',lw=2,fontsize=18):
        print('band_plot start ...')
        band=np.load(self.wkdir+'band.npz')

        dftpp.band_plot(self,band['spin'],band['Ek'],\
        0,self._kdiv_qe(kdiv),klabel,Ebound,lw,fontsize,self.wkdir)
        
    def fatband_read(self,Ef,projout='projwfc.fat.out',projprefix='fatband'):
        print('fatband_read start ...')
        # extract file name
        flist=os.listdir(self.wkdir)
        
        SOC_file=[]
        up_file=[]
        dn_file=[]
        for fname in flist:
            if fname==projprefix:
                SOC_file=fname
            elif fname==projprefix+'.projwfc_up':
                up_file=fname
            elif fname==projprefix+'.projwfc_down':
                dn_file=fname
        
        if (SOC_file!=[]) & (up_file==[]) & (dn_file==[]): # only SOC
            SOC=1
            spin=1
            proj_fname=[SOC_file]
        elif (SOC_file==[]) & (up_file!=[]) & (dn_file==[]): # only up
            SOC=0
            spin=1
            proj_fname=[up_file]
        elif (SOC_file==[]) & (up_file!=[]) & (dn_file!=[]): # up & dn
            SOC=0
            spin=2
            proj_fname=[up_file,dn_file]
        else:
            print('SOC_file='+SOC_file+', up_file='+up_file+', dn_file='+dn_file)
            print('Error: projwfc files does not reflect spin correctly!')
            print('       Note: projwfc_up should not appear when SOC is included!')
            sys.exit()  
        print(' => file to read:')
        print(proj_fname)
        # construct band parameters ============================
        # extract header information
        print('reading fatband data ...')
        with open(self.wkdir+proj_fname[0],'r') as file:
            header=[file.readline() for n in range(0,2)][1].split()
            tot_at=int(header[-2])
            tot_spec=int(header[-1])
            [file.readline() for n in range(tot_at+tot_spec+5)]
            header=file.readline().split()
            tot_state=int(header[0])*spin
            tot_k=int(header[1])
            tot_ban=int(header[2])*spin
        print(' => tot_at={0}, tot_spec={1}'.format(tot_at, tot_spec))            
        print(' => spin={0}, tot_state={1}, tot_k={2}, tot_ban={3}'\
        .format(spin, tot_state, tot_k, tot_ban))
        
        # read projwfc.dat (usually very large file) ===========
        Ek_weight=np.zeros((tot_k,tot_ban,tot_state))
        for spin_n in range(1,spin+1):
            # read weight from .projwfc_
            print(' => reading Ek_weight from '+proj_fname[spin_n-1])
            
            file=open(self.wkdir+proj_fname[spin_n-1],'r')            
            [file.readline() for n in range(0,9+tot_at+tot_spec)]

            print(' for spin={0}, total state={1}'.format(spin_n,int(tot_state/spin)))
            print('   ',end='')
            for n in range(0,int(tot_state/spin)):
                file.readline()
                weight=np.array(\
                [float(file.readline().split()[2])for m in range(0,tot_k*int(tot_ban/spin))]\
                ).reshape((tot_k,int(tot_ban/spin)))
                
                if spin_n==1:
                    Ek_weight[:,0:int(tot_ban/spin),n]=weight
                else:
                    Ek_weight[:,int(tot_ban/spin):,n+int(tot_state/spin)]=weight
                
                
                print('%3d ' % (n),end='')
                if (divmod(n,10)[1]==9):
                    print()
                    print('   ',end='')
                
            print()
            file.close()
            
        # read band structure ==================================
        print('=> read band eigenvalues from '+projout)
        with open(self.wkdir+projout) as file:
            flines=file.readlines()
            
        ban_ln=self.grep(flines,'==== e')
        if len(ban_ln)!=tot_k*tot_ban:
            print('Error: number of bands ={0} in projwfc.out =/= tot_k*tot_ban*spin'.format(len(ban_ln)))           
            sys.exit()
        
        Ek=np.zeros((tot_k,tot_ban))
        for kn in range(0,tot_k*spin):
            ban_val=[float(flines[ban_ln_n].replace('e(','e(  ').split()[4])\
                    for ban_ln_n in ban_ln[0+kn*int(tot_ban/spin):\
                    int(tot_ban/spin)+kn*int(tot_ban/spin)]]
            if kn+1<=tot_k:
                Ek[kn,0:int(tot_ban/spin)]=ban_val
            else:
                Ek[kn-tot_k,int(tot_ban/spin):]=ban_val
        Ek=Ek-Ef

        # read state info from projwfc.out =====================
        print('=> reading state info from '+projout)
        with open(self.wkdir+projout,'r') as file:
            txt=[file.readline() for n in range(0,tot_state+100)]
        
        if self.grep(txt,'m_j=')!=[]:
            SOC=1
        else:
            SOC=0
            
        count=0
        state_info=[]
        for spin_n in range(1,spin+1):
            for line in txt: 
                if (line.find('state #')!= -1):
                    if SOC==1:
                        tmp=line.replace('m_j=',' ').replace('j=',' ').replace('l=',' ')\
                        .replace('(',' ').replace(')',' ').replace('#',' ').replace(':',' ').split()  
                        s=dict()
                        s['label']=count
                        s['num']=int(tmp[3])
                        s['spec']=tmp[4]
                        s['wfc']=int(tmp[7])
                        s['l']=float(tmp[-3])
                        s['j']=float(tmp[-2])
                        s['mj']=float(tmp[-1])             
                        state_info.append('{label:3d} => {spec:>2s}-{num:<2d} orb-{wfc:<2d} ({num:2d} / {l} / {j:3.1f} / {mj:+3.1f})\n'.format(**s))
                    elif SOC==0:
                        tmp= line.replace('#',' ').replace(':',' ').replace('(',' ')\
                        .replace(')',' ').replace('l=',' ').replace('m=',' ').split()
                         
                        s=dict()
                        s['label']=count
                        s['num']=int(tmp[3])
                        s['spec']=tmp[4]
                        s['wfc']=int(tmp[7])
                        s['l']=int(tmp[-2])
                        s['ml']=int(tmp[-1]) 
                        s['ms']=int(np.sign(1.5-spin_n))
                        state_info.append(\
                        '{label:3d} => {spec:>2s}-{num:<2d} orb-{wfc:<2d} ({num:2d} / {l} / {ml:d} / {ms:+d})\n'.format(**s))

                    count+=1
        
        # save results 
        print('=> save fatband information in fatband.npz')
        np.savez(self.wkdir+'fatband.npz',spin=spin,Ek=Ek,Ek_weight=Ek_weight,state_info=state_info)
        
        # print state_info
        print(' => check fatband-state.dat for state_info:')
        print(np.array(state_info))
        with open(self.wkdir+'fatband-state.dat','w') as file:
            if SOC is 0:
                file.write('  # => atom  orb-#  (at / l / ml/ ms)\n')
            elif SOC is 1:
                file.write('  # => atom  orb-#  (at /   l /   j /   mj)\n')
            file.writelines(state_info)
        
    def fatband_plot(self,state_grp,kdiv='default',klabel='default',\
    Ebound='default',ini_fig_num=1,marker_size=30,colorcode='b',fontsize=18):
        print('fatband_plot start ...')
        fatband=np.load(self.wkdir+'fatband.npz')
            
        # plot band stucture
        dftpp.band_plot(self,1,fatband['Ek'],\
        0,self._kdiv_qe(kdiv),klabel,Ebound,2,fontsize)
        
        # plot fatband structure
        dftpp.fatband_plot(self,fatband['Ek'],fatband['Ek_weight'],0,\
        fatband['state_info'],self.state_grp_trans(fatband['state_info'],state_grp),\
        self._kdiv_qe(kdiv),klabel,Ebound,ini_fig_num,marker_size,colorcode,fontsize,\
        self.wkdir)
    
    def pdos_read(self,Ef):
        print('pdos_read start ...')
        # search read files
        flist=[fname for fname in os.listdir(self.wkdir) if fname.find('_atm#')!=-1]
        
        # get state info
        state_info=[]
        lorb={'s':0,'p':1,'d':2,'f':3}
        state_count=0
        for fname in flist:
            # get header
            with open(self.wkdir+fname,'r') as file:
                header=file.readline()
                if header.find('up') is not -1:
                    spin=2
                    soc=0
                elif header.find('_') is not -1:
                    spin=1
                    soc=1
                else:
                    spin=1
                    soc=0
            # get state_info
            if (soc is 0):
                s_ind=fname.replace('(',' ').replace(')',' ').replace('#',' ').split()
                s=dict()
                s['num']=int(s_ind[1])
                s['spec']=s_ind[2]
                s['wfc']=int(s_ind[4])
                s['l']=lorb[s_ind[-1]]
                for ml in range(1,2*s['l']+1+1):
                    s['ml']=ml
                    for ms in range(1,spin+1):
                        s['ms']=int(np.sign(1.5-ms))
                        s['label']=state_count
                        state_info.append(\
                        '{label:>3d} => {spec:>2s}-{num:<2d} orb-{wfc:<2d} ({num:>2d} / {l:d} / {ml:d} / {ms:+d})\n'.format(**s))
                        state_count+=1
            else:
                s_ind=fname.replace('#',' ').replace('_',' ').replace('(',' ')\
                .replace(')',' ').replace('j',' ').split()
                s=dict()
                s['num']=int(s_ind[2])
                s['spec']=s_ind[3]
                s['wfc']=int(s_ind[5])
                s['l']=lorb[s_ind[6]]
                s['j']=float(s_ind[7])
                for mj in np.linspace(-s['j'],s['j'],int(2*s['j']+1)):
                    s['mj']=mj
                    s['label']=state_count
                    state_info.append(\
                    '{label:>3d} => {spec:>2s}-{num:<2d} orb-{wfc:<2d} ({num:>2d} / {l:d} / {j:3.1f} / {mj:+3.1f})\n'.format(**s))
                    state_count+=1
                    

        # check energy div (to make pdos as a matrix with well-defined size)
        with open(self.wkdir+flist[0],'r') as file:
            flines=file.readlines()
            tot_E=len(flines[1:])
            E=np.array([float(fline_n.split()[0]) for fline_n in flines[1:]])-Ef
            

        # read pdos
        tot_pdos=len(state_info)
        pdos=np.zeros((tot_E,tot_pdos))
        pdos_count=0
        for fname in flist:
            with open(self.wkdir+fname,'r') as file:
                flines=file.readlines()
                nldos=len([colname_n for colname_n in flines[0].split() if colname_n.find('ldos')!=-1])
                npdos=len([colname_n for colname_n in flines[0].split() if colname_n.find('pdos')!=-1])              
                for n, pdos_val in enumerate(flines[1:]):
                    pdos_En=[float(val) for val in pdos_val.split()[nldos+1:]]
                    pdos[n,pdos_count:pdos_count+len(pdos_En)]=pdos_En
                # print(pdos[:,pdos_count:pdos_count+len(pdos_En)])
                pdos_count+=len(pdos_En)
        
        # save to file
        with open(self.wkdir+'pdos-state.dat','w') as file:
            if (soc is 0):
                file.write('  # => at-#  orb-#  (at / l /ml / ms)\n')
            else:
                file.write('  # => at-#  orb-#  (at / l /   j /   mj)\n')
            file.writelines(state_info)
        np.savez(self.wkdir+'pdos.npz',E=E,state_info=state_info,pdos=pdos)
        
    def pdos_plot(self,state_grp,Ebound,lw=3,fontsize=18):
        pdos=np.load(self.wkdir+'pdos.npz')
        dftpp.spectral_plot(self,pdos['E'],pdos['pdos'],pdos['state_info'],\
        self.state_grp_trans(pdos['state_info'],state_grp),\
        savefig_path=self.wkdir+'pdos.png',Ebound=Ebound,\
        lw=lw,fontsize=fontsize)

