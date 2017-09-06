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

        # check if V<0, if yes, swap a1 and a2
        V=np.dot(np.cross(a_vec[0,:],a_vec[1,:]),a_vec[2,:])
        print('system volume={0} A^3'.format(V))
        if V < 0:
            print('Warning: V<0, a1/b1 and a2/b2 swapped, else abinit will report error!')
            a_vec=a_vec[[1,0,2],:]
            sublat=sublat[:,[1,0,2]]
            kpath=kpath[:,[1,0,2]]

        kdiv=self.kdiv(a_vec,kpath,kdense)
        if type(atom[0])==str:
            atom=[ self.__grep__(ptable,atn)[0] for atn in atom]          
        # convert atomic name format to atomic number format
        spec=[]
        [spec.append(at_n) for at_n in atom if spec.count(at_n)==0]
        # -------------------------------------------------
        # output input file for abinit
        file=open(self.wkdir+'abinit.in','w')
        file.write('# << title: {0} / task: ground state >>\n'.format(prefix))
        file.write('\n')
        file.write('# System ========================================\n')
        file.write('ndtset 4         # of running dataset\n')
        file.write('jdtset 1 2 3 4   # index of running dataset\n')
        file.write('autoparal 1      # level of prallal\n')
        file.write('pawovlp 10       # allowed paw overlap,>15 dagnerous!\n')
        file.write('\n')
        file.write('# Strucure ======================================\n')
        file.write('chkprim 0   # warnin on non-primitive cell\n')
        file.write('nsym 0      # auto symmetry finder\n')
        file.write('ntypat %i    # number of species\n' % len(set(atom)))
        file.write('znucl       # Z of each species\n  ')
        [file.write('%i ' % s) for s in spec]
        file.write('\n')
        file.write('natom %i    # number of atoms\n' % len(atom))
        file.write('typat       # species of each atom\n')
        for n, atom_n in enumerate(atom):
            for m, spec_m in enumerate(spec): 
                if spec_m==atom_n:
                    file.write('%3i ' % (m+1)) 
            if (divmod(n,10)[1]==9) or (n+1==len(atom)):
                file.write('\n')
        file.write('acell 3*1.0 angstrom    # lattice constant\n')
        file.write('rprim                   # lattice vectors\n')
        for n in range(0,3):
            file.write('  %12.8f %12.8f %12.8f\n' % tuple(a_vec[n,:].tolist()))       
        file.write('xred                    # atom in reduced coordinates\n')    
        for n in range(0,len(atom)):
            file.write('  {0[0]:12.8f} {0[1]:12.8f} {0[2]:12.8f}  #  {1:>2s} ({2:>3d})\n'.\
                format(sublat[n,:],ptable[atom[n]],(n+1)))
            #file.write('%12.8f %12.8f %12.8f\n' % tuple(sublat[n,:].tolist()))
        file.write('\n')
        file.write('# Magnetism ======================================\n') 
        if mag=='on':
            file.write('# -----\n')
            file.write('# Please put the magnetic moments of each atoms in spinat.\n')
            file.write('# If you want to fix the moment during SCF, you can use\n')
            file.write('# magconon=1 for magnetic constraint calculations. The \n')
            file.write('# strength of the constraint is controlled by magcon_lambda.\n') 
            file.write('# If scf cannot be achieved, make magcon_lambda smaller. If\n')
            file.write('# the self-consistent moments are far from what you want, make\n')
            file.write('# it larger. From 0 to hundreds are all reasonable.\n')
            file.write('# -----\n\n')
        if soc is 'on':
            file.write('nspinor 2            # spinor\n')
            file.write('nspden  4            # spin density\n')
            file.write('nsppol  1            # spin polarizaton\n')
            file.write('pawspnorb 1          # turn on paw SOC\n')
            file.write('spnorbscl 1.0        # scaling of SOC\n')
        elif (mag is 'on'):
            file.write('nspinor 1            # spinor\n')
            file.write('nspden  2            # spin density\n')
            file.write('nsppol  2            # spin polarzation\n')   
        elif (mag is 'off'):
            file.write('nspinor 1       # spinor\n')
            file.write('nspden  1       # spin density\n')
            file.write('nsppol  1       # spin polarzation\n')
        
        if (mag is 'on'):
            file.write('magconon 0           # moment constraint, 0:off, 1:dir only, 2:full vec\n')
            file.write('magcon_lambda 1.0    # constraint multipler, 0 ~ hundreds\n') 
            file.write('spinat               # magmom of each atom\n')
            #[file.write('   0.0  0.0  0.0     # atom:{0:>3d}\n'.format((n+1))) for n in range(0,len(atom))]            
            [file.write('  0.0  0.0  0.0  #  {0:>2s} ({1:>3d})\n'.\
                format(ptable[atom[n]],(n+1))) for n in range(0,len(atom))]
        file.write('\n')
        if dftu=='on':
            file.write('# Correlation =====================================\n')
            file.write('# -----\n')
            file.write('# Please fill the correct U, J and orbitals. Note that\n')
            file.write('# local PBE0 (only if SOC is off) is an alternative approach\n')
            file.write('# for parameter-free DFT+U which include the exact\n')
            file.write('# exchange within the PAW sphere. So do not use simtaneously.\n')
            file.write('# -----\n\n')
            file.write('dmatpuopt  1              # output density matrix\n\n')
            file.write('# DFT+U\n')
            file.write('usepawu 0                 # DFT+U,0:off,1:FLL,2:AMF\n')
            file.write('lpawu                     # l for each species in DFT+U\n  ')
            [file.write('-1 ') for n in range(0,len(spec))]
            file.write('\n')
            file.write('upawu                     # u for each species in DFT+U\n  ')
            [file.write('0.0  ') for n in range(0,len(spec))]
            file.write('eV \n')
            file.write('jpawu                     # j for each species in DFT+U\n  ')
            [file.write('0.0  ') for n in range(0,len(spec))]
            file.write('eV \n\n') 
            if soc=='off':
                file.write('# Local PBE0\n')
                file.write('useexexch 0               # Local PBE0, 0:off, 1:on\n')
                file.write('exchmix 0.25              # PBE0 coeff, default: 0.25\n')
                file.write('lexexch                   # orbitals to apply local PBE0\n')
                [file.write('-1 ') for n in range(0,len(spec))]
                file.write('\n\n')
        file.write('# Plane Waves ================================\n')
        file.write('ecut 16.5      # energy cut in Ha\n')
        file.write('pawecutdg 41.5 # ecut double grid, suggest: 2.5~3.0 ecut\n')
        file.write('occopt 7       # occupation, 7 should work for all\n')
        file.write('diemac 1e+6    # dielectric const, 1e+6 should work for all\n')
        file.write('\n')
        file.write('# DS1: SCF ===================================\n')
        file.write('# -----\n')
        file.write('# Abinit requires kmesh to satisify the symmetry. The simplest\n')
        file.write('# way to do it is to let Abinit automatically look for kmesh by\n')
        file.write('# constructing a supercell in real space. To do it, first run\n')
        file.write('# kptrlen & prtkpt to search for best supercell kptrlatt.\n')
        file.write('# put the best values of kptrlatt, nshiftk, shiftk from\n') 
        file.write('# abinit.out. erase kptrlen & prtkpt. Then run scf with correct\n')
        file.write('# kptrlatt, nshiftk, shiftk. \n')
        file.write('# -----\n\n')
        file.write('nstep1 100          # max scf steps\n')
        file.write('tolvrs1 1e-3        # convergence of potential\n\n')
        file.write('# run this part first, then erase them after kptrlatt is obtained\n')
        file.write('kptrlen1 50          # search for best kptrlatt\n')
        file.write('prtkpt1 1\n\n')
        file.write('# then run this part with correct paramaters\n')
        file.write('#kptrlatt1            # kmesh in real space (obtained by kptrlen)\n')
        file.write('#   0  0  0\n')
        file.write('#   0  0  0\n')
        file.write('#   0  0  0\n')
        file.write('#nshiftk1 1          # number of k shift\n')
        file.write('#shiftk1             # k-shift\n')
        file.write('#   0.0  0.0  0.0\n\n')
        file.write('\n')
        file.write('#DS2: Band Calculation ==========================\n')
        file.write('iscf2 -2           # -2: nscf run\n')
        file.write('getden2 1          # get charge density from n-th dataset\n')
        file.write('tolwfr2 1.0d-8     # nscf wf convergence\n')
        file.write('kptopt2 -%i\n' % len(kdiv))
        file.write('ndivk2\n  ')
        [file.write('{0:d}  '.format(val)) for val in kdiv] 
        file.write('\nkptbounds2\n')
        [file.write('{0[0]:12.8f}  {0[1]:12.8f}  {0[2]:12.8f}\n'.format(kpath_n)) for kpath_n in kpath]
        file.write('\n')
        file.write('#DS3: Fatband Calculation =======================\n')
        file.write('# -----\n')
        file.write('# To analyze fatband using DFTtoolbox, keep pawfatband=2\n')
        file.write('# to get higest resolution up to (L,M).\n')
        file.write('# -----\n\n')
        file.write('iscf3 -2           # -2: nscf run\n')
        file.write('getden3 1          # get charge density from n-th dataset\n')
        file.write('tolwfr3 1.0d-8     # nscf wf convergence\n')
        file.write('kptopt3 -%i\n' % len(kdiv))
        file.write('ndivk3\n  ')
        [file.write('{0:d}  '.format(val)) for val in kdiv] 
        file.write('\nkptbounds3\n')
        [file.write('{0[0]:12.8f}  {0[1]:12.8f}  {0[2]:12.8f}\n'.format(kpath_n)) for kpath_n in kpath]   
        file.write('pawfatbnd3 2       # fat band (L,M), DFTtoolbox default, dont change !\n')
        file.write('natsph3 {0:3d}        # total projected atom\n'.format(len(atom)))
        file.write('iatsph3            # indexs of projected atoms\n  ') 
        for n in range(0,len(atom)):
            file.write('{0:>2d}  '.format((n+1)))
            if divmod(n+1,10)[1] is 0:
                file.write('\n  ')
        file.write('\n\n')
        file.write('#DS4: PDOS Calculation ==========================\n')
        file.write('# -----\n')
        file.write('# To analyze PDOS using DFTtoolbox, keep prtdos=3,\n')
        file.write('# prtdosm=2, pawprtdos=2 to get highest resolution up to (L,M)\n')
        file.write('# -----\n\n')
        file.write('iscf4 -3            # -3: dos nscf run if occopt=3~7\n')
        file.write('getden4 1           # get charge density from n-th dataset\n')
        file.write('tolwfr4 1.0d-8      # nscf wf convergence\n')
        file.write('prtdos4 3           # calculate projected dos, must 3 for DFTtoolbox\n')
        file.write('prtdosm4 2          # print pdos to m-component, must 1 or 2 for DFTtoolbox\n')
        file.write('pawprtdos4 2        # calculate PAW part only, must 2 for DFTtoolbox\n')
        file.write('natsph4 {0:3d}         # total projected atom\n'.format(len(atom)))
        file.write('iatsph4             # indexs of projected atoms\n  ') 
        for n in range(0,len(atom)):
            file.write('{0:>2d}  '.format((n+1)))
            if (divmod(n+1,10)[1] is 0):
                file.write('\n  ')
        file.write('\n# For following variables, please copy from DS1 !!!!!!!!!!!\n')
        file.write('#kptrlatt4            # kmesh in real space (obtained by kptrlen)\n')
        file.write('#   0  0  0\n')
        file.write('#   0  0  0\n')
        file.write('#   0  0  0\n')
        file.write('#nshiftk4 1          # number of k shift\n')
        file.write('#shiftk4             # k-shift\n')
        file.write('#   0.0  0.0  0.0\n\n')
        file.close()
        
        with open(self.wkdir+'abinit.file','w') as file:
            file.write('abinit.in\n')
            file.write('abinit.out\n')
            file.write('abinit_i\n')
            file.write('abinit_o\n')
            file.write('abinit_t\n')
            [file.write('/PSP_DIR/{0}-PSP_FILE_NAME\n'.format(ptable[s])) for s in spec]
        print('abinit.in & abinit.file are generated, please check all varaibles!')
        print('To run it, just type:')
        print('  $ > abinit < abinit.file > abinit.log')
    def relax(self,prefix):
        ptable=self.ptable()
        #spec=sorted(tuple(set(atom)))         
        atom, a_vec, sublat=self.getxsf(self.wkdir,prefix)
        if type(atom[0])==str:
            atom=[ self.__grep__(ptable,atn)[0] for atn in atom]     
        # check if V<0, if yes, swap a1 and a2
        V=np.dot(np.cross(a_vec[0,:],a_vec[1,:]),a_vec[2,:])
        print('system volume={0} A^3'.format(V))
        if V < 0:
            print('Warning: V<0, a1/b1 and a2/b2 swapped, else abinit will report error!')
            a_vec=a_vec[[1,0,2],:]
            sublat=sublat[:,[1,0,2]]    
        # convert atomic name format to atomic number format
        spec=[]
        [spec.append(at_n) for at_n in atom if spec.count(at_n)==0]

        # output input file for abinit
        file=open(self.wkdir+'abinit.in','w')
        file.write('# << title: {0} / task: relax >>\n'.format(prefix))
        file.write('\n')
        file.write('# System ========================================\n')
        file.write('ndtset 2         # of running dataset\n')
        file.write('jdtset 1 2       # index of running dataset\n')
        file.write('autoparal 0      # level of prallal\n')
        file.write('pawovlp 10       # allow paw overlap,>15 dagnerous!\n')
        file.write('\n')
        file.write('# Strucure ======================================\n')
        file.write('chkprim 0   # warnin on non-primitive cell\n')
        file.write('nsym 0      # auto symmetry finder\n')
        file.write('ntypat %i    # number of species\n' % len(set(atom)))
        file.write('znucl       # Z of each species\n  ')
        [file.write('%i ' % s) for s in spec]
        file.write('\n')
        file.write('natom %i    # number of atoms\n' % len(atom))
        file.write('typat       # species of each atom\n')
        for n, atom_n in enumerate(atom):
            for m, spec_m in enumerate(spec): 
                if spec_m==atom_n:
                    file.write('%3i ' % (m+1)) 
            if (divmod(n,10)[1]==9) or (n+1==len(atom)):
                file.write('\n')
        file.write('acell 3*1.0 angstrom    # lattice constant\n')
        file.write('rprim                   # lattice vectors\n')
        for n in range(0,3):
            file.write('  %12.8f %12.8f %12.8f\n' % tuple(a_vec[n,:].tolist()))       
        file.write('xred                    # atom in reduced coordinates\n')    
        for n in range(0,len(atom)):
            file.write('  {0[0]:12.8f} {0[0]:12.8f} {0[0]:12.8f}  #  {1:>2s} ({2:>3d})\n'.\
                format(sublat[n,:],ptable[atom[n]],(n+1)))
        file.write('\n')
        file.write('# Plane Waves ================================\n')
        file.write('ecut 16.5      # energy cut in Ha\n')
        file.write('pawecutdg 41.5 # ecut double grid, suggest: 2.5~3.0 ecut\n')
        file.write('occopt 7       # occupation, 7 should work for all\n')
        file.write('diemac 1e+6    # dielectric const, 1e+6 should work for all\n')
        file.write('tolvrs 1e-3    # potential convergence for scf\n')
        file.write('\n')
        file.write('# Kmesh ======================================\n')   
        file.write('# -----\n')
        file.write('# Abinit requires kmesh to satisify the symmetry. The simplest\n')
        file.write('# way to do it is to let Abinit automatically look for kmesh by\n')
        file.write('# constructing a supercell in real space. To do it, first run\n')
        file.write('# kptrlen & prtkpt to search for best supercell kptrlatt.\n')
        file.write('# put the best values of kptrlatt, nshiftk, shiftk from\n') 
        file.write('# abinit.out. erase kptrlen & prtkpt. Then run scf with correct\n')
        file.write('# kptrlatt, nshiftk, shiftk. \n')
        file.write('# -----\n\n')   
        file.write('# run it frist to obtain kptrlatt, then comment out\n')
        file.write('kptrlen1 50          # search for best kptrlatt\n')
        file.write('prtkpt1 1\n\n')
        file.write('# uncomment the following part after kptrlatt is obtained\n')
        file.write('#kptrlatt1            # kmesh in real space (obtained by kptrlen)\n')
        file.write('#   0  0  0\n')
        file.write('#   0  0  0\n')
        file.write('#   0  0  0\n')
        file.write('#nshiftk 1       # how many k-shift\n')
        file.write('#shiftk           # shift vector\n')
        file.write('#    0.0 0.0 0.0\n')
        file.write('\n')
        file.write('# DS1: In-Cell Relax =============================\n')
        file.write('-----\n')
        file.write('# It is strongly recommend to use a two-step optimization to\n')
        file.write('# make system convergent easier. First keep the cell fixed\n')
        file.write('# and relax atoms within the cell. Second, relax the cell using\n')
        file.write('# optcell=1 (lattice constants only) or, more aggresively,\n')
        file.write('# optcell=2 (relax the whole system).\n')
        file.write('-----\n\n')
        file.write('optcell1 0       # 0: cell fixed, 1:acell only, 2:full relax\n')
        file.write('ionmov1 3        # BFGS algorithm, work for must cases\n')
        file.write('ntime1 30        # max optimization steps\n')
        file.write('tolmxf1 1e-4     # force convergence\n')
        file.write('ecutsm1 0.5      # total energy smooth, enlarge if change sharply\n')
        file.write('dilatmx1 1.0     # estimated lattice expansion, enlarge if charge sharply\n')
        file.write('\n')
        file.write('# DS2: Cell Relax =================================\n')
        file.write('getcell2 1       # get cell from n-th dataset\n')
        file.write('getxred2 1       # get xred from n-th dataset\n')
        file.write('optcell2 1       # 0: cell fixed, 1:acell only, 2:full relax\n')
        file.write('ionmov2 3        # BFGS algorithm, work for must cases\n')
        file.write('ntime2 30        # max optimization steps\n')
        file.write('tolmxf2 1e-4     # force convergence\n')
        file.write('ecutsm2 0.5      # total energy smooth, enlarge if change sharply\n')
        file.write('dilatmx2 1.1     # estimated lattice expansion, enlarge if expand sharply\n')
        file.close()
        with open(self.wkdir+'abinit.file','w') as file:
            file.write('abinit.in\n')
            file.write('abinit.out\n')
            file.write('abinit_i\n')
            file.write('abinit_o\n')
            file.write('abinit_t\n')
            [file.write('/PSP_DIR/{0}-PSP_FILE_NAME\n'.format(ptable[s])) for s in spec]
        print('abinit.in & abinit.file are generated, please check all varaibles!')
        print('To run it, just type:')
        print('  $ > abinit < abinit.file > abinit.log')        
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
    __kdiv_qe__(self,kdiv)
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
        
    def __kdiv_conv__(self,kdiv):
        # convert kdiv from abinit format to postprocess format
        if kdiv=='default':
            kdiv_conv='default'
        else:
            kdiv_conv=[ sum(kdiv[0:n+1]) for n,val in enumerate(kdiv)]
            kdiv_conv.insert(0,0)
            kdiv_conv[-1]=kdiv_conv[-1]-1            
        
        return kdiv_conv
        
        
    def band_read(self,Ef,dataset,spin=1):
        print('band_read start ...')
        flist=os.listdir(self.wkdir)
        fn=self.grep(flist,'_DS'+str(dataset)+'_EIG')
        [fn.pop(n) for n, val in enumerate(fn) if flist[val].find('.nc')!=-1]
        
        if len(fn)!=1:
            print('Error: '+self.wkdir+' contains no/more DS'+str(dataset)+'_EIG files!')
            sys.exit()
        else:
            with open(self.wkdir+flist[fn[0]],'r') as file:
                txt=file.readlines()
        
        # read total bands and total k-points
        tot_k=int(txt[0][32:36])
        tot_ban=int(txt[1][17:20])
        
        spn_ln=self.grep(txt,'SPIN')
        if len(spn_ln)==2:
            spin=2
            # delete header lines
            txt.pop(spn_ln[0])
            txt.pop(spn_ln[1]-1)
        else:
            spin=1
            # delete header lines
            txt.pop(0)
            
        # count size of the each band data block (include k-point)
        ban_mod8=divmod(tot_ban,8)
        if ban_mod8[1]==0:
            block_lines=ban_mod8[0]+1
        else:
            block_lines=ban_mod8[0]+1+1
        
        # read band data
        k_point=np.zeros((spin*tot_k,3))
        Ek=np.zeros((spin*tot_k,tot_ban))
        
        for n in range(0,spin*tot_k):
            k_point[n,:]=[float(val) for val in txt[n*block_lines][41:65].split()]
            
            ban_val=[]
            for val in txt[n*block_lines+1:(n+1)*block_lines]:                
                ban_val+=val.split()
            
            Ek[n,:]=[float(val) for val in ban_val]
        
        # reshape k_point and Ek to correct size        
        if spin==2:
            k_point=k_point[0:tot_k,:]
            Ek=np.concatenate((Ek[0:tot_k,:],Ek[tot_k:,:]),1)
        
        Ek=(Ek-Ef)*27.2
        np.savez(self.wkdir+'band-DS'+str(dataset)+'.npz',Ek=Ek,k_point=k_point,spin=spin)
        print(' => band data saved to band-DS'+str(dataset)+'.npz')
        
    def band_plot(self,dataset,kdiv='default',klabel='default',Ebound='default',lw=2,fontsize=18):
        print('band_plot start ...')
        band=np.load(self.wkdir+'band-DS'+str(dataset)+'.npz')
        
        # abinit kdiv has a bug in band calculation (not in fatband). 
        # the last segment has one more k-point
        if kdiv!='default':
            kdiv[-1]=kdiv[-1]+1
            
        dftpp.band_plot(self,band['spin'].tolist(),band['Ek'],\
        0,self.__kdiv_conv__(kdiv),klabel,Ebound,lw,fontsize,self.wkdir)
        
        shutil.move(self.wkdir+'band.png',self.wkdir+'band-DS'+str(dataset)+'.png')
        print(' => band plot saved to band-DS'+str(dataset)+'.png')
        
    def fatband_read(self,dataset):
        print('fatband_read start ...')
        flist=[fname for fname in os.listdir(self.wkdir) if fname.find('DS'+str(dataset)+'_FATBANDS_') is not -1]
        flist.sort()
        # get state_info from file names =====================
        state_info=[]
        for n, fname in enumerate(flist):
            label=fname[fname.find('_at')+3:].replace('_',' ').split()
            s=dict()
            s['label']=n
            s['num']=int(label[0])
            s['name']=label[1]
            s['ms']=int(np.sign(1.5-int(label[2][-1])))
            s['l']=int(label[3][-1])
            s['ml']=int(label[4][-2:])
            state_info.append('{label:2d} => {name:2s} ({num:3d} / {l} / {ml:+} / {ms:+})\n'.format(**s))
        
        # check whether size of the  calculation ==========
        if self.grep(flist,'_is2_')!=[]:
            spin=2
        else:
            spin=1
            
        with open(self.wkdir+flist[0]) as file:
                txt=file.readlines()
                
        #line numbers where ban data begins
        bn=self.grep(txt,'# BAND number') 
        
        # define size parameters
        tot_ban=len(bn) # band number w/o spin 
        tot_k=bn[1]-bn[0]-3
        tot_state=len(state_info) # tot state include spin
        
        # Read eigenvalues =====================================        
        Ek=np.zeros((tot_k,spin*tot_ban))
        for spin_n in range(1,spin+1):
            with open(self.wkdir+flist[self.grep(flist,'_is'+str(spin_n)+'_')[0]]) as file:
                txt=file.readlines()
            
            # read eigenvalues 
            print(' => reading Ek')
            for ban_ind,ban_ln in enumerate(bn):
                Ek[:,ban_ind+tot_ban*(spin_n-1)]=\
                [float(line.split()[1]) for line in txt[ban_ln+1:ban_ln+tot_k+1]]
        
        # read fatband data ====================================
        print(' => reading Ek_weight')
        Ek_weight=np.zeros((tot_k,spin*tot_ban,tot_state))
        print('     tot_state=%i' % tot_state)
        print('     ',end='')
        for state_n,fname in enumerate(flist):
            # show reading information
            if divmod(state_n,10)[1]==9:
                print('%4i ' % state_n)
                print('     ',end='')
            elif state_n==tot_state-1:
                print('%4i ' % state_n)
            else:
                print('%4i ' % state_n,end='')
                
            with open(self.wkdir+fname) as file:
                txt=file.readlines()
            
            # read fatband data
            for ban_n, ban_ln in enumerate(bn):
                if fname.find('_is1_')!=-1:
                    Ek_weight[:,ban_n,state_n]=\
                    [float(line.split()[2]) for line in txt[ban_ln+1:ban_ln+tot_k+1]]
                elif fname.find('_is2_')!=-1:
                    Ek_weight[:,ban_n+tot_ban,state_n]=\
                    [float(line.split()[2]) for line in txt[ban_ln+1:ban_ln+tot_k+1]]
 
        # save results ===========================================
        print(' => save band weight and state info as fatband-DS'+str(dataset)+'.npz')
        np.savez(self.wkdir+'fatband-DS'+str(dataset)+'.npz',Ek=Ek,Ek_weight=Ek_weight,state_info=state_info)
        
        # print state_info
        print(' => check fatband-DS'+str(dataset)+'-state.dat for state_info')
        with open(self.wkdir+'fatband-DS'+str(dataset)+'-state.dat','w') as file:
            file.write(' # => at (atn / l / ml / ms)\n')
            file.writelines(state_info)
        print(np.array(state_info))
        
    def fatband_plot(self,dataset,state_grp,kdiv='default',klabel='default',\
    Ebound='default',ini_fig_num=1,marker_size=30,colorcode='b',fontsize=18):
        print('fatband_plot start ...')
        fatband=np.load(self.wkdir+'fatband-DS'+str(dataset)+'.npz')
        
        # plot band stucture
        dftpp.band_plot(self,1,fatband['Ek'],\
        0,self.__kdiv_conv__(kdiv),klabel,Ebound,2,fontsize)    
        
        # plot fatband     
        dftpp.fatband_plot(self,fatband['Ek'],fatband['Ek_weight'],0,\
        fatband['state_info'],self.state_grp_trans(fatband['state_info'],state_grp),\
        self.__kdiv_conv__(kdiv),klabel,Ebound,ini_fig_num,marker_size,colorcode,fontsize\
        ,self.wkdir)
        
        # rename to output png file to abinit filename
        for n in range(0,len(state_grp)):
            shutil.move(self.wkdir+'fatband-'+str(n)+'.png',\
            self.wkdir+'fatband-DS'+str(dataset)+'-'+str(n)+'.png')   
            
        print(' => fatband plot saved to fatband-DS'+str(dataset)+'-x.png')    

    def pdos_read(self,dataset):
        # Note: only can read pawprtdos=0/2, prtdosm=0
        print('pdos_read start ...')
        flist=[ fname for fname in os.listdir(self.wkdir) if fname.find('_o_DS{0}_DOS_AT'.format(dataset))!=-1]
        tot_atom=len(flist)-1
        # read header
        E_range=[]
        Ef=[]
        tot_E=0
        spin=1
        with open(self.wkdir+flist[0],'r') as file:
            for n in range(0,30):
                line=file.readline()
                if line.find('between')!=-1:
                    E_range=[float(line.split()[2]),float(line.split()[4]),float(line.split()[9])]
                elif line.find('Spin-up DOS')!=-1:
                    spin=2
                elif line.find('Fermi energy')!=-1:
                    Ef=float(line.split()[-1]) 
                elif line.find('covering the interval ')!=-1:
                    tot_E=int(line.split()[2])
        if (E_range==[]) or (Ef==[]) or (tot_E==0):
            print('Error: header read error!')
            sys.exit()
        else:
            print(' => Ef={0}, E_range={1}~{2}, tot_E={3}'.format(Ef,E_range[0],E_range[1],tot_E))
        
        # generate E
        E=(np.linspace(E_range[0],E_range[1],tot_E)-Ef)*27.2 
        
        # read PDOS, note: integral DOS means accumulated DOS from the band bottom
        state_info=[]
        state_count=0
        for atn, fname in enumerate(flist):
            # construct state info
            s=dict()
            s['num']=int(fname[-4:])
            for ms in range(1,spin+1):
                s['ms']=int(np.sign(1.5-ms))
                for l in range(0,4):
                    s['l']=l
                    for ml in range(-l,+l+1):
                        s['ml']=ml
                        s['label']=state_count                   
                        state_info.append('{label:4d} => ({num:3d} / {l} / {ml:+} / {ms:+})\n'.format(**s)) 
                        state_count+=1
                
            with open(self.wkdir+fname,'r') as file:
                flines=file.readlines()
            
            data_start=self.grep(flines,'# energy(Ha)')
            for spin_n in range(1,spin+1):
                data=np.array([val.split()[11:27] for val in \
                flines[data_start[spin_n-1]+1:data_start[spin_n-1]+tot_E+1] ]).astype(np.float)
                if 'pdos' in locals():
                    pdos=np.concatenate((pdos,data),1)
                else:
                    pdos=data
        np.savez(self.wkdir+'pdos-DS'+str(dataset)+'.npz',E=E,pdos=pdos,state_info=state_info) 
        with open(self.wkdir+'pdos-DS'+str(dataset)+'-state.dat','w') as file:
            file.write('   # => (atn / l / ml / ms)\n')
            file.writelines(state_info)
            print(np.array(state_info))
            
    def pdos_plot(self,dataset,state_grp,Ebound,lw=3,fontsize=18):
        print('starting pdos_plot ...')

        pdos=np.load(self.wkdir+'pdos-DS'+str(dataset)+'.npz')
        dftpp.spectral_plot(self,pdos['E'],pdos['pdos'],pdos['state_info'],\
        self.state_grp_trans(pdos['state_info'],state_grp),\
        savefig_path=self.wkdir+'pdos-DS'+str(dataset)+'.png',Ebound=Ebound,\
        lw=lw,fontsize=fontsize)
        
