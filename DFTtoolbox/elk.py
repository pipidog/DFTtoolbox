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
        atom, a_vec, sublat=self.getxsf(self.wkdir,prefix,'on')
        klabel, kpath=self.getkpf(self.wkdir,prefix) 
        kdiv=self.kdiv(a_vec,kpath,kdense)
        if type(atom[0])==str:
            atom=[ self.__grep__(ptable,atn)[0] for atn in atom]          
        # convert atomic name format to atomic number format
        spec=[]
        [spec.append(at_n) for at_n in atom if spec.count(at_n)==0]
        # -------------------------------------------------
        # begin generate input files
        file=open(self.wkdir+'elk.in','w')
        file.write('! << title:{0} / task: ground state>>\n\n'.format(prefix))
        file.write('! system ===================\n')
        file.write('tasks\n')
        file.write('   0              ! scf\n')
        file.write('  10              ! pdos\n')
        file.write('  20              ! band\n')
        file.write('  21              ! fatbnd\n\n')
        file.write('sppath          ! species file path\n')
        file.write('  \'/YOUR-ELK-ROOT/species/\'\n\n')
        file.write('xctype          ! exchange functional (3:LSDA, 20:PBE, 22:PBESOL)\n')
        file.write('  20 0 0\n\n')
        file.write('! structure ===============\n')
        file.write('primcell       ! auto find primitive cell for calculation\n')
        file.write('  .true.\n\n')
        file.write('scale          ! lattice const\n')
        file.write('  1.89\n\n')
        file.write('avec           ! lattice vectors\n')
        [file.write('  {0[0]:12.8f} {0[1]:12.8f} {0[2]:12.8f}\n'.format(a_n)) for a_n in a_vec]
        file.write('\n')
        file.write('atoms          ! atom infomation\n')
        file.write(   '%3i          :nspecies\n' % len(spec))
        for spec_n in spec:
            sublat_label=np.nonzero(np.array(atom)==spec_n)[0]
            file.write('  %2s.in      :spfname\n' % ptable[spec_n])
            if mag is 'on':
                file.write('  %2i         :natoms; atposl, bfcmt below\n' % len(sublat_label))
                for label in sublat_label:
                    file.write('%12.8f %12.8f %12.8f   0.0   0.0   0.0\n' % tuple(sublat[label,:].tolist())) 
            else:
                file.write('  %2i         :natoms; atposl below\n' % len(sublat_label))
                for label in sublat_label:
                    file.write('%12.8f %12.8f %12.8f \n' % tuple(sublat[label,:].tolist())) 
        file.write('\n')
        if (mag is 'on') or (soc is 'on'):
            file.write('! spin ================\n')
        if mag is 'on':
            file.write('spinpol        ! spin polarization\n') 
            file.write('  .true.\n\n')
            file.write('reducebf       ! reduction of MT b-field in each iteration\n')
            file.write('  0.5\n\n')
            file.write('cmagz          ! fix moment in z-axis\n')
            file.write('  .true.\n\n') 
        if soc is 'on':
            file.write('spinorb        ! spin-orbit coupling\n')
            file.write('  .true.\n\n')
            file.write('socscf         ! scaling of soc\n') 
            file.write('   1.0\n\n')
        if (dftu is 'on'):           
            file.write('! DFT+U ================\n')
            file.write('! hints:\n')
            file.write('! dftu => 0:off, 1:FLL, 2:AMF, FLL should work for most cases\n')
            file.write('! inpdftu => keep 1 for using U and J as input\n')
            file.write('! is => species number, l => orbital, U, J => in Ha\n\n') 
            file.write('dft+u\n')
            file.write('  1   1             : dftu, inpdftu\n')
            file.write('  0   0  0.0  0.0   : is, l, U, J\n\n')
        file.write('! kmesh =================\n')
        file.write('ngridk       ! BZ grid\n  ')
        [file.write('{0} \n'.format(int(np.round(45/la.norm(a_n))))) for a_n in a_vec]
        file.write('\n')
        file.write('vkloff       ! shift of grid\n')
        file.write('  0.5 0.5 0.5\n\n')
        file.write('nempty       ! add few empty states to imporve convergence\n')
        file.write('  10\n\n')
        file.write('! convergence ===========\n')
        file.write('epsengy      ! energy criterion (default:1e-4)\n')
        file.write('  1e-4\n\n')
        file.write('maxscl       ! max scf loop\n')
        file.write('  100\n\n')
        file.write('! band structure ========\n')
        file.write('! band path: ')
        [file.write('{0}-'.format(kl)) for kl in klabel] 
        file.write('\n')
        file.write('plot1d\n')
        file.write('    %i  %i : nvp1d, npp1d \n'% (kpath.shape[0],sum(kdiv)))
        for n, kpath_n in enumerate(kpath):
            if n==0:
                file.write('{0[0]:12.8f} {0[1]:12.8f} {0[2]:12.8f} {1}\n'.format(kpath_n,':vlv1d'))              
            else:
                file.write('{0[0]:12.8f} {0[1]:12.8f} {0[2]:12.8f}\n'.format(kpath_n))
        file.close()
        
    def relax(self,prefix):        
        # collecing all input infomation ----------------
        ptable=self.ptable()
        #spec=sorted(tuple(set(atom)))         
        atom, a_vec, sublat=self.getxsf(self.wkdir,prefix,'on')
        if type(atom[0])==str:
            atom=[ self.__grep__(ptable,atn)[0] for atn in atom]          
        # convert atomic name format to atomic number format
        spec=[]
        [spec.append(at_n) for at_n in atom if spec.count(at_n)==0]
        # begin generate input files
        file=open(self.wkdir+'elk.in','w')
        file.write('! << title:{0} / task: relaxation >>\n\n'.format(prefix))
        file.write('! system ===================\n')
        file.write('tasks\n')
        file.write('   2              ! relax\n\n')
        file.write('sppath          ! species file path\n')
        file.write('  \'/YOUR-ELK-ROOT/species/\'\n\n')
        file.write('xctype          ! exchange functional (3:LSDA, 20:PBE, 22:PBESOL)\n')
        file.write('  20 0 0\n\n')
        file.write('! structure ===============\n')
        file.write('primcell       ! auto find primitive cell for calculation\n')
        file.write('  .true.\n\n')
        file.write('scale          ! lattice const\n')
        file.write('  1.89\n\n')
        file.write('avec           ! lattice vectors\n')
        [file.write('  {0[0]:12.8f} {0[1]:12.8f} {0[2]:12.8f}\n'.format(a_n)) for a_n in a_vec]
        file.write('\n')
        file.write('atoms          ! atom infomation\n')
        file.write(   '%3i          :nspecies\n' % len(spec))
        for spec_n in spec:
            sublat_label=np.nonzero(np.array(atom)==spec_n)[0]
            file.write('  %2s.in      :spfname\n' % ptable[spec_n])
            file.write('  %2i         :natoms; atposl below\n' % len(sublat_label))
            for label in sublat_label:
                file.write('%12.8f %12.8f %12.8f \n' % tuple(sublat[label,:].tolist())) 
        file.write('\n')    
        file.write('! kmesh =================\n')
        file.write('ngridk       ! BZ grid\n  ')
        [file.write('{0} '.format(int(np.round(45/la.norm(a_n))))) for a_n in a_vec]
        file.write('\n')
        file.write('vkloff       ! shift of grid\n')
        file.write('  0.5 0.5 0.5\n\n')
        file.write('nempty       ! add few empty states to imporve convergence\n')
        file.write('  10\n\n')
        file.write('! convergence ===========\n')
        file.write('epsengy      ! energy criterion (default:1e-4)\n')
        file.write('  1e-4\n\n')
        file.write('maxscl       ! max scf loop\n')
        file.write('  100\n\n')
        file.write('! relaxation ================\n')
        file.write('latvopt         ! 1: full relax, 2: fix Volume\n')
        file.write('  1\n\n')
        file.write('epsforce     ! criterion for force convergence\n')
        file.write('  1e-3\n\n')
        file.close()
        
class postproc(dftpp):
    def __init__(self,wkdir):
        if (wkdir[-1]!='/') or (wkdir[-1]!='\\'):
            wkdir+='/'
        self.wkdir=wkdir
               
    def band_read(self):
        # read bandline
        with open(self.wkdir+'BANDLINES.OUT') as file:
            banvort=[float(fline.split()[0]) for n, fline in enumerate(file) if divmod(n+1,3)[1]==1]

        # extract kdiv
        kdiv=[]
        with open(self.wkdir+'BAND.OUT') as file:
            for n, fline_n in enumerate(file):
                klen=float(fline_n.split()[0])
                if klen in set(banvort):
                    kdiv.append(n)
                    if klen==banvort[-1]:
                        break   
        tot_k=kdiv[-1]+1
        
        # read band data
        with open(self.wkdir+'BAND.OUT') as file:    
            flines=file.readlines()
            
        tot_ban=int(len(flines)/(kdiv[-1]+2))
        Ek=np.zeros((tot_k,tot_ban))
        for ban_n in range(0,tot_ban):
            data=flines[ban_n*(tot_k+1):(ban_n+1)*(tot_k+1)-1]
            # check read results
            if (float(data[0].split()[0])!=banvort[0]) or (float(data[-1].split()[0])!=banvort[-1]):
                print('Error: band read mismatch!')
                sys.exit()
            else:
                Ek[:,ban_n]=[ float(val.split()[1]) for val in data]
        Ek=Ek*27.2        
        np.savez(self.wkdir+'band.npz',Ek=Ek,kdiv=kdiv)
        
    def band_plot(self,klabel,Ebound,lw=2,fontsize=18):
        print('band_plot start ...')
        band=np.load(self.wkdir+'band.npz')
        
        dftpp.band_plot(self,1,band['Ek'],\
        0,band['kdiv'],klabel,Ebound,lw,fontsize,self.wkdir)
        
    def fatband_read(self):
        print('fatabnd_read start ...')
        flist=[fname for fname in os.listdir(self.wkdir) if fname.find('BAND_S') is not -1]
        tot_at=len(flist)
        
        # read bandline
        with open(self.wkdir+'BANDLINES.OUT') as file:
            banvort=[float(fline.split()[0]) for n, fline in enumerate(file) if divmod(n+1,3)[1]==1]

        # extract kdiv
        kdiv=[]
        with open(self.wkdir+flist[0]) as file:
            for n, fline_n in enumerate(file):
                klen=float(fline_n.split()[0])
                if klen in set(banvort):
                    kdiv.append(n)
                    if klen==banvort[-1]:
                        break   
        tot_k=kdiv[-1]+1
        
        # read band data
        with open(self.wkdir+flist[0]) as file:    
            flines=file.readlines()
            
        tot_ban=int(len(flines)/(kdiv[-1]+2))
        Ek=np.zeros((tot_k,tot_ban))
        for ban_n in range(0,tot_ban):
            data=flines[ban_n*(tot_k+1):(ban_n+1)*(tot_k+1)-1]
            # check read results
            if (float(data[0].split()[0])!=banvort[0]) or (float(data[-1].split()[0])!=banvort[-1]):
                print('Error: band read mismatch!')
                sys.exit()
            else:
                Ek[:,ban_n]=[ float(val.split()[1]) for val in data]
        Ek=Ek*27.2          
        
        # read species from input file
        with open(self.wkdir+'elk.in','r') as file:
            spec=[line.replace('.',' ').replace(':',' ').split()[0] \
            for line in file if line.find(':spfname') is not -1]

        # read fatband Weight 
        Ek_weight=np.zeros((tot_k,tot_ban,tot_at*4))
        state_info=[]
        state_count=0
        for atn, fname in enumerate(flist):
            # generate state_info
            site=[spec[int(fname[6:8])-1],int(fname[10:14])]
            for l in range(0,4):
                state_info.append(\
                '{0:>3d} => {1:>2s}{2:<} ({3:>3d} / {4} / a / a)\n'.\
                format(state_count,site[0],site[1],atn+1,l))
                state_count+=1
                
            with open(self.wkdir+fname) as file:    
                flines=file.readlines()
                
            for ban_n in range(0,tot_ban):
                data=flines[ban_n*(tot_k+1):(ban_n+1)*(tot_k+1)-1]
                # check read results
                if (float(data[0].split()[0])!=banvort[0]) or (float(data[-1].split()[0])!=banvort[-1]):
                    print('Error: band read mismatch!')
                    sys.exit()
                else:                
                    tmp=np.array([val.split()[3:] for val in data]).astype(np.float)
                    for l in range(0,4):
                        Ek_weight[:,ban_n,atn*4+l]=tmp[:,l]
                        
        print(np.array(state_info))
        with open(self.wkdir+'fataband-state.dat','w') as file:
            file.write('  # => atn ( at / l / ml / ms)\n')
            file.writelines(state_info)
        np.savez(self.wkdir+'fatband.npz',Ek=Ek,kdiv=kdiv,Ek_weight=Ek_weight,state_info=state_info)
        
    def fatband_plot(self,state_grp,klabel='default',\
    Ebound='default',ini_fig_num=1,marker_size=30,colorcode='b',fontsize=18):
        print('fatband_plot start ...')
        fatband=np.load(self.wkdir+'fatband.npz')
        
        # plot band stucture
        dftpp.band_plot(self,1,fatband['Ek'],0,fatband['kdiv'],klabel,Ebound,2,fontsize)

        # plot fatband structure
        dftpp.fatband_plot(self,fatband['Ek'],fatband['Ek_weight'],0,\
        fatband['state_info'],self.state_grp_trans(fatband['state_info'],state_grp),\
        fatband['kdiv'],klabel,Ebound,\
        ini_fig_num,marker_size,colorcode,fontsize,self.wkdir)
        
    def pdos_read(self):
        print('starting pdos_read')
        flist=[fname for fname in os.listdir(self.wkdir) if fname.find('PDOS_S') is not -1]
        tot_at=len(flist)
        
        # read E
        E=[]
        with open(self.wkdir+flist[0]) as file:
            for n, fline in enumerate(file):
                val=fline.split()
                if len(val)==0:
                    break
                else:
                    E.append(val[0])
        E=np.array(E).astype(np.float)
        tot_E=len(E)
        
        # read species from input file
        with open(self.wkdir+'elk.in','r') as file:
            spec=[line.replace('.',' ').replace(':',' ').split()[0] \
            for line in file if line.find(':spfname') is not -1]
        # check spin
        with open(self.wkdir+flist[0]) as file:
            flines=file.readlines()
            if len(flines)==(tot_E+1)*16:
                spin=1
            elif len(flines)==(tot_E+1)*32:
                spin=2
            else:
                print('Error: cannot determine spin in pdos!')
                print(len(flines),(tot_E+1)*16,(tot_E+1)*32)
                sys.exit()
        
        # read pdos
        state_info=[]
        state_count=0
        for atn, fname in enumerate(flist):
            site=[spec[int(fname[6:8])-1],int(fname[10:14])]
            for ms in range(1,spin+1):
                for l in range(0,4):
                    for ml in range(-l,+l+1):                    
                        state_info.append(\
                        '{0:>3d} => {1:>2s}{2:<} ({3:>3d} / {4} / {5:+} / {6:+})\n'.\
                        format(state_count,site[0],site[1],atn+1,l,ml,int(np.sign(1.5-ms))))
                        state_count+=1
                    
            with open(self.wkdir+fname) as file:
                flines=file.readlines()
            for lm in range(0,16*spin):
                tmp=np.array([val.split() for val in \
                flines[(tot_E+1)*lm:(tot_E+1)*(lm+1)-1]]).astype(np.float)
                
                if (tmp[0,0]!=E[0]) or (tmp[-1,0]!=E[-1]):    
                    print('pdos read mismatch!')
                    print(tmp[0,0],tmp[-1,0])
                    print(E[0],E[-1])
                    sys.exit()
                else:
                    tmp=tmp[:,1]
                    tmp.shape=(tot_E,1)
                    if 'pdos' in locals():
                        pdos=np.concatenate((pdos,tmp),1)
                    else:
                        pdos=tmp
        # check pdos size
        if pdos.shape[1]!=len(state_info):
            print('Error: size of pdos and state_info mismatch!')
            sys.exit()
            
        pdos=np.abs(pdos)
        E*=27.2    
        print(np.array(state_info))
        with open(self.wkdir+'pdos-state.dat','w') as file:
            file.writelines(state_info)
        np.savez(self.wkdir+'pdos.npz',E=E,pdos=pdos,state_info=state_info)
        
    def pdos_plot(self,state_grp,Ebound,lw=3,fontsize=18):
        print('starting pdos_plot ...')
        pdos=np.load(self.wkdir+'pdos.npz')
        dftpp.spectral_plot(self,pdos['E'],pdos['pdos'],pdos['state_info'],\
        self.state_grp_trans(pdos['state_info'],state_grp),savefig_path=self.wkdir+'pdos.png',Ebound=Ebound,\
        lw=lw,fontsize=fontsize)
        
