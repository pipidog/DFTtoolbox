There are four subfolders in this folder:
* /init/ : contains an example tells you how to generate an input file.
 simply edit the file qe_init.py. Then run the file using Python. 
 
 Let's set task='ground, it will ask DFTtoolbox to generate all necessary
 files for four basic ground state calculations: scf, band, fatband, pdos. 
 In QE, these calcuations are achieved by several files:
 
 pw.scf.in (scf calcuations), projwfc.pdos.in (PDOS calcuations)
 pw.bands.in (band calcuations), projwfc.fat.in (fatband calcuations)
 
 These input file is not ready to use. There are several variables you need
 to edit. For example, if you'r doing a magnetic calculation, you will 
 still need to input the magnetic moments in "starting_magnetization". 
 Also, if you include DFT+U, you will still need to tell QE what are your 
 U, J, etc. (Note: QE already assume U will only apply to the outer d or 
 f orbital, so you don't need to assign it). So check all the variables 
 carefully to tweak own calculation.

 Similarly, by setting task='relax', you will get an input for structral 
 optimization. Typically, structral optimization doesn't consider SOC and 
 magnetism, so there are not much to tweak. Again, check the whole input
 file carefully to make sure it is what you want. 
 
 Note that, QE allows users to arbitarily name their input and output files.
 It makes DFTtoolbox difficult to find the correct files. Therefore, I 
 would strongly recommend to keep all file name default by DFTtoolbox. 
 As a result, you should run your calcuations in the following way:
 
    pw.x < pw.scf.in > pw.scf.out
    projwfc.x < projwfc.dos.in > projwfc.dos.out  (must come after scf, if you want PDOS)
    pw.x < pw.bands.in > pw.bands.out
    projwfc.x < projwfc.fat.in > projwfc.fat.out
    
If you really want to use different names, you can still do it. However, you 
will need to feed DFTtoolbox the names. Not a big deal, just more inputs. 


* /lda/, /lsda/, /soc/:
 In these folder, there is an example of FeO calculation. All the output that
 are needed for DFTtoolbox are also inclued. Simplly edit the file qe_pp.py, 
 and run the file using Python. DFTtoolbox will automatically read the output
 files, save them in .npz format, and plot it in png format. 
 
 Band character analysis is an important feature of DFTtoolbox. Everytime when
 fatband_read() or pdos_read() were performed, they will generate a file call
 fatband-DSx-state.dat and PDOS-DSx-state.dat. In these file, it list all the 
 states in the follwing format, e.g. in the lsda folder:
 
  6 => Fe-1  orb-4  ( 1 / 1 / 2 / +1)
 
 The first number is the label of the state. Fe-1 means the "1st" atom is "Fe". 
 orb-4 means it is the "4-th" orbital in Fe-1. In QE, the band projection is 
 determined by your pseudopotential. Sometimes if you include, say two s-orbitals,
 you will find there are two orbitals with L=0. However, you will still have
 different orbital number. After obtital number, it comes with four numbers 
 (1 / 1 / 2 / +1). The first one is the number of atom. So 1 means
 the 1-st atom in your input file. The other three numbers are quantum numbers,
 i.e., in LDSA case, (L, Ml, Ms). Similar idea apply for lda or SOC case. 
 
 Note that, for the meaning of Ml, you should check QE's official website:
 http://www.quantum-espresso.org/wp-content/uploads/Doc/INPUT_PROJWFC.html
 see "orbital order" section. 
 
 * How to set state_gpr variable?
 There is a variable "state_grp" in all abinit_pp.py. This variable is used ask
 a filter for fatband and PDOS plot. It uses the following format, e.g:
   '1:3/2/1/-1'
 It means to select site 1 to 3 with L=2, Ml=1 and Ms=-1 (only site part allows
 1:3 convention). If you want to select all L=2 at site 1:3, then use
   '1:3/2/a/a', 
 a means 'all'. You can also try, say,
   '1:3/a/a/1'
 it will select all spin-up at site 1:3.    
   
 Sometimes we want to combine several state in a single plot or a single lines, use
 the following convention:
  [['1:3/2/a/a','4:6/1/-1/+1'],['7:9/1/a/a']]
 Then in fatband, the weights of states '1:3/2/a/a' and '4:6/1/-1/+1' will be added
 and display in a single fatband plot. state '7:9/1/a/a' will be another fatband plot.
 Note that, even if you just need a sinlge fatband plot, you still have to use double
 list notation [['1:3/2/a/a']] because DFTtoolbox always considers the fisrt sublist
 as a single plot in fatband or a single curve in PDOS.  