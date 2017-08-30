There are four subfolders in this folder:
* /init/ : contains an example tells you how to generate an input file.
 simply edit the file abinit_init.py. Then run the file using Python. 
 
 Let's set task='ground, it will ask DFTtoolbox to generate all necessary
 files for four basic ground state calculations: scf, band, fatband, pdos. 
 This input file is not ready to use. There are several variables you need
 to edit. For example, if you'r doing a magnetic calculation, you will 
 still need to input the magnetic moments in "bfcmt". Also, if you include
 DFT+U, you will still need to tell elk what are your correlated orbitals
 , U, J, etc. So tweak all the variables carefully to build own calculation.

 Similarly, by setting task='relax', you will get an input for structral 
 optimization. Typically, structral optimization doesn't consider SOC and 
 magnetism, so there are not much to tweak. Again, check the whole input
 file carefully to make sure it is what you want. 
 
 
 Another file abinit.file is also generate. This is the file that abinit will
 need. You will need to put path of your PAW pseudopotential.

 To make elk run, simple type:
 $> elk > elk.log

* /lda/, /lsda/:
 In these folder, there is an example of FeO calculation. All the output that
 are needed for DFTtoolbox are also inclued. Simplly edit the file elk_pp.py, 
 and run the file using Python. DFTtoolbox will automatically read the output
 files, save them in .npz format, and plot it in png format. 
 
 Band character analysis is an important feature of DFTtoolbox. Everytime when
 fatband_read() or pdos_read() were performed, they will generate a file call
 fatband-state.dat and PDOS-state.dat. In these file, it list all the 
 states in the follwing format, e.g. in the lsda folder:
 
  2 =>  O1 (  1 / 2 / a / a)
 
 The first number is the label of the state.O1 means the first O in "atoms" section
 in elk.in. Then comes with four numbers (1 / 2 / a / a). The first one means it
 is the first atom of the system. The other three numbers are quantum numbers,
 i.e., in LDSA case, (L, Ml, Ms). Similar idea apply for lda or SOC case. 
 Note that, because elk only supports fatband resoultion up to L, so Ml and Ms
 are 'a' means 'all'. If you want higher resoultion up to (Ml, Ms), you have
 to calculate PDOS (see PDOS-state.dat). It's nothing related to DFTtoolbox but
 a limition of elk.  
 
 * How to set state_gpr variable?
 There is a variable "state_grp" in all abinit_pp.py. This variable is used as
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