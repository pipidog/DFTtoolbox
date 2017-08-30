from DFTtoolbox.qe import postproc
import os
# Parameter ========================================================
run_task=[1,2,3,4]
wkdir=os.path.dirname(os.path.realpath(__file__))
# band_read & fatband_read
Ef=13.6197
#band_plot
kdiv=[15,7,5,15,13,9,5,10,9,1]
klabel=['$\Gamma$','X','W','K','$\Gamma$','L','U','W','L','K']
Ebound=[-5,5]
#fatband_plot
state_grp=[['1:1/2/1.5/a'],['2:2/1/1.5/a']]
# Main ================================================================
pp=postproc(wkdir)
for task in run_task:
    if task==1: #'band_read':
        pp.band_read(Ef=Ef,bandfile='pw.bands.out')
    elif task==2: #'band_plot':
        pp.band_plot(kdiv=kdiv,klabel=klabel,Ebound=Ebound)
    elif task==3: #'fatband_read':
        pp.fatband_read(Ef=Ef,projout='projwfc.fat.out',projprefix='fatband')
    elif task==4: #'fatband_plot':
        pp.fatband_plot(state_grp=state_grp,kdiv=kdiv,klabel=klabel,Ebound=Ebound)
    elif task==5: #pdos_read:
    	pp.pdos_read(Ef=Ef)
    elif task==6:
        pp.pdos_plot(state_grp=state_grp,Ebound=Ebound)