from DFTtoolbox.abinit import postproc
import os

# Parameters =======================================
run_task=[3,4,5,6]
wkdir=os.path.dirname(os.path.realpath(__file__))
dataset_ban=2
dataset_fatban=3
dataset_pdos=4
# band_read
Ef=0.33929
# band_plot & fatband_plot
kdiv=[15, 7, 5, 15, 13, 9, 5, 10, 9]
klabel=['$\Gamma$','X','W','K','$\Gamma$','L','U','W','L','K']
Ebound=[-5,5]
state_grp=[['1:1/2/a/a'],['2:2/1/a/a']]
# Main =============================================
print(wkdir)
pp=postproc(wkdir)
for task in run_task:
    if task is 1: #'band_read':
        pp.band_read(Ef=Ef,dataset=dataset_ban)
    elif task is 2 : #'band_plot':
        pp.band_plot(dataset=dataset_ban,kdiv=kdiv,klabel=klabel,Ebound=Ebound)
    elif task is 3: #'fatband_read':
        pp.fatband_read(dataset=dataset_fatban)
    elif task is 4: #'fatband_plot':
        pp.fatband_plot(dataset=dataset_fatban,kdiv=kdiv,\
        klabel=klabel,Ebound=Ebound,state_grp=state_grp)
    elif task is 5: # pdos_read
    	pp.pdos_read(dataset=dataset_pdos)
    elif task is 6: # pdos_plot
        pp.pdos_plot(dataset=dataset_pdos,state_grp=state_grp,Ebound=Ebound)

