from DFTtoolbox.elk import postproc
import os
# Parameters =======================================
run_task=[1,2,3,4,5,6]
wkdir=os.path.dirname(os.path.realpath(__file__))
klabel=['$\Gamma$','X','W','K','$\Gamma$','L','U','W','L','K']
Ebound=[-5,5]
state_grp=[['1:1/1/a/a'],['2:2/2/a/a']]

# Main =============================================
print(wkdir)
pp=postproc(wkdir)
for task in run_task:
    if task is 1:
        pp.band_read()
    elif task is 2:
        pp.band_plot(klabel,Ebound)
    elif task is 3:
        pp.fatband_read()
    elif task is 4:
        pp.fatband_plot(state_grp,klabel,Ebound)
    elif task is 5:
        pp.pdos_read()
    elif task is 6:
        pp.pdos_plot(state_grp,Ebound)