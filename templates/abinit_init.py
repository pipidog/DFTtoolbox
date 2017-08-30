from DFTtoolbox.abinit import init
import os
# parameters ====================================
wkdir=os.path.dirname(os.path.realpath(__file__))
task='relax'   # 'relax' / 'ground'
prefix='FeO'   # prefix of xsf and kpf files
soc='on'       # 'on'/'off' of spin-orbit coupling (for gs only)
mag='on'       # 'on'/'off' magnetic calculation (for gs only)
dftu='on'      # 'on'/'off' DFT+U calculation    (for gs only)
kdense=20      # how many k-points per 1/A in band path? (default:20, for gs only)
# run ===========================================
c=init(wkdir)
if task=='ground':
	c.ground(prefix,soc,mag,dftu,kdense)
elif task=='relax':
	c.relax(prefix)