import numpy as np
import matplotlib.pyplot as plt

plt.subplots_adjust(hspace=0.4)

filenames = ['vd0','vd01']

tt_data = [np.loadtxt('cl_tt_'+filename+'.d') for filename in filenames]
ee_data = [np.loadtxt('cl_ee_'+filename+'.d') for filename in filenames]
bb_data = [np.loadtxt('cl_bb_'+filename+'.d') for filename in filenames]
te_data = [np.loadtxt('cl_te_'+filename+'.d') for filename in filenames]

colors = ['black', 'red', 'green', 'blue']

plt.subplot(221)
for i, d in enumerate(tt_data):
	plt.semilogx(d[:,0], d[:,1], color=colors[i])
	plt.semilogx(d[:,0], d[:,2], color=colors[i], linestyle='--')
	plt.semilogx(d[:,0], d[:,3], color=colors[i], linestyle=':')
	plt.title('TT')
	plt.xlabel(r'$\ell$')
	plt.ylabel(r'$D_{\ell}$')
	plt.grid(True)

plt.subplot(222)
for i, d in enumerate(ee_data):
	plt.loglog(d[:,0], d[:,1], color=colors[i])
	plt.loglog(d[:,0], d[:,2], color=colors[i], linestyle='--')
	plt.loglog(d[:,0], d[:,3], color=colors[i], linestyle=':')
	plt.title('EE')
	plt.xlabel(r'$\ell$')
	plt.ylabel(r'$D_{\ell}$')
	plt.grid(True)

plt.subplot(223)
for i, d in enumerate(ee_data):
	plt.loglog(d[:,0], d[:,1], color=colors[i], linestyle='--')
	plt.loglog(d[:,0], d[:,2], color=colors[i], linestyle=':')
	plt.title('BB')
	plt.xlabel(r'$\ell$')
	plt.ylabel(r'$D_{\ell}$')
	plt.grid(True)

plt.show()