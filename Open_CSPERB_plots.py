import matplotlib.pyplot as plt
import numpy as N
import pickle
from sys import path

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib import colors
from scipy.interpolate import interp2d
plt.rcParams['mathtext.default'] = 'rm'
plt.rcParams['mathtext.fontset'] = 'stix'

def flux_limits(mf, Ts, flux_limits_file):
	# Net flux limits in W.m2 for Ts in K and for mass flows of 1 to 5 kg/s in OD 60.3 mm 740H pipes with 1.2 mm wall thickness.
	data = N.loadtxt(flux_limits_file, delimiter=',')
	mfs = data[0,1:]
	Tdata = data[1:,0]
	fluxlimdata = data[1:,1:]
	fit = interp2d(mfs, Tdata, fluxlimdata)
	flux_lim = fit(mf, Ts)
	return flux_lim[:,0]
	
def flux_limits_V(V, Ts, flux_limits_file):
	# Net flux limits in W.m2 for Ts in K and for mass flows of 1 to 5 kg/s in OD 60.3 mm 740H pipes with 1.2 mm wall thickness.
	data = N.loadtxt(flux_limits_file, delimiter=',')
	Vs = data[0,1:]
	Tdata = data[1:,0]
	fluxlimdata = data[1:,1:]
	fit = interp2d(Vs, Tdata, fluxlimdata)
	flux_lim = fit(V, Ts)
	return flux_lim[:,0]
	
def tower_receiver_plots(files, efficiency=True, maps_3D=True, flux_map=True, flow_paths=True, saveloc=None, billboard=False, flux_limits_file=None):
	fileo = open(files,'r')
	data = pickle.load(fileo)
	fileo.close()
	if saveloc == None:
		saveloc = files

	height = data['height']
	n_banks = data['n_banks']
	n_elems = data['n_elems']
	n_tubes = data['n_tubes']
	D_tubes_o = data['D_tubes_o']
	D_tubes_i = data['D_tubes_i'] 
	eff_abs = data['eff_abs'] 
	abs_t = data['abs_t']
	eff_ems = data['eff_ems']
	ems_t = data['ems_t']
	k_t = data['k_t']
	if billboard:
		width = data['width']
		ahr = data['wh']
		ahr_map = data['wh_map']
	else:
		radius = data['radius']
		ahr = data['ahr']
		ahr_map = data['ahr_map']
	fp = data['fp']
	areas = data['areas']
	areas_fp = data['areas_fp']
	HC = data['HC']
	T_in = data['T_in']
	T_out = data['T_out']
	h_conv_ext = data['h_conv_ext']
	h = data['h'] 
	m = data['m']
	flux_in = data['flux_in']
	q_net = data['q_net']
	q_rad = data['q_rad']
	q_ref = data['q_ref']
	q_conv = data['q_conv_ext']
	T_amb = data['T_amb']
	T_HC = data['T_HC']
	T_w_int = data['T_w_int']
	T_ext = data['T_ext']
	h_conv_int = data['h_conv_int']
	V = data['V']
	fluxmap = data['fluxmap']
	HC = data['HC']
	n_tubes = data['n_tubes']
	Dp = data['Dp'] 
	pipe_lengths = data['pipe_lengths']
	#maxheight = N.amax(ahr[:,1,:])
	#maxrad = N.amax(ahr[:,2,:])
	
	plt.rc('font', size=8)
	eff_rec=N.sum(q_net)/N.sum(fluxmap*areas[ahr_map])
	print ('Receiver efficiency: '+str(eff_rec))
	# Qin, eff_abs,eff_ems,T_ext_mean,h_ext,q_refl,q_emi,q_conv,eff_rec
	#print T_ext
	#print N.average(T_ext),N.sqrt(N.sqrt(N.average(T_ext**4))),h_conv_ext
	results=[N.sum(fluxmap*areas[ahr_map])/1e6,eff_abs,eff_ems,N.average(T_ext),N.sqrt(N.sqrt(N.average(T_ext**4))),h_conv_ext,N.sum(q_ref)/1e6,N.sum(q_rad)/1e6,N.sum(q_conv)/1e6,eff_rec]
	
	if efficiency:
		#print 'Radius:', radius,' m'
		eta_abs=(N.sum(fluxmap*areas[ahr_map])-N.sum(q_ref))/N.sum(fluxmap*areas[ahr_map])
		eta_th=N.sum(q_net)/(N.sum(fluxmap*areas[ahr_map])-N.sum(q_ref))
		eta_th_abs=N.sum(q_net)/N.sum(fluxmap*areas[ahr_map])
		q_in=N.sum(fluxmap*areas[ahr_map])/1e6
		q_refl=N.sum(q_ref)/1e6
		q_emi=N.sum(q_rad)/1e6
		q_conv=N.sum(q_conv)/1e6
		q_net_1=N.sum(q_net)/1e6
		Q_peak=N.amax(fluxmap)*eff_abs/1e3
		wall_t_max=N.amax(T_ext)
		for f in range(len(fp)):
			mfp = m[f]*N.ones(len(n_tubes[fp[f]]))
			nts = N.array(n_tubes[fp[f]])
			flux_lims = flux_limits(mfp/nts, T_HC[f], flux_limits_file)/1e3
			flux_fp = q_net[fp[f]]/areas[fp[f]]/1e3
			flux_fp = N.hstack((flux_fp[0], flux_fp))
			

		p_w = []
		vel=[]
		pressure_d=[]
		for f in xrange(len(fp)):
			vels = m[f]/n_tubes[fp[f]]/(N.pi*(D_tubes_i/2.)**2.*HC.rho(T_HC[f][1:]))
			p_w.append(m[f]*N.sum(Dp[fp[f]]/(HC.rho(T_HC[f][1:])+HC.rho(T_HC[f][:-1]))/2.))
			pressure_d.append(N.add.accumulate(Dp[fp[f]])[-1]/1e5)
			#print 'Heat transfer coefficient (kW/m2/K): ', h_conv_int[fp[f]]/1e3
			vel.append(N.amin(vels))
			vel.append(N.amax(vels))
		#print vel
		Fraction=[]		
		for f in range(len(fp)):
			mfp = m[f]*N.ones(len(n_tubes[fp[f]]))
			nts = N.array(n_tubes[fp[f]])
			flux_lims = flux_limits(mfp/nts, T_HC[f], flux_limits_file)/1e3
			flux_fp = q_net[fp[f]]/areas[fp[f]]/1e3
			flux_fp = N.hstack((flux_fp[0], flux_fp))
			Fraction.append(N.around(N.amax(N.hstack(flux_fp/flux_lims)), decimals=3))

	if maps_3D:

		fig = plt.figure(figsize=(6,3.4), dpi=1000)
		plt.subplots_adjust(left=0., bottom=0.06, top =1.02, right=1.01, wspace=0)
	
		angres_per_bank = 1
		ax = fig.add_subplot(1,2,1, projection='3d', aspect=1.3)

		flux = q_net/areas/1e6 #N.hstack(fluxmap)[N.hstack(ahr_map)]/1e6 
		# Set-up the colormap:
		cNorm  = colors.Normalize(vmin=0., vmax=1.2)
		mp = cm.ScalarMappable(norm=cNorm, cmap=cm.plasma)

		for i in xrange(len(flux)):
			ts = ahr[i,0,:]
			hs = ahr[i,1,:]
			rs = ahr[i,2,0]
			dt = (ts[1]-ts[0])/float(angres_per_bank)
			for j in xrange(angres_per_bank):
				thetas = [ts[0]+float(j)*dt, ts[0]+(j+1.)*dt]

				x = rs*N.cos(thetas)

				y = rs*N.sin(thetas)

				z0 = N.ones(2)*hs[0]
				z1 = N.ones(2)*hs[1]

				x = N.vstack((x, x))
				y = N.vstack((y, y))
				z = N.vstack((z0, z1))

				face = mp.to_rgba(flux[i])
				surf = ax.plot_surface(x, y, z, shade=False, antialiased=False, color=face, edgecolor=face, linewidth=0.5, alpha=1)

		rcard = N.amax(radius)+0.5
		zcard = N.amax(height)+0.5

		ax.scatter(xs=[0,rcard,0,-rcard], ys=[rcard,0,-rcard,0], zs=[zcard,zcard,zcard,zcard], s=150, c=N.tile([1,1,1,0.5],(4,1)), lw=0, zorder=9999.)
		ax.text(-rcard,0,zcard, s='N', zorder=10000., ha='center', va='center')
		ax.text(0,rcard,zcard, s='E', zorder=10000., ha='center', va='center')
		ax.text(rcard,0,zcard, s='S', zorder=10000., ha='center', va='center')
		ax.text(0,-rcard,zcard, s='W', zorder=10000., ha='center', va='center')

		ax.set_xlim(-radius, radius)
		ax.set_ylim(-radius, radius)
		ax.set_zlim(0, height)

		ax.view_init(elev=25, azim=95)

		mp = cm.ScalarMappable(cmap=cm.plasma)
		mp.set_array(N.hstack([N.amin(flux),N.amax(flux)]))
		cbar = plt.colorbar(mp, orientation='horizontal', pad=0.1, shrink=0.8)
		cbar.set_label(r'$\dot{q}^{\prime \prime}_\mathrm{abs}$ [MW.m$^{-2}$]')

		plt.figtext(0.25, 0.01, '(a)')
		plt.savefig(open(saveloc+'_3D_maps.png','w'), dpi=400)
		plt.clf()
		ax = fig.add_subplot(1,2,2, projection='3d', aspect=1.3)

		flux = T_ext[N.hstack(ahr_map)]
		# Set-up the colormap:
		cNorm  = colors.Normalize(vmin=N.amin(flux), vmax=N.amax(flux))
		mp = cm.ScalarMappable(norm=cNorm, cmap=cm.hot)

		for i in xrange(len(flux)):
			ts = ahr[i,0,:]
			hs = ahr[i,1,:]
			rs = ahr[i,2,0]
			for j in xrange(angres_per_bank):
				thetas = [ts[0]+float(j)/angres_per_bank*(ts[1]-ts[0]), ts[0]+(j+1.)/angres_per_bank*(ts[1]-ts[0])]

				x = rs*N.cos(thetas)

				y = rs*N.sin(thetas)

				z0 = N.ones(2)*hs[0]
				z1 = N.ones(2)*hs[1]

				x = N.vstack((x, x))
				y = N.vstack((y, y))
				z = N.vstack((z0, z1))

				face = mp.to_rgba(T_ext[i])
				surf = ax.plot_surface(x, y, z, shade=False, antialiased=False, color=face, edgecolor=face, linewidth=0.5,  alpha=1)

		ax.scatter(xs=[0,rcard,0,-rcard], ys=[rcard,0,-rcard,0], zs=[zcard,zcard,zcard,zcard], s=150, c=N.tile([1,1,1,0.5],(4,1)), lw=0, zorder=9999.)
		ax.text(-rcard,0,zcard, s='N', zorder=10000., ha='center', va='center')
		ax.text(0,rcard,zcard, s='E', zorder=10000., ha='center', va='center')
		ax.text(rcard,0,zcard, s='S', zorder=10000., ha='center', va='center')
		ax.text(0,-rcard,zcard, s='W', zorder=10000., ha='center', va='center')

		ax.set_xlim(-radius, radius)
		ax.set_ylim(-radius, radius)
		ax.set_zlim(0, height)

		ax.view_init(elev=25, azim=95)

		mp = cm.ScalarMappable(cmap=cm.hot)
		mp.set_array(N.hstack([N.amin(T_ext-273.15),N.amax(T_ext-273.15)]))
		cbar = plt.colorbar(mp, orientation='horizontal', pad=0.1, shrink=0.8)
		cbar.set_label(r'$T_\mathrm{abs}$ [$^{\circ}$C]')

		plt.figtext(0.75, 0.01, '(b)')

		plt.savefig(open(saveloc+'_3D_maps.png','w'), dpi=400)
		plt.clf()
		plt.close(fig)

	if flux_map:
		fig = plt.figure(figsize=(3,3.4), dpi=1000)
		plt.subplots_adjust(left=0., bottom=0.06, top =1.02, right=1.01, wspace=0)
	
		angres_per_bank = 1
		ax = fig.add_subplot(1,1,1, projection='3d', aspect=1.3)
		banks = N.shape(fluxmap)[1]
		fluxmap = N.concatenate((fluxmap[:,banks/4:], fluxmap[:,:banks/4]), axis=1) # correction to start from east
		flux = N.hstack(fluxmap)/1e6
		
		# Set-up the colormap:
		if isinstance(flux_map,list) and len(flux_map)==2:
			vmin, vmax = flux_map
		else:
			vmin=0.
			vmax=1.2

		cNorm  = colors.Normalize(vmin, vmax)
		mp = cm.ScalarMappable(norm=cNorm, cmap=cm.viridis)

		for i in xrange(len(flux)):
			ts = ahr[i,0,:]
			hs = ahr[i,1,:]
			rs = ahr[i,2,0]
			dt = (ts[1]-ts[0])/float(angres_per_bank)
			for j in xrange(angres_per_bank):
				thetas = [ts[0]+float(j)*dt, ts[0]+(j+1.)*dt]

				x = rs*N.cos(thetas)

				y = rs*N.sin(thetas)

				z0 = N.ones(2)*hs[0]
				z1 = N.ones(2)*hs[1]

				x = N.vstack((x, x))
				y = N.vstack((y, y))
				z = N.vstack((z0, z1))

				face = mp.to_rgba(flux[i])
				surf = ax.plot_surface(x, y, z, shade=False, antialiased=False, color=face, edgecolor=face, linewidth=0.5, alpha=1)

		rcard = N.amax(radius)+0.5
		zcard = N.amax(height)+0.5

		ax.scatter(xs=[0,rcard,0,-rcard], ys=[rcard,0,-rcard,0], zs=[zcard,zcard,zcard,zcard], s=150, c=N.tile([1,1,1,0.5],(4,1)), lw=0, zorder=9999.)
		ax.text(-rcard,0,zcard, s='N', zorder=10000., ha='center', va='center')
		ax.text(0,rcard,zcard, s='E', zorder=10000., ha='center', va='center')
		ax.text(rcard,0,zcard, s='S', zorder=10000., ha='center', va='center')
		ax.text(0,-rcard,zcard, s='W', zorder=10000., ha='center', va='center')

		ax.set_xlim(-radius, radius)
		ax.set_ylim(-radius, radius)
		ax.set_zlim(0, height)

		ax.view_init(elev=25, azim=95)

		mp = cm.ScalarMappable(cmap=cm.viridis)
		mp.set_array(N.hstack([vmin, vmax]))
		cbar = plt.colorbar(mp, orientation='horizontal', pad=0.1, shrink=0.8)
		cbar.set_label(r'$\dot{q}^{\prime \prime}_\mathrm{in}$ [MW.m$^{-2}$]')

		plt.savefig(open(saveloc+'_3D_flux.png','w'), dpi=400)
		plt.clf()

		fig = plt.figure(figsize=(3,3.4), dpi=1000)
		plt.subplots_adjust(left=0., bottom=0.10, top =1.02, right=1.01, wspace=0)
	
		angres_per_bank = 1
		ax = fig.add_subplot(1,1,1, projection='3d', aspect=1.3)
		banks = N.shape(fluxmap)[1]
		fluxmap = N.concatenate((fluxmap[:,banks/4:], fluxmap[:,:banks/4]), axis=1) # correction to start from east
		flux = T_ext[N.hstack(ahr_map)]-273.15
		
		# Set-up the colormap:
		if isinstance(flux_map,list) and len(flux_map)==2:
			vmin, vmax = flux_map
		else:
			vmin=520.
			vmax=780.

		cNorm  = colors.Normalize(vmin, vmax)
		mp = cm.ScalarMappable(norm=cNorm, cmap=cm.hot)

		for i in xrange(len(flux)):
			ts = ahr[i,0,:]
			hs = ahr[i,1,:]
			rs = ahr[i,2,0]
			dt = (ts[1]-ts[0])/float(angres_per_bank)
			for j in xrange(angres_per_bank):
				thetas = [ts[0]+float(j)*dt, ts[0]+(j+1.)*dt]

				x = rs*N.cos(thetas)

				y = rs*N.sin(thetas)

				z0 = N.ones(2)*hs[0]
				z1 = N.ones(2)*hs[1]

				x = N.vstack((x, x))
				y = N.vstack((y, y))
				z = N.vstack((z0, z1))

				face = mp.to_rgba(flux[i])
				surf = ax.plot_surface(x, y, z, shade=False, antialiased=False, color=face, edgecolor=face, linewidth=0.5, alpha=1)

		rcard = N.amax(radius)+0.5
		zcard = N.amax(height)+0.5

		ax.scatter(xs=[0,rcard,0,-rcard], ys=[rcard,0,-rcard,0], zs=[zcard,zcard,zcard,zcard], s=150, c=N.tile([1,1,1,0.5],(4,1)), lw=0, zorder=9999.)
		ax.text(-rcard,0,zcard, s='N', zorder=10000., ha='center', va='center')
		ax.text(0,rcard,zcard, s='E', zorder=10000., ha='center', va='center')
		ax.text(rcard,0,zcard, s='S', zorder=10000., ha='center', va='center')
		ax.text(0,-rcard,zcard, s='W', zorder=10000., ha='center', va='center')

		ax.set_xlim(-radius, radius)
		ax.set_ylim(-radius, radius)
		ax.set_zlim(0, height)

		ax.view_init(elev=25, azim=95)

		mp = cm.ScalarMappable(cmap=cm.hot)
		mp.set_array(N.hstack([vmin, vmax]))
		cbar = plt.colorbar(mp, orientation='horizontal', pad=0.1, shrink=0.8)
		cbar.ax.tick_params(labelsize=16) 
		cbar.set_label(r'${\it{T}_\mathrm{w}}$ ($^{\circ}$C)',size=16)

		plt.savefig(open(saveloc+'_3D_T.png','w'), dpi=400)
		plt.clf()
		
		fig = plt.figure(figsize=(3,3.4), dpi=1000)
		plt.subplots_adjust(left=0., bottom=0.10, top =1.02, right=1.01, wspace=0)
	
		angres_per_bank = 1
		ax = fig.add_subplot(1,1,1, projection='3d', aspect=1.3)
		banks = N.shape(fluxmap)[1]
		fluxmap = N.concatenate((fluxmap[:,banks/4:], fluxmap[:,:banks/4]), axis=1) # correction to start from east
		flux = q_net/areas/1e6
		
		# Set-up the colormap:
		if isinstance(flux_map,list) and len(flux_map)==2:
			vmin, vmax = flux_map
		else:
			vmin=0.
			vmax=1.2

		cNorm  = colors.Normalize(vmin, vmax)
		mp = cm.ScalarMappable(norm=cNorm, cmap=cm.plasma)

		for i in xrange(len(flux)):
			ts = ahr[i,0,:]
			hs = ahr[i,1,:]
			rs = ahr[i,2,0]
			dt = (ts[1]-ts[0])/float(angres_per_bank)
			for j in xrange(angres_per_bank):
				thetas = [ts[0]+float(j)*dt, ts[0]+(j+1.)*dt]

				x = rs*N.cos(thetas)

				y = rs*N.sin(thetas)

				z0 = N.ones(2)*hs[0]
				z1 = N.ones(2)*hs[1]

				x = N.vstack((x, x))
				y = N.vstack((y, y))
				z = N.vstack((z0, z1))

				face = mp.to_rgba(flux[i])
				surf = ax.plot_surface(x, y, z, shade=False, antialiased=False, color=face, edgecolor=face, linewidth=0.5, alpha=1)

		rcard = N.amax(radius)+0.5
		zcard = N.amax(height)+0.5

		ax.scatter(xs=[0,rcard,0,-rcard], ys=[rcard,0,-rcard,0], zs=[zcard,zcard,zcard,zcard], s=150, c=N.tile([1,1,1,0.5],(4,1)), lw=0, zorder=9999.)
		ax.text(-rcard,0,zcard, s='N', zorder=10000., ha='center', va='center')
		ax.text(0,rcard,zcard, s='E', zorder=10000., ha='center', va='center')
		ax.text(rcard,0,zcard, s='S', zorder=10000., ha='center', va='center')
		ax.text(0,-rcard,zcard, s='W', zorder=10000., ha='center', va='center')

		ax.set_xlim(-radius, radius)
		ax.set_ylim(-radius, radius)
		ax.set_zlim(0, height)

		ax.view_init(elev=25, azim=95)

		mp = cm.ScalarMappable(cmap=cm.plasma)
		mp.set_array(N.hstack([vmin,vmax]))
		cbar = plt.colorbar(mp, orientation='horizontal', pad=0.1, shrink=0.8)
		cbar.ax.tick_params(labelsize=16) 
		cbar.set_label(r'${\dot{\it{Q}}}_\mathrm{net}$ (MW.m$^{-2}$)',size=16)

		plt.savefig(open(saveloc+'_3D_Q-net.png','w'), dpi=400)
		plt.clf()

	if flow_paths:
		'''
		fig = plt.figure(figsize=(8*2,(2.+(len(fp)-1)*1.1)/1.7), dpi=1000)
		bot=0.08
		top=0.93
		
		plt.subplots_adjust(left=0.5, bottom=bot, right=0.95, top = top)
		for f in xrange(len(fp)):
			plt.subplot(int(len(fp)/4),4,f+1)
			if len(fp)>1:
				plt.text(x=1, y=660, s='Flow path %s'%str(f+1), va='top')
			bank_lengths = pipe_lengths[f]
			bank_lengths_2 = (bank_lengths[1:]+bank_lengths[:-1])/2.
			plt.plot(bank_lengths_2, T_ext[fp[f]]-273.15, label=r'T${_\mathrm{abs}}$', color='0')
			plt.plot(bank_lengths_2, T_w_int[fp[f]]-273.15, label=r'T${_\mathrm{wall,int}}$', color='0.4')
			plt.plot(bank_lengths, T_HC[f]-273.15, label=r'T${_\mathrm{HC}}$', color='0.6')
			plt.xlim(xmax=bank_lengths[-1])

			plt.ylabel('T [$^{\circ}$C]')
		plt.xlabel('Flow path length [m]')

		plt.subplot(len(fp),1,1)
		plt.legend(loc=3,ncol=3, borderaxespad=0, bbox_to_anchor=(0.,1.05))

		#plt.savefig(open(saveloc+'_Temp_fp.png','w'), dpi=400)
		#plt.close(fig)
		'''
		bot=0.08
		top=0.93
		if n_banks/len(fp)==2:
			fig = plt.figure(figsize=(4*2,(2.+(len(fp)-1)*1.1)/1.7), dpi=1000)
		else: 
			fig = plt.figure(figsize=(10,8), dpi=1000)
		plt.subplots_adjust(left=0.15, bottom=0.12, right=0.98, top = top)
		Success=[]
		Positive=[]
		safety_factor=0.9
		A_over=N.array([])
		for f in xrange(len(fp)):
			bank_lengths = pipe_lengths[f]
			bank_lengths_2 = (bank_lengths[1:]+bank_lengths[:-1])/2.
			if n_banks/len(fp)==2:
				ax=plt.subplot(int(len(fp)/2),2,f+1)
				size=12
			else:
				ax=plt.subplot(int(len(fp)/4),4,f+1)
				size=24
			if len(fp)>1:
				plt.text(x=bank_lengths[-5], y=1000, s='fp %s'%str(f+1), ha='right',fontsize=size)
			
			ax.plot(bank_lengths_2, q_net[fp[f]]/areas[fp[f]]/1e3, label=r'${\dot{q}^{\prime \prime}_\mathrm{abs}}$', color='0.6')
			flux_lims = flux_limits_V(V[f], T_HC[f], flux_limits_file)/1e3
			safe_flux_lims = safety_factor*flux_lims
			ax.plot(bank_lengths, flux_lims, color='r',linewidth=1.,label=r'${\dot{q}^{\prime \prime}_\mathrm{limit}}$')
			ax.plot(bank_lengths, safe_flux_lims, color='r',linestyle='--',linewidth=1.,label=r'${\dot{q}^{\prime \prime}_\mathrm{safe}}$')
			if n_banks/len(fp)==2:
				plt.vlines(x=height,ymin=0,ymax=1500,linestyles='--',linewidth=0.5)
			D=safety_factor*flux_lims[:-1]-q_net[fp[f]]/areas[fp[f]]/1e3 # difference between flux limits and net flux
			#print D
			# for front tubes
			Index=N.where(D[:n_elems]<0)
			if Index[0]!=[]:
				if abs(Index[0][0]-n_elems/2-1)>=abs(Index[0][-1]-n_elems/2-1):
					boundary_1=Index[0][0]
					boundary_2=2*(n_elems/2-1)-boundary_1
				else:
					boundary_2=Index[0][-1]
					boundary_1=2*(n_elems/2-1)-boundary_2
				D_choose=D[boundary_1:boundary_2]
				index_positive=N.where(D_choose>=0)
				index_negative=N.where(D_choose<0)
				Success.append(sum(D_choose[index_positive])+5.*sum(D_choose[index_negative])>0.)
				if index_negative!=[]:
					Positive.append(False)
					A_over=N.append(A_over,-sum(D_choose[index_negative]))
			else:
				Positive.append(True)
				Success.append(True)
				A_over=N.append(A_over,0.)			
			
			if n_banks/len(fp)==2:
				Index=N.where(D[n_elems:]<0)
				if Index[0]!=[]:
					if abs(Index[0][0]-(n_elems/2-1))>=abs(Index[0][-1]-(n_elems/2-1)):
						boundary_1=Index[0][0]
						boundary_2=2*(n_elems/2-1)-boundary_1
					else:
						boundary_2=Index[0][-1]
						boundary_1=2*(n_elems/2-1)-boundary_2
					D_choose=D[boundary_1+n_elems:boundary_2+n_elems]
					index_positive=N.where(D_choose>=0)
					index_negative=N.where(D_choose<0)
					Success.append(sum(D_choose[index_positive])+5.*sum(D_choose[index_negative])>0.)
					if index_negative!=[]:
						Positive.append(False)
						A_over=N.append(A_over,-sum(D_choose[index_negative]))				
				else:
					Positive.append(True)
					Success.append(True)
					A_over=N.append(A_over,0.)
			#plt.hlines(y=0., xmin=0., xmax=N.amax(bank_lengths_2), lw=0.5)
			'''
			if n_banks/len(fp)==2:
				if f%2==0:
					plt.ylabel('Flux (kW.m$^{-2}$)',fontsize=size)
				#if f==0:
					#plt.legend(loc=3,ncol=3, borderaxespad=0, bbox_to_anchor=(0.,1.05),fontsize=size)
				if f>5:
					plt.xlabel('Flow path length (m)',fontsize=size)
			else:
				if f%4==0:
					plt.ylabel('Flux (kW.m$^{-2}$)',fontsize=size)
				#if f==0:
					#plt.legend(loc=3,ncol=3, borderaxespad=0, bbox_to_anchor=(0.,1.05),fontsize=size)
				if f>11:
					plt.xlabel('Flow path length (m)',fontsize=size)
			plt.ylim(-100,1800)
			'''
			if f==13:
				plt.xlabel('                    $\it{L}_{fp}$ $(m)$',fontsize=size)
			if f==8:
				plt.ylabel('                    $\dot{\it{Q}}$ $(kW.m^{-2})$',fontsize=size)
			if f==0 or f==4 or f==8:
				ax.tick_params(axis='y', which='major', labelsize=size,direction='in')
				ax.tick_params(axis='x', which='major', bottom=False,top=False,labelbottom=False)
			elif f==12:
				ax.tick_params(axis='both', which='major', labelsize=size,direction='in')
			elif f==13 or f==14 or f==15:
				ax.tick_params(axis='x', which='major', labelsize=size,direction='in')
				ax.tick_params(axis='y', which='major', left=False,right=False,labelleft=False)
			else:
				ax.tick_params(axis='x', which='major', bottom=False,top=False,labelbottom=False)
				ax.tick_params(axis='y', which='major', left=False,right=False,labelleft=False)
			ax.locator_params(nbins=5)
		#plt.xlabel('Flow path length [m]')
		plt.savefig(open(saveloc+'_flux_fp.png','w'), dpi=400)
		plt.close(fig)
		plt.clf()
		aiming_results=[Success,Positive,A_over]
	return results,aiming_results

def flow_path_plot(files='/home/charles/Documents/Boulot/These/Sodium receiver_CMI/ref_case_result_1', fps=[0], saveloc=None, flux_limits_file=None):

	import matplotlib.gridspec as gridspec

	fileo = open(files,'r')
	data = pickle.load(fileo)
	fileo.close()
	if saveloc == None:
		saveloc = files

	ahr = data['ahr']
	radius = data['radius']
	height = data['height']
	n_banks = data['n_banks']
	n_elems = data['n_elems']
	n_tubes = data['n_tubes']
	D_tubes_o = data['D_tubes_o']
	D_tubes_i = data['D_tubes_i'] 
	eff_abs = data['eff_abs'] 
	abs_t = data['abs_t']
	eff_ems = data['eff_ems']
	ems_t = data['ems_t']
	k_t = data['k_t']
	ahr_map = data['ahr_map']
	fp = data['fp']
	areas = data['areas']
	areas_fp = data['areas_fp']
	HC = data['HC']
	T_in = data['T_in']
	T_out = data['T_out']
	h_conv_ext = data['h_conv_ext']
	h = data['h'] 
	m = data['m']
	flux_in = data['flux_in']
	q_net = data['q_net']
	q_rad = data['q_rad']
	q_ref = data['q_ref']
	q_conv = data['q_conv_ext']
	T_amb = data['T_amb']
	T_HC = data['T_HC']
	T_w_int = data['T_w_int']
	T_ext = data['T_ext']
	h_conv_int = data['h_conv_int']
	V = data['V']
	fluxmap = data['fluxmap']
	HC = data['HC']
	n_tubes = data['n_tubes']
	Dp = data['Dp'] 
	pipe_lengths = data['pipe_lengths']

	maxheight = N.amax(ahr[:,1,:])
	maxrad = N.amax(ahr[:,2,:])

	plt.rc('font', size=8)

	for flow_path in fps: 

		fig = plt.figure(figsize=(6.,6.))
		ncols = 2
		nrows = 5
		gs = gridspec.GridSpec(ncols=ncols, nrows=nrows, width_ratios=[.9,2.], left=0.01, top=.98, wspace=0.25, right=0.98, hspace=0.08, bottom=0.08)

		# Text:
		plt.figtext(x=0.02, y=0.95, s='Flow-path %s'%str(flow_path))
		plt.figtext(x=0.02, y=0.92, s='Net power gain: %s MW${_\mathrm{th}}$'%str(N.around(N.sum(q_net[fp[flow_path]]/1e6), decimals=1)))
		plt.figtext(x=0.02, y=0.89, s=r'Efficiency: %s%%'%str(N.around(N.sum(q_net[fp[flow_path]])/N.sum(flux_in[flow_path]*areas_fp[flow_path]), decimals=3)*100.))

		# 3D receiver snipet with only flow_path in dense:
		rec3Dsnipet = fig.add_subplot(gs[1:3,0], projection='3d', aspect=1.4)

		angres_per_bank = 1
		banks = N.shape(fluxmap)[1]
		fluxmap_fp = N.concatenate((fluxmap[:,banks/4:], fluxmap[:,:banks/4]), axis=1) # correction to start from east
		flux = N.hstack(fluxmap_fp)/1e6
		
		# Set-up the colormap:
		vmin=0.
		vmax=1.2

		cNorm  = colors.Normalize(vmin, vmax)
		mp = cm.ScalarMappable(norm=cNorm, cmap=cm.viridis)

		for i in xrange(len(flux)):
			if (i == N.hstack(fp[flow_path])).any():
				alpha = 1
				lw = 0.5
			else:
				alpha=0.1
				lw = 0
			ts = ahr[i,0,:]
			hs = ahr[i,1,:]
			rs = ahr[i,2,0]
			dt = (ts[1]-ts[0])/float(angres_per_bank)
			for j in xrange(angres_per_bank):
				thetas = [ts[0]+float(j)*dt, ts[0]+(j+1.)*dt]

				x = rs*N.cos(thetas)

				y = rs*N.sin(thetas)

				z0 = N.ones(2)*hs[0]
				z1 = N.ones(2)*hs[1]

				x = N.vstack((x, x))
				y = N.vstack((y, y))
				z = N.vstack((z0, z1))

				face = mp.to_rgba(flux[i])
				surf = rec3Dsnipet.plot_surface(x, y, z, shade=False, antialiased=False, color=face, edgecolor=face, linewidth=lw, alpha=alpha)

		rcard = N.amax(radius)+0.5
		zcard = N.amax(height)+0.5

		rec3Dsnipet.scatter(xs=[0,rcard,0,-rcard], ys=[rcard,0,-rcard,0], zs=[zcard,zcard,zcard,zcard], s=150, c=N.tile([1,1,1,0.7],(4,1)), lw=0, zorder=9999.)
		rec3Dsnipet.text(-rcard,0,zcard, s='N', zorder=10000., ha='center', va='center')
		rec3Dsnipet.text(0,rcard,zcard, s='E', zorder=10000., ha='center', va='center')
		rec3Dsnipet.text(rcard,0,zcard, s='S', zorder=10000., ha='center', va='center')
		rec3Dsnipet.text(0,-rcard,zcard, s='W', zorder=10000., ha='center', va='center')

		rec3Dsnipet.set_xlim(-radius/1.2, radius/1.2)
		rec3Dsnipet.set_ylim(-radius/1.2, radius/1.2)
		rec3Dsnipet.set_zlim(0, height)

		rec3Dsnipet.set_xticks([],[])
		rec3Dsnipet.set_yticks([],[])
		rec3Dsnipet.set_zticks([],[])

		# Get rid of colored axes planes
		# First remove fill
		rec3Dsnipet.xaxis.pane.fill = False
		rec3Dsnipet.yaxis.pane.fill = False
		rec3Dsnipet.zaxis.pane.fill = False

		rec3Dsnipet.set_axis_off()

		rec3Dsnipet.view_init(elev=25, azim=95)
	
		mp = cm.ScalarMappable(cmap=cm.viridis)
		mp.set_array(N.hstack([vmin, vmax]))
		cbar = plt.colorbar(mp, orientation='horizontal', pad=-0.05, shrink=0.8)
		cbar.set_label(r'$\dot{q}^{\prime \prime}_\mathrm{in}$ [MW.m$^{-2}$]')

		# T ext 3D:
		# 3D receiver snipet with only flow_path in dense:
		rec3Dsnipet = fig.add_subplot(gs[3:5,0], projection='3d', aspect=1.4)
		rec3Dsnipet.set_facecolor('None')
		angres_per_bank = 1
		banks = N.shape(fluxmap)[1]
		fluxmap_fp = N.concatenate((fluxmap[:,banks/4:], fluxmap[:,:banks/4]), axis=1) # correction to start from east
		flux = N.hstack(T_ext)-273.15
		
		# Set-up the colormap:
		vmin=N.amin(flux)
		vmax=N.amax(flux)

		cNorm  = colors.Normalize(vmin, vmax)
		mp = cm.ScalarMappable(norm=cNorm, cmap=cm.hot)

		for i in xrange(len(flux)):
			if (i == N.hstack(fp[flow_path])).any():
				alpha = 1
				lw = 0.5
			else:
				alpha=0.1
				lw = 0
			ts = ahr[i,0,:]
			hs = ahr[i,1,:]
			rs = ahr[i,2,0]
			dt = (ts[1]-ts[0])/float(angres_per_bank)
			for j in xrange(angres_per_bank):
				thetas = [ts[0]+float(j)*dt, ts[0]+(j+1.)*dt]

				x = rs*N.cos(thetas)

				y = rs*N.sin(thetas)

				z0 = N.ones(2)*hs[0]
				z1 = N.ones(2)*hs[1]

				x = N.vstack((x, x))
				y = N.vstack((y, y))
				z = N.vstack((z0, z1))

				face = mp.to_rgba(flux[i])
				surf = rec3Dsnipet.plot_surface(x, y, z, shade=False, antialiased=False, color=face, edgecolor=face, linewidth=lw, alpha=alpha)

		rcard = N.amax(radius)+0.5
		zcard = N.amax(height)+0.5

		rec3Dsnipet.scatter(xs=[0,rcard,0,-rcard], ys=[rcard,0,-rcard,0], zs=[zcard,zcard,zcard,zcard], s=150, c=N.tile([1,1,1,0.7],(4,1)), lw=0, zorder=9999.)
		rec3Dsnipet.text(-rcard,0,zcard, s='N', zorder=10000., ha='center', va='center')
		rec3Dsnipet.text(0,rcard,zcard, s='E', zorder=10000., ha='center', va='center')
		rec3Dsnipet.text(rcard,0,zcard, s='S', zorder=10000., ha='center', va='center')
		rec3Dsnipet.text(0,-rcard,zcard, s='W', zorder=10000., ha='center', va='center')

		rec3Dsnipet.set_xlim(-radius/1.2, radius/1.2)
		rec3Dsnipet.set_ylim(-radius/1.2, radius/1.2)
		rec3Dsnipet.set_zlim(0, height)

		rec3Dsnipet.set_xticks([],[])
		rec3Dsnipet.set_yticks([],[])
		rec3Dsnipet.set_zticks([],[])

		# Get rid of colored axes planes
		# First remove fill
		rec3Dsnipet.xaxis.pane.fill = False
		rec3Dsnipet.yaxis.pane.fill = False
		rec3Dsnipet.zaxis.pane.fill = False

		rec3Dsnipet.set_axis_off()

		rec3Dsnipet.view_init(elev=25, azim=95)
	
		mp = cm.ScalarMappable(cmap=cm.hot)
		mp.set_array(N.hstack([vmin, vmax]))
		cbar = plt.colorbar(mp, orientation='horizontal', pad=-0.05, shrink=0.8)
		cbar.set_label(r'${T_\mathrm{wall, ext}}$ [${^\circ}$C]')


		# flux fp plot:
		fluxplot = fig.add_subplot(gs[0,1])
		bank_lengths = pipe_lengths[flow_path]
		bank_lengths_2 = (bank_lengths[1:]+bank_lengths[:-1])/2.

		flux_lims = flux_limits(N.ones(len(fp[flow_path]))/n_tubes[fp[flow_path]]*m[flow_path], T_HC[flow_path], flux_limits_file)/1e3
		plt.plot(bank_lengths, flux_lims, color='r', label=r'${\dot{q}_\mathrm{LIMIT}^{\prime \prime}}$')
		plt.plot(bank_lengths_2, flux_in[flow_path]/1e3, label=r'${\dot{q}_\mathrm{in}^{\prime \prime}}$', color='0')
		plt.plot(bank_lengths_2, q_net[fp[flow_path]]/areas[fp[flow_path]]/1e3, label=r'${\dot{q}^{\prime \prime}_\mathrm{net}}$', color='0.6')


		plt.legend()

		xticks = plt.gca().get_xticks()
		plt.xticks(xticks, [])
		plt.tick_params(direction='in')
		plt.xlim(N.amin(bank_lengths), N.amax(bank_lengths))
		plt.ylabel('Flux [kW.m$^{-2}$]')

		# Temperatures:
		Tplot = fig.add_subplot(gs[1,1])
		plt.plot(bank_lengths_2, T_ext[fp[flow_path]]-273.15, label=r'T${_\mathrm{abs}}$', color='0')
		plt.plot(bank_lengths_2, T_w_int[fp[flow_path]]-273.15, label=r'T${_\mathrm{wall,int}}$', color='0.4')
		plt.plot(bank_lengths, T_HC[flow_path]-273.15, label=r'T${_\mathrm{HC}}$', color='0.6')

		plt.legend()

		xticks = plt.gca().get_xticks()
		plt.xticks(xticks, [])
		plt.tick_params(direction='in')
		plt.xlim(N.amin(bank_lengths), N.amax(bank_lengths))

		plt.ylabel('T [$^{\circ}$C]')

		# Fraction of allowable:
		Fracplot = fig.add_subplot(gs[2,1])
		
		plt.plot(bank_lengths_2, q_net[fp[flow_path]]/areas[fp[flow_path]]/1e3/((flux_lims[:-1]+flux_lims[1:])/2.), color='g')
		plt.ylabel(r'${\frac{\dot{q}_\mathrm{net}^{\prime \prime}}{\dot{q}_\mathrm{LIMIT}^{\prime \prime}}}$')
		
		xticks = plt.gca().get_xticks()
		plt.xticks(xticks, [])
		plt.tick_params(direction='in')
		plt.xlim(N.amin(bank_lengths), N.amax(bank_lengths))
		plt.ylim(-0.1,1.1)

		# Pressure drop:
		pdroplot = fig.add_subplot(gs[3,1])
		Dploc = -N.add.accumulate(Dp[fp[flow_path]])/1e5
		plt.plot(bank_lengths_2, Dploc, color='0')

		xticks = plt.gca().get_xticks()
		plt.xticks(xticks, [])
		plt.tick_params(direction='in')
		plt.xlim(N.amin(bank_lengths), N.amax(bank_lengths))

		plt.ylabel('${p(\ell)-p(0)}$ [bar]')

		# Velocity:
		vplot = fig.add_subplot(gs[4,1])
		plt.plot(bank_lengths_2, V[fp[flow_path]], color='0.5')

		plt.xlim(N.amin(bank_lengths), N.amax(bank_lengths))

		plt.xlabel('Flow path length (${\ell}$) [m]')
		plt.ylabel('${v_\mathrm{Na}}$ [m.s$^{-1}$]')


		plt.savefig(open(saveloc+'_fp_%s.png'%str(flow_path),'w'), dpi=400)
		plt.clf()
		plt.close(fig)
		
		

def flow_path_plot_billboard(files='/home/charles/Documents/Boulot/These/Sodium receiver_CMI/ref_case_result_1', fps=[0], saveloc=None, flux_limits_file=None):

	import matplotlib.gridspec as gridspec

	fileo = open(files,'r')
	data = pickle.load(fileo)
	fileo.close()
	if saveloc == None:
		saveloc = files

	wh = data['wh']
	width = data['width']
	height = data['height']
	n_banks = data['n_banks']
	n_elems = data['n_elems']
	n_tubes = data['n_tubes']
	D_tubes_o = data['D_tubes_o']
	D_tubes_i = data['D_tubes_i'] 
	eff_abs = data['eff_abs'] 
	abs_t = data['abs_t']
	eff_ems = data['eff_ems']
	ems_t = data['ems_t']
	k_t = data['k_t']
	wh_map = data['wh_map']
	fp = data['fp']
	areas = data['areas']
	areas_fp = data['areas_fp']
	HC = data['HC']
	T_in = data['T_in']
	T_out = data['T_out']
	h_conv_ext = data['h_conv_ext']
	h = data['h'] 
	m = data['m']
	flux_in = data['flux_in']
	q_net = data['q_net']
	q_rad = data['q_rad']
	q_ref = data['q_ref']
	q_conv = data['q_conv_ext']
	T_amb = data['T_amb']
	T_HC = data['T_HC']
	T_w_int = data['T_w_int']
	T_ext = data['T_ext']
	h_conv_int = data['h_conv_int']
	V = data['V']
	fluxmap = data['fluxmap']
	HC = data['HC']
	n_tubes = data['n_tubes']
	Dp = data['Dp'] 
	pipe_lengths = data['pipe_lengths']

	maxheight = N.amax(wh[:,1,:])
	maxwidth = N.amax(wh[:,0,:])

	plt.rc('font', size=8)

	for flow_path in fps: 

		fig = plt.figure(figsize=(6.,6.))
		ncols = 2
		nrows = 5
		gs = gridspec.GridSpec(ncols=ncols, nrows=nrows, width_ratios=[.9,2.], left=0.01, top=.98, wspace=0.25, right=0.98, hspace=0.08, bottom=0.08)

		# Text:
		plt.figtext(x=0.02, y=0.95, s='Flow-path %s'%str(flow_path))
		plt.figtext(x=0.02, y=0.92, s='Net power gain: %s MW${_\mathrm{th}}$'%str(N.around(N.sum(q_net[fp[flow_path]]/1e6), decimals=1)))
		plt.figtext(x=0.02, y=0.89, s=r'Efficiency: %s%%'%str(N.around(N.sum(q_net[fp[flow_path]])/N.sum(flux_in[flow_path]*areas_fp[flow_path])*100., decimals=1)))

		# receiver snipet with only flow_path in dense:
		rec3Dsnipet = fig.add_subplot(gs[1:3,0], projection='3d', aspect=1.4)

		angres_per_bank = 1
		banks = N.shape(fluxmap)[1]
		fluxmap_fp = fluxmap
		flux = N.hstack(fluxmap_fp)/1e6
		
		# Set-up the colormap:
		vmin=0.
		vmax=1.2

		cNorm  = colors.Normalize(vmin, vmax)
		mp = cm.ScalarMappable(norm=cNorm, cmap=cm.viridis)

		for i in xrange(len(flux)):
			if (i == N.hstack(fp[flow_path])).any():
				alpha = 1
				lw = 0.5
			else:
				alpha=0.1
				lw = 0
			ws = wh[i,0,:]
			hs = wh[i,1,:]
			for j in xrange(angres_per_bank):
				x = ws-width/2.
				y = ws*0

				z0 = N.ones(2)*hs[0]
				z1 = N.ones(2)*hs[1]

				x = N.vstack((x, x))
				y = N.vstack((y, y))
				z = N.vstack((z0, z1))

				face = mp.to_rgba(flux[i])
				surf = rec3Dsnipet.plot_surface(x, y, z, shade=False, antialiased=False, color=face, edgecolor=face, linewidth=lw, alpha=alpha)

		rec3Dsnipet.set_xlim(-width/2.2, width/2.2)
		rec3Dsnipet.set_ylim(-width/2.2, width/2.2)
		rec3Dsnipet.set_zlim(0, height)

		rec3Dsnipet.set_xticks([],[])
		rec3Dsnipet.set_yticks([],[])
		rec3Dsnipet.set_zticks([],[])

		# Get rid of colored axes planes
		# First remove fill
		rec3Dsnipet.xaxis.pane.fill = False
		rec3Dsnipet.yaxis.pane.fill = False
		rec3Dsnipet.zaxis.pane.fill = False

		rec3Dsnipet.set_axis_off()

		rec3Dsnipet.view_init(elev=-25, azim=100)
	
		mp = cm.ScalarMappable(cmap=cm.viridis)
		mp.set_array(N.hstack([vmin, vmax]))
		cbar = plt.colorbar(mp, orientation='horizontal', pad=-0.05, shrink=0.8)
		cbar.set_label(r'$\dot{q}^{\prime \prime}_\mathrm{in}$ [MW.m$^{-2}$]')

		# T ext 3D:
		# 3D receiver snipet with only flow_path in dense:
		rec3Dsnipet = fig.add_subplot(gs[3:5,0], projection='3d', aspect=1.4)
		rec3Dsnipet.set_facecolor('None')
		angres_per_bank = 1
		banks = N.shape(fluxmap)[1]
		fluxmap_fp = N.concatenate((fluxmap[:,banks/4:], fluxmap[:,:banks/4]), axis=1) # correction to start from east
		flux = N.hstack(T_ext)-273.15
		
		# Set-up the colormap:
		vmin=N.amin(flux)
		vmax=N.amax(flux)

		cNorm  = colors.Normalize(vmin, vmax)
		mp = cm.ScalarMappable(norm=cNorm, cmap=cm.hot)

		for i in xrange(len(flux)):
			if (i == N.hstack(fp[flow_path])).any():
				alpha = 1
				lw = 0.5
			else:
				alpha=0.1
				lw = 0
			ws = wh[i,0,:]
			hs = wh[i,1,:]
			for j in xrange(angres_per_bank):
				x = ws-width/2.
				y = ws*0

				z0 = N.ones(2)*hs[0]
				z1 = N.ones(2)*hs[1]

				x = N.vstack((x, x))
				y = N.vstack((y, y))
				z = N.vstack((z0, z1))

				face = mp.to_rgba(flux[i])
				surf = rec3Dsnipet.plot_surface(x, y, z, shade=False, antialiased=False, color=face, edgecolor=face, linewidth=lw, alpha=alpha)

	
		rec3Dsnipet.set_xlim(-width/2.2, width/2.2)
		rec3Dsnipet.set_ylim(-width/2.2, width/2.2)
		rec3Dsnipet.set_zlim(0, height)

		rec3Dsnipet.set_xticks([],[])
		rec3Dsnipet.set_yticks([],[])
		rec3Dsnipet.set_zticks([],[])

		# Get rid of colored axes planes
		# First remove fill
		rec3Dsnipet.xaxis.pane.fill = False
		rec3Dsnipet.yaxis.pane.fill = False
		rec3Dsnipet.zaxis.pane.fill = False

		rec3Dsnipet.set_axis_off()

		rec3Dsnipet.view_init(elev=-25, azim=100)
	
		mp = cm.ScalarMappable(cmap=cm.hot)
		mp.set_array(N.hstack([vmin, vmax]))
		cbar = plt.colorbar(mp, orientation='horizontal', pad=-0.05, shrink=0.8, ticks=N.linspace(520,760,7))
		cbar.set_label(r'${T_\mathrm{wall, ext}}$ [${^\circ}$C]')



		# flux fp plot:
		fluxplot = fig.add_subplot(gs[0,1])
		bank_lengths = pipe_lengths[flow_path]
		bank_lengths_2 = (bank_lengths[1:]+bank_lengths[:-1])/2.

		flux_lims = flux_limits(N.ones(len(fp[flow_path]))/n_tubes[fp[flow_path]]*m[flow_path], T_HC[flow_path], flux_limits_file)/1e3
		plt.plot(bank_lengths, flux_lims, color='r', label=r'${\dot{q}_\mathrm{LIMIT}^{\prime \prime}}$')
		#plt.plot(bank_lengths_2, flux_in[flow_path]/1e3, label=r'${\dot{q}_\mathrm{in}^{\prime \prime}}$', color='0')
		plt.plot(bank_lengths_2, q_net[fp[flow_path]]/areas[fp[flow_path]]/1e3, label=r'${\dot{q}^{\prime \prime}_\mathrm{net}}$', color='0.6')


		plt.legend()

		xticks = plt.gca().get_xticks()
		plt.xticks(xticks, [])
		plt.tick_params(direction='in')
		plt.xlim(N.amin(bank_lengths), N.amax(bank_lengths))
		plt.yticks(N.linspace(0,1200,7))

		plt.ylabel('Flux [kW.m$^{-2}$]')

		# Temperatures:
		Tplot = fig.add_subplot(gs[1,1])
		plt.plot(bank_lengths_2, T_ext[fp[flow_path]]-273.15, label=r'T${_\mathrm{abs}}$', color='0')
		plt.plot(bank_lengths_2, T_w_int[fp[flow_path]]-273.15, label=r'T${_\mathrm{wall,int}}$', color='0.4')
		plt.plot(bank_lengths, T_HC[flow_path]-273.15, label=r'T${_\mathrm{HC}}$', color='0.6')

		plt.legend()

		xticks = plt.gca().get_xticks()
		plt.xticks(xticks, [])
		plt.tick_params(direction='in')
		plt.xlim(N.amin(bank_lengths), N.amax(bank_lengths))

		plt.yticks(N.linspace(520,770,6))

		plt.ylabel('T [$^{\circ}$C]')

		# Fraction of allowable:
		Fracplot = fig.add_subplot(gs[2,1])
		
		plt.plot(bank_lengths_2, q_net[fp[flow_path]]/areas[fp[flow_path]]/1e3/((flux_lims[:-1]+flux_lims[1:])/2.), color='g')
		plt.ylabel(r'${\frac{\dot{q}_\mathrm{net}^{\prime \prime}}{\dot{q}_\mathrm{LIMIT}^{\prime \prime}}}$')
		
		xticks = plt.gca().get_xticks()
		plt.xticks(xticks, [])
		plt.tick_params(direction='in')
		plt.xlim(N.amin(bank_lengths), N.amax(bank_lengths))
		plt.ylim(-0.1,1.1)

		# Pressure drop:
		pdroplot = fig.add_subplot(gs[3,1])
		Dploc = -N.add.accumulate(Dp[fp[flow_path]])/1e5
		plt.plot(bank_lengths_2, Dploc, color='0')

		xticks = plt.gca().get_xticks()
		plt.xticks(xticks, [])
		plt.tick_params(direction='in')
		plt.xlim(N.amin(bank_lengths), N.amax(bank_lengths))

		plt.ylabel('${p(\ell)-p(0)}$ [bar]')

		# Velocity:
		vplot = fig.add_subplot(gs[4,1])
		plt.plot(bank_lengths_2, V[fp[flow_path]], color='0.5')

		plt.xlim(N.amin(bank_lengths), N.amax(bank_lengths))

		plt.xlabel('Flow path length (${\ell}$) [m]')
		plt.ylabel('${v_\mathrm{Na}}$ [m.s$^{-1}$]')


		plt.savefig(open(saveloc+'_fp_%s.png'%str(flow_path),'w'), dpi=400)
		plt.clf()
		plt.close(fig)

def pipes_differential_heating(flux_table_detailed, receiver):

	fileo = open(receiver,'r')
	data = pickle.load(fileo)
	areas = data['areas']
	fp = data['fp']
	n_banks = data['n_banks']
	n_elems = data['n_elems']
	n_tubes = data['n_tubes']
	wh_map = data['wh_map']

	#print areas, fp, n_tubes, wh_map, n_elems
	fluxmap = N.loadtxt(flux_table_detailed, skiprows=7, delimiter=',', usecols=N.arange(1,n_banks*int(n_tubes[0])+1))[::-1]
	fluxmap = N.array(fluxmap*1000.)
	
	Q_tubes = N.sum(fluxmap*areas[0]/n_tubes[0], axis=0)
	plt.figure()
	plt.subplot(211)

	plt.plot(range(len(Q_tubes)), Q_tubes/1e3)
	plt.vlines(N.linspace(0,n_banks*n_tubes[0], n_banks+1)-.5, 0, 50)

	plt.xlim(0,n_banks*n_tubes[0]-1)
	plt.ylim(0,50)
	plt.ylabel('Incident energy (kW)')

	plt.subplot(212)
	T_HC = data['T_HC']
	q_net = data['q_net']

	banks = [[3,2,1,0], [4,5,6,7]]
	for f in range(len(fp)):
		for b in range(n_banks/2):
			T_inlets = T_HC[f][b*n_elems]
			T_outlet = T_HC[f][(b+1)*n_elems-1]
			#q_net_avg = N.sum(q_net[fp[f]][b*n_elems:(b+1)*n_elems-1])/n_tubes[0]
			start = int(banks[f][b]*n_tubes[0])
			end = int(start+n_tubes[0])
			if start>end:
				start, end = end, start

			Q_tubes_avg = N.sum(Q_tubes[start:end])/n_tubes[0]

			T_guess = T_inlets+(T_outlet-T_inlets)*(1.+(Q_tubes[start:end]-Q_tubes_avg)/Q_tubes_avg)-T_outlet
			plt.plot(range(start, end), T_guess)
	plt.vlines(N.linspace(0,n_banks*n_tubes[0], n_banks+1)-.5, -15, 15)

	plt.xlabel('Pipes')
	plt.ylabel('Approximated outlet ${\Delta}$T (K)')
	plt.xlim(0,n_banks*n_tubes[0]-1)
	plt.ylim(-15, 15)

	plt.savefig(path[0]+'/differential_energy.png', dpi=400)
	plt.clf()
	
