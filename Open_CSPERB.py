import numpy as N
from sys import path

def eval_v_max(e_net_fp, HC, T_in, T_out, W_abs, n_b, D_tube_o, D_tube_in, pipe_spacing):
	
	n_t = N.floor((W_abs/(n_b)/(D_tube_o+pipe_spacing)))
	Dh = HC.h(T_out)-HC.h(T_in)
	vmax = e_net_fp/(Dh*n_t*HC.rho(T_out)*N.pi/4.*D_tube_in**2.)
	vmin = HC.mu(T_out)/D_tube_o/HC.rho(T_out)*1e4
	#print e_net_fp/(2.44*(Dh*HC.rho(T_out)*N.pi/4.*D_tube_in**2.))*(D_tube_o+pipe_spacing)/N.pi
	
	return vmax, n_t

def determine_fp(total_power_incident, HC, T_in, T_out, D_tube_o, D_tube_in, n_b_max, W_abs, v_lim_max, v_lim_min=0., prism=True, bank_eff=1., min_fp=1, n_b_min=1, pipe_spacing=1e-3, even_fp=False):
	'''
	Screen possible flow configurations based on simplified heat transfer assumptions that lead to a conservative flow-velocity estimate:
	- Even distribution of flux on the banks
	- user defined pipe bank absorbed thermal efficiency
	Arguments:
	total_power_incident: incident power on the receiver
	HC: heat carrier instance
	T_in, T_out: outlet and inlet HC temperature
	n_b_max: maximum number of pipe banks considered
	D_tube_o, D_tube_i: outer and inner tube diameters
	W_abs: total width of the absorber, this is the cylinder perimeter for a cylindrical receiver.
	prism: if True, the actual absorber width is calculated from the width of cylindrical receiver banks rather than the full circular base perimeter
	v_lim_max: flow velocity limit in the receiver
	bank_eff: adjustment factor to account for bank absorbed-themal efficiency
	min_fp: minimum number of flow-paths to consider
	n_b_min: minimum number of banks to consider
	pipe_spacing: space between pipes in the banks.
	even_fp: imposes even number of flow-poaths if true.
	Returns:
	n_bs: Number of banks evaluated
	n_fps: simplest flow paths for the corresponding number of banks
	v_max: maximum velicity estimated
	n_ts: number of tubes per bank
	'''
	n_bs = []
	n_fps = []
	v_maxs = []
	n_ts = []
	# Screen the number of banks:
	for n_b in N.arange(n_b_min, n_b_max+1):
		# Check the possible flow-path configurations for this number of banks:
		n_fp = []
		for i in N.arange(min_fp,n_b+1):
			if ((n_b)%i)==0:
				if even_fp:
					if i%2:
						continue
				n_fp.append(i)
		n_fp = N.array(n_fp)
		print ('For %s banks, %s flow path configurations evaluated: %s'%(str(n_b),str(len(n_fp)), str(n_fp)))

		e_net_fp = bank_eff*total_power_incident/n_fp
		if prism:
			radius = W_abs/(2.*N.pi)
			bank_half_angle = N.pi/n_b
			W_abs_conf = n_b*2.*radius*N.sin(bank_half_angle)
		else:
			W_abs_conf = W_abs 
		v_max, n_t = eval_v_max(e_net_fp, HC, T_in, T_out, W_abs_conf, n_b, D_tube_o, D_tube_in, pipe_spacing)
		# Find the simplest configuration that respects the limit in velocity. The simplest configuration is the one that has the least flow-paths for the number of banks specified and that respects velocity constraints.
		idx_valid = N.logical_and(v_max<v_lim_max,v_max>v_lim_min)

		if idx_valid.any():
			idx_min_fp = N.argmin(n_fp[idx_valid])
			n_fp_min, v_max = str(n_fp[idx_valid][idx_min_fp]), v_max[idx_valid][idx_min_fp]

			print ('For %s banks, %s flow-paths is the simplest configuration and it has %s m/s maximum velocity'%(n_b, n_fp_min, v_max))
			n_bs.append(n_b)
			n_fps.append(n_fp_min)
			v_maxs.append(v_max)
			n_ts.append(n_t)
		else:
			print ('For %s banks, no valid flow-paths. Velocities:'%(n_b), v_max)
	return n_bs, n_fps, v_maxs, n_ts

class Bill_receiver():
	def __init__(self, width, height, n_banks, n_elems, D_tubes_o, D_tubes_i, abs_t, ems_t,  k_coating=None, D_coating_o=0.):

		wid = N.linspace(0., width, n_banks+1)
		hei = N.linspace(0., height, n_elems+1)

		self.width = width
		self.height = height

		self.n_banks = n_banks
		self.n_elems = n_elems

		# Rectangular mesh array. Convention is to take it starting East and counter-clockwise and starting from the bottom of the receiver.
		wh = N.zeros((n_banks*n_elems, 2, 2))

		wh[:,0,0], wh[:,0,1] = N.tile(wid[:-1], n_elems), N.tile(wid[1:], n_elems)
		wh[:,1,0], wh[:,1,1] = N.repeat(hei[:-1], n_banks), N.repeat(hei[1:], n_banks)

		self.wh = wh
		idxs = N.arange(N.shape(self.wh)[0], dtype=int)
		idxs = N.reshape(idxs, (n_elems, n_banks))
		self.wh_map = idxs # map the wh binning scheme indices on a 2D array fitting the fluxmap dimensions and orientation and respecting the wh order. Used for flow_path indexing.

		self.D_tubes_i = D_tubes_i
		self.D_tubes_o = D_tubes_o
		self.D_coating_o = D_coating_o

		self.areas = (self.wh[:,0,1]-self.wh[:,0,0])*(self.wh[:,1,1]-self.wh[:,1,0])
		self.abs_t = abs_t
		self.ems_t = ems_t

		self.eff_abs = abs_t/(2./N.pi*(1.-abs_t)+abs_t)
		self.eff_ems = ems_t/(2./N.pi*(1.-ems_t)+ems_t)

		self.k_coating = k_coating

	def flow_path(self, option='CTE2bu', sp_file=None, load=None):
		'''
		Flow-path options:
		- CTE2xy: Center to Edge, dual flow-paths.
			   x = 't'op or 'b'ottom for the injection point. Dual flow-paths.
			   y = 'u'nidirectional banks or 'a'lternating bank flow directions
		'''

		fluxmap = N.loadtxt(sp_file, skiprows=7, delimiter=',', usecols=N.arange(1,self.n_banks+1))[::-1]
		self.fluxmap = N.array(fluxmap*1000.)
		if load != None:
			self.fluxmap = fluxmap*load
		flatmap = N.hstack(self.fluxmap)

		if option[:4] == 'CTE2':

			n_fp=2
			passes = self.n_banks/n_fp

			self.fp = []
			self.flux_fp = []
			self.areas_fp = []

			for f in range(n_fp):

				flux_fp = N.zeros(N.shape(self.wh)[0]/n_fp)
				fp = N.zeros(N.shape(self.wh)[0]/n_fp, dtype=int)

				for p in range(passes):
					if f==0:
						start = int(N.floor(self.n_banks/2)-p-1)
					else:
						start = int(N.floor(self.n_banks/2)+p)

					end = start+self.n_elems

					fploc = self.wh_map[:,start]
					fluxloc = self.wh_map[:,start]

					# Check if alternatig directions
					if option[-1] == 'a':
						# Reverse bank direction if odd pass
						if p%2:
							fploc = fploc[::-1]
							fluxloc = fluxloc[::-1]

					# Reverse bank if top injection.
					if option[-2]=='t':
						fploc = fploc[::-1]
						fluxloc = fluxloc[::-1]
	
					fp[p*self.n_elems:(p+1)*self.n_elems] = fploc
					flux_fp[p*self.n_elems:(p+1)*self.n_elems] = flatmap[fluxloc]

				self.fp.append(fp)
				self.flux_fp.append(flux_fp)
				self.areas_fp.append(self.areas[fp])

	def balance(self, HC, material, T_in, T_out, T_amb, h_conv_ext, filesave='/home/charles/Documents/Boulot/ASTRI/ASTRI_Sandia_prototype/ref_case_result', load=1., air_velocity=5.):


		self.T_in = T_in
		self.T_out = T_out
		self.air_velocity = air_velocity

		self.q_net = N.zeros(len(self.wh))
		self.q_rad = N.zeros(len(self.wh))
		self.q_ref = N.zeros(len(self.wh))
		self.q_conv = N.zeros(len(self.wh))
		self.T_w_int = N.zeros(len(self.wh))
		self.T_ext = N.ones(len(self.wh))*self.T_in
		self.h_conv_int = N.zeros(len(self.wh))
		self.V = N.zeros(len(self.wh))
		self.n_tubes = N.zeros(len(self.wh))
		self.Dp = N.zeros(len(self.wh))
		self.pipe_lengths = []

		if h_conv_ext == 'WSVH':
			from Convection_loss import cyl_conv_loss_coeff_WSVH
			self.h_conv_ext = cyl_conv_loss_coeff_WSVH(self.height, 2.*self.radius, self.air_velocity, (self.T_in+self.T_out)/2., T_amb)
		if h_conv_ext == 'SK':
			from Convection_loss import cyl_conv_loss_coeff_SK
			self.h_conv_ext = cyl_conv_loss_coeff_SK(self.height, 2.*self.radius, self.D_coating_o/2., self.air_velocity, (self.T_in+self.T_out)/2., T_amb)
		else:
			self.h_conv_ext = h_conv_ext

		convergence_tot = N.ones(len(self.wh))
		while (convergence_tot>1e-6).any():
			self.h = []
			self.m = []
			self.T_HC = []
			T_ext_old = N.hstack(self.T_ext)
			for f in xrange(len(self.fp)):
				areas = self.areas_fp[f]
				fp = self.fp[f]

				n_elems_fp = len(fp)
				flux_fp = load*self.flux_fp[f]

				whloc = self.wh[fp]
				elem_lengths = N.abs(whloc[:,1,1]-whloc[:,1,0])

				h = N.ones(n_elems_fp+1)
				h[0] = HC.h(T_in)
				h[-1] = HC.h(T_out)
			
				# Initilaisation of teh flow-path
				q_ref = (1.-self.eff_abs)*flux_fp*areas
				q_net = self.eff_abs*flux_fp*areas
				T_ext = T_in+(T_out-T_in)*q_net/N.sum(q_net)
				T_HC = N.ones(len(q_net)+1)*T_in

				n_tubes = N.floor(areas/elem_lengths/self.D_coating_o)

				self.n_tubes[fp] = n_tubes

				convergence = N.ones(len(areas))

				while (convergence>1e-4).any():
					q_rad = self.eff_ems*areas*5.67e-8*(T_ext**4.-T_amb**4.)
					q_conv = areas*self.h_conv_ext*(T_ext-T_amb)
					q_net = (self.eff_abs*flux_fp*areas-q_rad-q_conv)

					m = N.sum(q_net)/(h[-1]-h[0])

					for i in xrange(len(areas)):

						h[i+1] = h[i]+q_net[i]/m
						T_next = T_HC[i]*(1.+q_net[i]/N.sum(q_net))
						h_next = HC.h(T_next)
						conv_h = 1.
						while conv_h>1e-7:
							k = HC.k((T_HC[i]+T_next)/2.)
							conduction = N.pi*(self.D_tubes_i/2.)**2./elem_lengths[i]*(-k*(T_HC[i]-T_HC[i+1]))
							h[i+1] = h[i]+(q_net[i]+conduction)/m
							T_next = T_HC[i]+(T_next-T_HC[i])*(h[i+1]-h[i])/(h_next-h[i])
							h_next = HC.h(T_next)
							conv_h = abs(h[i+1]-h_next)/h[i+1]

							T_HC[i+1] = T_next

					T_w_int = (T_HC[1:]+T_HC[:-1])/2.
					conv_T_int = N.ones(len(T_w_int))
					while (conv_T_int>1e-7).any():

						h_conv_int = HC.h_conv_tube(m/n_tubes, (T_HC[:-1]+T_HC[1:])/2., T_w_int, self.D_tubes_i)
						T_w_int_new = (T_HC[:-1]+T_HC[1:])/2.+q_net/n_tubes/(h_conv_int*N.pi*self.D_tubes_i/2.*elem_lengths)

						conv_T_int = N.abs(T_w_int-T_w_int_new)/T_w_int

						T_w_int = (T_w_int_new+T_w_int)/2.

					if self.k_coating == None:
						R_cond = N.log(self.D_tubes_o/self.D_tubes_i)/(material.k(T_w_int)*N.pi*elem_lengths)
					else:
						R_cond = (N.log(self.D_tubes_o/self.D_tubes_i)/material.k(T_w_int)+N.log(self.D_coating_o/self.D_tubes_o)/self.k_coating)/(N.pi*elem_lengths)

					T_ext_new = T_w_int+q_net/n_tubes*R_cond

					V = HC.V_tube(m/n_tubes, (T_HC[:-1]+T_HC[1:])/2., self.D_tubes_i)
					convergence = N.abs(T_ext-T_ext_new)/T_ext
					T_ext = N.copy(T_ext_new)

				# Pressure drops:
				#dists = N.sqrt((x[1:]-x[:-1])**2.+(y[1:]-y[:-1])**2.)
				#elem_lengths[:-1] += dists

				self.pipe_lengths.append(N.add.accumulate(N.hstack([0,elem_lengths])))

				self.h.append(h)
				self.m.append(m)

				self.T_HC.append(T_HC)

				self.q_net[fp] = q_net
				self.q_rad[fp] = q_rad
				self.q_ref[fp] = q_ref
				self.q_conv[fp] = q_conv

				self.T_w_int[fp] = T_w_int
				self.T_ext[fp] = T_ext
				self.h_conv_int[fp] = h_conv_int
				self.Dp[fp] = HC.p_drop(m/n_tubes, (T_HC[:-1]+T_HC[1:])/2., self.D_tubes_i, elem_lengths)
				self.V[fp] = V

			if h_conv_ext == 'WSVH':
				from Convection_loss import cyl_conv_loss_coeff_WSVH
				self.h_conv_ext = cyl_conv_loss_coeff_WSVH(self.height, 2.*self.radius, self.air_velocity, N.average(T_ext), T_amb)
			if h_conv_ext == 'SK':
				from Convection_loss import cyl_conv_loss_coeff_SK
				self.h_conv_ext = cyl_conv_loss_coeff_SK(self.height, 2.*self.radius, self.D_coating_o/2., self.air_velocity, N.average(T_ext), T_amb)
			else:
				self.h_conv_ext = h_conv_ext
			convergence_tot = N.abs(N.hstack(self.T_ext)-T_ext_old)/N.hstack(self.T_ext)

		import pickle

		data = {'wh': self.wh, 'width':self.width, 'height':self.height, 'n_banks':self.n_banks, 'n_elems':self.n_elems, 'D_tubes_o':self.D_tubes_o, 'D_tubes_i':self.D_tubes_i, 'eff_abs':self.eff_abs, 'abs_t':self.abs_t, 'eff_ems':self.eff_ems, 'ems_t':self.ems_t, 'k_t':material.k(self.T_w_int), 'wh_map':self.wh_map, 'fp':self.fp, 'areas':self.areas, 'areas_fp':self.areas_fp, 'HC':HC, 'T_in':self.T_in, 'T_out':self.T_out, 'h_conv_ext':self.h_conv_ext, 'h':self.h, 'm':self.m, 'flux_in':self.flux_fp, 'q_net':self.q_net, 'q_rad':self.q_rad, 'q_ref':self.q_ref, 'q_conv_ext':self.q_conv, 'T_amb':T_amb, 'T_HC':self.T_HC, 'T_w_int':self.T_w_int, 'T_ext':self.T_ext, 'h_conv_int':self.h_conv_int, 'V': self.V, 'fluxmap':self.fluxmap, 'n_tubes':self.n_tubes, 'Dp':self.Dp, 'pipe_lengths':self.pipe_lengths}

		file_o = open(filesave, 'w')
		pickle.dump(data, file_o)
		file_o.close()


class Cyl_receiver():

	def __init__(self, radius, height, n_banks, n_elems, D_tubes_o, D_tubes_i, abs_t, ems_t,  k_coating=None, D_coating_o=0., bank_geom='curved'):

		ang = N.linspace(0., 2.*N.pi, n_banks+1) # from east, counter-clockwise
		hei = N.linspace(0., height, n_elems+1)
		if len(hei>(n_elems+1)):
			hei = hei[:n_elems+1]

		rad = radius

		self.radius = radius
		self.height = height

		self.n_banks = n_banks
		self.n_elems = n_elems
		# Cylindrical mesh array. Convention is to take it starting East and counter-clockwise and starting from the bottom of the receiver.
		ahr = N.zeros((n_banks*n_elems, 3, 2))
		ahr[:,0,0], ahr[:,0,1] = N.tile(ang[:-1], n_elems), N.tile(ang[1:], n_elems)
		ahr[:,1,0], ahr[:,1,1] = N.repeat(hei[:-1], n_banks), N.repeat(hei[1:], n_banks)
		ahr[:,2,0], ahr[:,2,1] = rad, rad

		self.ahr = ahr
		idxs = N.arange(N.shape(self.ahr)[0], dtype=int)
		idxs = N.reshape(idxs, (n_elems, n_banks))
		self.ahr_map = idxs # map the ahr binning scheme indices on a 2D array fitting the fluxmap dimensions and orientation and respecting the ahr order. Used for flow_path indexing.

		self.D_tubes_i = D_tubes_i
		self.D_tubes_o = D_tubes_o
		self.D_coating_o = D_coating_o
		if bank_geom == 'curved':
			self.areas = (self.ahr[:,0,1]-self.ahr[:,0,0])*(self.ahr[:,1,1]-self.ahr[:,1,0])*self.ahr[:,2,0]
		elif bank_geom == 'flat':
			self.areas = 2.*self.ahr[:,2,0]*N.sin((self.ahr[:,0,1]-self.ahr[:,0,0])/2.)*(self.ahr[:,1,1]-self.ahr[:,1,0])
		self.abs_t = abs_t
		self.ems_t = ems_t

		self.eff_abs = abs_t/(2./N.pi*(1.-abs_t)+abs_t)
		self.eff_ems = ems_t/(2./N.pi*(1.-ems_t)+ems_t)

		self.k_coating = k_coating
		

	def flow_path(self, option='SENWS', fluxmap_file='/home/charles/Documents/Boulot/ASTRI/Sodium receiver_CMI/flux-table.csv', load=None):
		'''
		Organises the fluxmap in a series of flow-paths to solve the energy balance problem.
		The fluxmap from SolarPILOT is given as a 2D array from South counter-clockwise.
		Flow-path options:
		- SENWS: South-East-North-West-South the single path around the receiver injecting in South at the top and running around the profile vertically in a counter-clockwise fashion.
		- SEN-SWN: is a double flow path option. First flow-path geos from South to North via East, the second from South to North, via West. Injected at the top as well.
		- mvit+x: Multiple Vertical flow-paths introduced at the top and conducting as many vertical passes as needed to cover the profile. The rotation is counter-clockwise, the number of flow-paths is interpreted from the number at the end of the string argument.
		- smvSit+x: Symmetrical multiple vertical flow-paths south introduced at the top. Same as previous but all inlet is introduced on the south face and progresses towards the north in two groups. One group (even flow-paths) goes counter-clockwise, the other (odd flow-paths) goes clockwise
		- cmvSit+x: Crossed multiple vertical flow-paths from south and introduced at the top. Same as previous but all inlet is introduced on the south face and progresses until filling the south facing half-cylinder. Afterwards, the flow-paths are "crossed" (central symmetry using the cylinder axis) and the rest of the progression goes from the west and east towards North before exiting the receiever.
		- cmvNit+x: Crossed multiple vertical flow-paths from North and introduced at the top. Same as previous but all inlet is introduced on the south face and progresses until filling the south facing half-cylinder. Afterwards, the flow-paths are "crossed" (central symmetry using the cylinder axis) and the rest of the progression goes from the west and east towards North before exiting the receiever.
		'''
		fluxmap = N.loadtxt(fluxmap_file, skiprows=7, delimiter=',', usecols=N.arange(1,self.n_banks+1))[::-1] # reverse the vertical direction to start the fluxmap from the bottom of the receiver
		self.fluxmap = N.array(fluxmap*1000.)
		if load != None:
			self.fluxmap = self.fluxmap*load
		flatmap = N.hstack(self.fluxmap)

		if option == 'SENWS':

			flux_fp = N.zeros(N.shape(self.ahr)[0])
			fp = N.zeros(N.shape(self.ahr)[0], dtype=int)
			
			for b in xrange(self.n_banks):
				# Reverse bank direction if odd bank.
				fploc = self.ahr_map[:,b]
				fluxloc = self.ahr_map[:,b]
				if b%2:
					fploc = fploc[::-1]
					fluxloc = fluxloc[::-1]
				fp[b*self.n_elems:(b+1)*self.n_elems] = fploc
				flux_fp[b*self.n_elems:(b+1)*self.n_elems] = flatmap[fluxloc]

			self.fp = [fp]
			self.flux_fp = [flux_fp]
			self.areas_fp = [self.areas[fp]]

		if option == 'SEN-SWN':
			n_fp = 2

			self.fp = []
			self.flux_fp = []
			self.areas_fp = []

			for f in xrange(n_fp):
				flux_fp = N.zeros(N.shape(self.ahr)[0]/n_fp)
				fp = N.zeros(N.shape(self.ahr)[0]/n_fp, dtype=int)
				for b in xrange(self.n_banks/n_fp):
					if f == 0:
						fploc = self.ahr_map[:,b+self.n_banks/n_fp]
						fluxloc = self.ahr_map[:,b+self.n_banks/n_fp]
					else: # reverse the rotation
						fploc = self.ahr_map[:,self.n_banks/n_fp-1-b]
						fluxloc = self.ahr_map[:,self.n_banks/n_fp-1-b]
					if b%2: # Reverse bank direction if odd bank.
						fploc = fploc[::-1]
						fluxloc = fluxloc[::-1]

					fp[b*self.n_elems:(b+1)*self.n_elems] = fploc
					flux_fp[b*self.n_elems:(b+1)*self.n_elems] = flatmap[fluxloc]

				self.fp.append(fp)
				self.flux_fp.append(flux_fp)
				self.areas_fp.append(self.areas[fp])

		if option[:3] == 'mvi':
			# Multiple Vertical flow-paths Itroduced at the Top (mvit+x) going down and up, etc. x is part of the string and specifies the number of flow paths.
			self.fp = []
			self.flux_fp = []
			self.areas_fp = []
			if option[3] == 't':
				top_injection = True
			elif option[3] == 'b':
				top_injection = False
			nf = int(option[4:])
			if (self.n_banks%nf) != 0:
				print ('Mismatch between the flow path and the discretisation.')
				stop
			vpasses = self.n_banks/nf
			for f in xrange(nf):
				fp = N.zeros(len(self.areas)/nf, dtype=N.int16)
				flux_fp = N.zeros(len(self.areas)/nf)
				for i in xrange(vpasses):
					# Start and end of the flow-path segment in the fluxmap referential:
					strt = f*self.n_banks/nf+i
					end = strt+self.n_banks*self.n_elems

					elems = N.arange(strt, end, self.n_banks)

					if (i%2)==0:				
						elems = elems[::-1]
					if top_injection == False:
						elems = elems[::-1]

					# Idices of the flow-path segment in the ahr referential.
					fp[i*self.n_elems: (i+1)*self.n_elems] = elems
					flux_fp[i*self.n_elems: (i+1)*self.n_elems] = flatmap[elems]

				self.fp.append(fp)
				self.flux_fp.append(flux_fp)
				self.areas_fp.append(self.areas[fp])

		if option[:5] == 'smvSi':
			# Symmetrical multiple vertical flow-paths south introduced at the top "t" or bottom "b". Same as previous but all inlet is introduced on the south face and progresses towards the north in two groups. One group (even flow-paths) goes counter-clockwise, the other (odd flow-paths) goes clockwise.
			self.fp = []
			self.flux_fp = []
			self.areas_fp = []
			if option[5] == 't':
				top_injection = True
			elif option[5] == 'b':
				top_injection = False
			nf = int(option[6:])

			if (self.n_banks%nf) != 0:
				print ('Mismatch between the flow path and the discretisation.')
				stop
			elif (nf%2) != 0:
				print ('Error, ', nf, ' flow-paths. The number of flow-paths must be even for "smvSi".')
				stop

			vpasses = self.n_banks/nf		
			for f in xrange(nf):
				fp = N.zeros(len(self.areas)/nf, dtype=N.int16)
				flux_fp = N.zeros(len(self.areas)/nf)						
				for i in xrange(vpasses):
					if f%2:
						strt = self.n_banks-(f+1+i*nf)/2
						end = strt+self.n_banks*self.n_elems
					else:
						strt = (f+i*nf)/2
						end = strt+self.n_banks*self.n_elems

					elems = N.arange(strt, end, self.n_banks)

					if (i%2) == 0:
						elems = elems[::-1] # if even pass, go down.
					if top_injection == False:
						elems = elems[::-1]

					fp[i*self.n_elems: (i+1)*self.n_elems] = elems
					flux_fp[i*self.n_elems: (i+1)*self.n_elems] = flatmap[elems]

				self.fp.append(fp)
				self.flux_fp.append(flux_fp)
				self.areas_fp.append(self.areas[fp])


		if option[:5] == 'cmvSi':
			# Crossed multiple vertical flow-paths from south and introduced at the top. Same as previous but all inlet is introduced on the south face and progresses until filling the south facing half-cylinder. Afterwards, the flow-paths are "crossed" (central symmetry using the cylinder axis) and the rest of the progression goes from the west and east towards North before exitin the receiever.

			self.fp = []
			self.flux_fp = []
			self.areas_fp = []
			if option[5] == 't':
				top_injection = True
			elif option[5] == 'b':
				top_injection = False
			nf = int(option[6:])

			if (self.n_banks%nf) != 0:
				print ('Mismatch between the flow path and the discretisation.')
				stop
			elif (nf%2) != 0:
				print ('Error, ', nf, ' flow-paths. The number of flow-paths must be even for "cmvSit".')
				stop
			vpasses = self.n_banks/nf
			
			for f in xrange(nf):
				fp = N.zeros(len(self.areas)/nf, dtype=N.int16)
				flux_fp = N.zeros(len(self.areas)/nf)
				half_pass = int(N.ceil(vpasses/2.))
				for i in range(vpasses):
					if i< half_pass:
						if f%2:
							strt = self.n_banks-1-(f-1)/2*half_pass-i
							end = strt+self.n_banks*self.n_elems
						else:
							strt = f/2*half_pass+i
							end = strt+self.n_banks*self.n_elems
					else:
						if f%2:
							strt = self.n_banks/2-((f-1)/2+1)*(vpasses-half_pass)+(i-half_pass)
							end = strt+self.n_banks*self.n_elems
						else:
							strt = self.n_banks/2+(f/2+1)*(vpasses-half_pass)-(i-half_pass)-1
							end = strt+self.n_banks*self.n_elems

					elems = N.arange(strt, end, self.n_banks)

					if (i%2) == 0:
						elems = elems[::-1] # if even pass, go down.
					if top_injection == False:
						elems = elems[::-1]

					fp[i*self.n_elems: (i+1)*self.n_elems] = elems
					flux_fp[i*self.n_elems: (i+1)*self.n_elems] = flatmap[elems]

				self.fp.append(fp)
				self.flux_fp.append(flux_fp)
				self.areas_fp.append(self.areas[fp])

		if option[:5] == 'cmvNi':
			# Crossed multiple vertical flow-paths from North and introduced at the top. Same as previous but all inlet is introduced on the North face and progresses until filling the North facing half-cylinder. Afterwards, the flow-paths are "crossed" (central symmetry using the cylinder axis) and the rest of the progression goes from the west and east towards South before exiting the receiever.

			self.fp = []
			self.flux_fp = []
			self.areas_fp = []
			Strt=[]
			if option[5] == 't':
				top_injection = True
			elif option[5] == 'b':
				top_injection = False
			nf = int(option[6:])

			if (self.n_banks%nf) != 0:
				print ('Mismatch between the flow path and the discretisation.')
				stop
			elif (nf%2) != 0:
				print ('Error, ', nf, ' flow-paths. The number of flow-paths must be even for "cmvNit".')
				stop
			vpasses = self.n_banks/nf
			
			for f in xrange(nf):
				# For each flow-path, fille the flow-path list of sequential elements in which the HC goes in order.
				fp = N.zeros(len(self.areas)/nf, dtype=N.int16)
				flux_fp = N.zeros(len(self.areas)/nf)
				# for each pass, one per bank of pipe, find the elemenst of the fluxmap that are being seen by the fluid.
				half_pass = int(N.ceil(vpasses/2.))
				for i in range(vpasses):
					if i< half_pass:
						if f%2:
							strt = self.n_banks/2-1-(f-1)/2*half_pass-i
							end = strt+self.n_banks*self.n_elems
						else:
							strt = self.n_banks/2+f/2*half_pass+i
							end = strt+self.n_banks*self.n_elems
					else:
						if f%2:
							strt = self.n_banks-((f-1)/2+1)*(vpasses-half_pass)+(i-half_pass)
							end = strt+self.n_banks*self.n_elems
						else:
							strt = (f/2+1)*(vpasses-half_pass)-(i-half_pass)-1
							end = strt+self.n_banks*self.n_elems
					elems = N.arange(strt, end, self.n_banks)
					Strt.append(strt)
					if (i%2.)==0:
						elems = elems[::-1] # if even pass, go down.
					if top_injection == False:
						elems = elems[::-1]

					fp[i*self.n_elems: (i+1)*self.n_elems] = elems
					flux_fp[i*self.n_elems: (i+1)*self.n_elems] = flatmap[elems]

				self.fp.append(fp)
				self.flux_fp.append(flux_fp)
				self.areas_fp.append(self.areas[fp])
		return Strt
		
	def balance(self, HC, material, T_in, T_out, T_amb, h_conv_ext, filesave='/home/charles/Documents/Boulot/ASTRI/Sodium receiver_CMI/ref_case_result', load=1., air_velocity=5.):

		self.h = []
		self.m = []
		self.T_HC = []
		self.T_in = T_in
		self.T_out = T_out
		self.air_velocity = air_velocity

		self.q_net = N.zeros(len(self.ahr))
		self.q_rad = N.zeros(len(self.ahr))
		self.q_ref = N.zeros(len(self.ahr))
		self.q_conv = N.zeros(len(self.ahr))
		self.T_w_int = N.zeros(len(self.ahr))
		self.T_ext = N.ones(len(self.ahr))*self.T_in
		self.h_conv_int = N.zeros(len(self.ahr))
		self.V = N.zeros(len(self.ahr))
		self.n_tubes = N.zeros(len(self.ahr))
		self.Dp = N.zeros(len(self.ahr))
		self.pipe_lengths = []

		if h_conv_ext == 'WSVH':
			from Convection_loss import cyl_conv_loss_coeff_WSVH
			self.h_conv_ext = cyl_conv_loss_coeff_WSVH(self.height, 2.*self.radius, self.air_velocity, (self.T_in+self.T_out)/2., T_amb)
		if h_conv_ext == 'SK':
			from Convection_loss import cyl_conv_loss_coeff_SK
			self.h_conv_ext = cyl_conv_loss_coeff_SK(self.height, 2.*self.radius, self.D_coating_o/2., self.air_velocity, (self.T_in+self.T_out)/2., T_amb)
		else:
			self.h_conv_ext = h_conv_ext
		convergence_tot = N.ones(len(self.ahr))
		while (convergence_tot>1e-4).any():
			self.m = []
			self.T_HC = []
			T_ext_old = N.hstack(self.T_ext)
			for f in xrange(len(self.fp)):
				areas = self.areas_fp[f]
				fp = self.fp[f]
				n_elems_fp = len(fp)
				flux_fp = load*self.flux_fp[f]
				ahrloc = self.ahr[fp]
				elem_lengths = N.abs(ahrloc[:,1,1]-ahrloc[:,1,0])
				x = N.cos(N.sum(ahrloc[:,0], axis=1)/2.)*(N.sum(ahrloc[:,2], axis=1)/2.)
				y = N.sin(N.sum(ahrloc[:,0], axis=1)/2.)*(N.sum(ahrloc[:,2], axis=1)/2.)

				h = N.ones(n_elems_fp+1)
				h[0] = HC.h(T_in)
				h[-1] = HC.h(T_out)
			
				q_ref = (1.-self.eff_abs)*flux_fp*areas
				q_net = self.eff_abs*flux_fp*areas  # ??
				T_ext = T_in+(T_out-T_in)*q_net/N.sum(q_net)
				T_HC = N.ones(len(q_net)+1)*T_in

				n_tubes = N.floor(areas/elem_lengths/self.D_coating_o)

				self.n_tubes[fp] = n_tubes

				convergence = N.ones(len(areas))

				while (convergence>1e-4).any():
					q_rad = self.eff_ems*areas*5.67e-8*(T_ext**4.-T_amb**4.)
					q_conv = areas*self.h_conv_ext*(T_ext-T_amb)
					q_net = (self.eff_abs*flux_fp*areas-q_rad-q_conv)

					m = N.sum(q_net)/(h[-1]-h[0])

					for i in xrange(len(areas)):

						h[i+1] = h[i]+q_net[i]/m
						T_next = T_HC[i]*(1.+q_net[i]/N.sum(q_net))
						h_next = HC.h(T_next)
						conv_h = 1.
						ite1=0
						while conv_h>1e-8 and ite1<100000:
							k = HC.k((T_HC[i]+T_next)/2.)
							conduction = N.pi*(self.D_tubes_i/2.)**2./elem_lengths[i]*(-k*(T_HC[i]-T_HC[i+1]))
							h[i+1] = h[i]+(q_net[i]+conduction)/m
							T_next = T_HC[i]+(T_next-T_HC[i])*(h[i+1]-h[i])/(h_next-h[i])
							h_next = HC.h(T_next)
							conv_h = abs(h[i+1]-h_next)/h[i+1]
							ite1+=1
						T_HC[i+1] = T_next
					
					T_w_int = T_HC[1:]
					conv_T_int = N.ones(len(T_w_int))
					while (conv_T_int>1e-8).any():

						h_conv_int = HC.h_conv_tube(m/n_tubes, (T_HC[:-1]+T_HC[1:])/2., T_w_int, self.D_tubes_i)
						T_w_int_new = (T_HC[:-1]+T_HC[1:])/2.+q_net/n_tubes/(h_conv_int*N.pi*self.D_tubes_i/2.*elem_lengths)

						conv_T_int = N.abs(T_w_int-T_w_int_new)/T_w_int

						T_w_int = (T_w_int_new+T_w_int)/2.

					if self.k_coating == None:
						R_cond = N.log(self.D_tubes_o/self.D_tubes_i)/(material.k(T_w_int)*N.pi*elem_lengths)
					else:
						R_cond = (N.log(self.D_tubes_o/self.D_tubes_i)/material.k(T_w_int)+N.log(self.D_coating_o/self.D_tubes_o)/self.k_coating)/(N.pi*elem_lengths)

					T_ext_new = T_w_int+q_net/n_tubes*R_cond
					V = HC.V_tube(m/n_tubes, (T_HC[:-1]+T_HC[1:])/2., self.D_tubes_i)
					
					convergence = N.abs(T_ext-T_ext_new)/T_ext
					T_ext = (T_ext+T_ext_new)/2.
				# Pressure drops:
				#dists = N.sqrt((x[1:]-x[:-1])**2.+(y[1:]-y[:-1])**2.)
				#elem_lengths[:-1] += dists
				#print convergence,T_ext
				self.pipe_lengths.append(N.add.accumulate(N.hstack([0,elem_lengths])))
				
				self.h.append(h)
				self.m.append(m)
				self.T_HC.append(T_HC)

				self.q_net[fp] = q_net
				self.q_rad[fp] = q_rad
				self.q_ref[fp] = q_ref
				self.q_conv[fp] = q_conv

				self.T_w_int[fp] = T_w_int
				self.T_ext[fp] = T_ext
				self.h_conv_int[fp] = h_conv_int
				self.Dp[fp] = HC.p_drop(m/n_tubes, (T_HC[:-1]+T_HC[1:])/2., self.D_tubes_i, elem_lengths)
				self.V[fp] = V
			
			#print N.average(T_ext),self.h_conv_ext
			'''
			# Update convective loss:
			T_conv_av = N.average(T_ext)
			if ~N.isnan(T_conv_av):
				if h_conv_ext == 'WSVH':
					self.h_conv_ext = cyl_conv_loss_coeff_WSVH(self.height, 2.*self.radius, self.air_velocity, T_conv_av, T_amb)
				if h_conv_ext == 'SK':
					self.h_conv_ext = cyl_conv_loss_coeff_SK(self.height, 2.*self.radius, self.D_coating_o/2., self.air_velocity, T_conv_av, T_amb)
			'''
			if N.isnan(N.average(T_ext)): ########################
				h_conv_ext = 20.
			if h_conv_ext == 'WSVH':
				from Convection_loss import cyl_conv_loss_coeff_WSVH
				self.h_conv_ext = cyl_conv_loss_coeff_WSVH(self.height, 2.*self.radius, self.air_velocity, N.average(T_ext), T_amb)
			if h_conv_ext == 'SK':
				from Convection_loss import cyl_conv_loss_coeff_SK
				self.h_conv_ext = cyl_conv_loss_coeff_SK(self.height, 2.*self.radius, self.D_coating_o/2., self.air_velocity, N.average(T_ext), T_amb)
			else:
				self.h_conv_ext = h_conv_ext
			convergence_tot = N.abs(N.hstack(self.T_ext)-T_ext_old)/N.hstack(self.T_ext)
		if N.isnan(N.hstack(self.V).any()):
			print ('Energy balance error')
		else:
			print ('Energy balance OK')
		import pickle
		data = {'ahr': self.ahr, 'radius':self.radius, 'height':self.height, 'n_banks':self.n_banks, 'n_elems':self.n_elems, 'D_tubes_o':self.D_tubes_o, 'D_tubes_i':self.D_tubes_i, 'eff_abs':self.eff_abs, 'abs_t':self.abs_t, 'eff_ems':self.eff_ems, 'ems_t':self.ems_t, 'k_t':material.k(self.T_w_int), 'ahr_map':self.ahr_map, 'fp':self.fp, 'areas':self.areas, 'areas_fp':self.areas_fp, 'HC':HC, 'T_in':self.T_in, 'T_out':self.T_out, 'h_conv_ext':self.h_conv_ext, 'h':self.h, 'm':self.m, 'flux_in':self.flux_fp, 'q_net':self.q_net, 'q_rad':self.q_rad, 'q_ref':self.q_ref, 'q_conv_ext':self.q_conv, 'T_amb':T_amb, 'T_HC':self.T_HC, 'T_w_int':self.T_w_int, 'T_ext':self.T_ext, 'h_conv_int':self.h_conv_int, 'V': self.V, 'fluxmap':self.fluxmap, 'n_tubes':self.n_tubes, 'Dp':self.Dp, 'pipe_lengths':self.pipe_lengths}

		file_o = open(filesave, 'w')
		pickle.dump(data, file_o)
		file_o.close()
		
if __name__=='__main__':
	from HC import *
	from Tube_materials import Inconel740H

	height = 24.
	diameter = 16.
	n_banks = 16
	D_tubes_o = 42.16e-3
	D_tubes_i = D_tubes_o-2.*1.2e-3
	'''
	determine_fp(total_power_incident=0.88*620.e6, HC=Na(), T_in=520+273.15, T_out=740+273.15, D_tube_o=D_tubes_o, D_tube_in=D_tubes_i, 
	  n_b_max=n_banks, W_abs=N.pi*diameter, v_lim_max=2.44, v_lim_min=0., prism=False, bank_eff=1., min_fp=1, n_b_min=n_banks, pipe_spacing=1e-3, even_fp=False)
	'''
	rec = Cyl_receiver(radius=0.5*diameter, height=height, n_banks=n_banks, n_elems=50, D_tubes_o=D_tubes_o, D_tubes_i=D_tubes_i, abs_t=0.98, ems_t=0.91, k_coating=1.2, D_coating_o=D_tubes_o+45e-6)
	fmap_file = path[0]+'/flux-table.csv'
	save_file=path[0]+'/flux-table'
	rec.flow_path(option='cmvNit8',fluxmap_file=fmap_file)
	HC = Na()
	rec.balance(HC=HC, material=Inconel740H(), T_in=520+273.15, T_out=740+273.15, T_amb=20.+273.15, h_conv_ext='SK', filesave=save_file,air_velocity=5.)
	
	
	from Open_CSPERB_plots import *
	flux_limits_file='%s/201015_N07740_thermoElasticPeakFlux_velocity/N07740_OD%s_WT1.20_peakFlux_vel.csv'%(path[0],round(D_tubes_o*1000,2))
	tower_receiver_plots(files=save_file, efficiency=False, maps_3D=False, flux_map=False, flow_paths=True,saveloc=None, billboard=False, flux_limits_file=flux_limits_file)
	
	
