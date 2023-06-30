import numpy as N

class Na():
	'''	
	Fink JK, Leibowitz L. Thermodynamic and Transport Properties of Sodium Liquid and Vapour. Reactor Engineering Division, Argonne National Laboratory; 1995.
	'''
	def __init__(self, force_Re_limits=False):
		self.Tmin = 371.
		self.Tmax = 1255.
		self.force_Re_limits = force_Re_limits

	def check_valid(self, T):
		if hasattr(T,'__len__'):
			Tlow = (T<self.Tmin).any()
			Thigh = (T>self.Tmax).any()
		else:
			Tlow = (T<self.Tmin)
			Thigh = (T>self.Tmax)
		if Tlow or Thigh:
			print ("Temperature of liquid sodium outside of correlations range")
			return False
		else:
			return True

	def rho(self, T):
		if self.check_valid(T):
			rho = 219.+275.32*(1.-T/2503.7)+511.58*N.sqrt(1.-T/2503.7) # kg/m3
			return rho
		else:
			return N.nan

	def mu(self, T):
		if self.check_valid(T):
			mu = N.exp(-6.4406-0.3958*N.log(T)+556.835/T)
			return mu
		else:
			return N.nan

	def k(self, T):
		if self.check_valid(T):
			k = 124.67-0.11382*T+5.5226e-5*T**2.-1.1842e-8*T**3.
			return k
		else:
			return N.nan

	def Cp(self, T):
		if self.check_valid(T):
			Cp = (1.6582-8.4790e-4*T+4.4541e-7*T**2.-2992.6*T**(-2.))*1e3
			return Cp
		else:
			return N.nan

	def h(self, T):
		if self.check_valid(T):
			h = (-365.77+1.6582*T-4.2395e-4*T**2.+1.4847e-7*T**3.+2992.6/T)*1e3 # kJ/kg
			return h
		else:
			return N.nan

	def s(self,T):
		if self.check_valid(T):
			s = (1.6582*N.log(T)-8.4790e-4*T+4.4541e-7/2.*T**2.+2992.6/2.*T**(-2.))*1e3
			return s
		else:
			return N.nan

	def V_tube(self, mf_HC_tube, T, D_tube_in):
		rho = self.rho(T)
		# Velocity in the tubes:
		V = mf_HC_tube/(N.pi*(D_tube_in/2.)**2.*rho)
		return V

	def h_conv_tube(self, mf_HC_tube, T, T_w, D_tube_in):
		rho = self.rho(T)
		mu = self.mu(T)
		Cp = self.Cp(T)
		k = self.k(T)
		# Velocity in the tubes:
		V = mf_HC_tube/(rho*N.pi*(D_tube_in/2.)**2.)
		# Reynolds:
		Re = D_tube_in*rho*V/mu
		# Prandtl:
		Pr = Cp*mu/k
		if (Pr<=0.1).all():
			if self.force_Re_limits == False:
				if N.logical_and((Re>=1e4), (Re<=1e6)).all():
					# Nusselt:
					Nu = 7.+0.025*(Re*Pr)**0.8
				else:
					print ('Re:', Re)
					print ('V:',V, '[m.s-1]')
					print ('Liquid sodium reynolds number outside of Lyon-Martinelli correlation range (1e4 <= Re <= 1e6)')
					return N.nan
			else:
				print ("Reynolds limits forced")
				print ('Maximum Reynolds:', N.amax(Re))
				# Nusselt:
				Nu = 7.+0.025*(Re*Pr)**0.8
				
		else:
			print (Pr)
			print ('Liquid sodium Prandtl number outside of Lyon-Martinelli correlation range(Pr <= 0.1)')
			return N.nan
		# Heat transfer coefficient:
		u = Nu*k/(D_tube_in)
		return u

	def h_conv_tube_V(self, V, T, T_w, D_tube_in):
		rho = self.rho(T)
		mu = self.mu(T)
		Cp = self.Cp(T)
		k = self.k(T)
		# Reynolds:
		Re = D_tube_in*rho*V/mu
		# Prandtl:
		Pr = Cp*mu/k
		if (Pr<=0.1).all():
			if N.logical_and((Re>=1e4), (Re<=1e6)).all():
				# Nusselt:
				Nu = 7.+0.025*(Re*Pr)**0.8
			else:
				print ('Liquid sodium reynolds number outside of Lyon-Martinelli correlation range (1e4 <= Re <= 1e6)')
				return N.nan
				
		else:
			print ('Liquid sodium Prandtl number outside of Lyon-Martinelli correlation range(Pr <= 0.1)')
			return N.nan
		# Heat transfer coefficient:
		u = Nu*k/(D_tube_in)
		return u

	def f_D(self, mf_HC_tube, T, D_tube_in, tube_roughness=45e-6): # Darcy Weissbach
		rho = self.rho(T)
		mu = self.mu(T)
		Cp = self.Cp(T)
		k = self.k(T)

		Pr = Cp*mu/k

		V = mf_HC_tube/(N.pi*(D_tube_in/2.)**2.*rho)
		Re = D_tube_in*rho*V/mu
		#f_D = (1.8*N.log10(Re)-1.5)**-2.
		S = N.log(Re/(1.816*N.log(1.1*Re/(N.log(1.+1.1*Re)))))
		f_D = (-2.*N.log10(tube_roughness/(3.71*D_tube_in)+2.18*S/Re))**(-2.) # Brkic using Lambert W-function approximation to solve Colebrook's implicit equation.
		return f_D

	def p_drop(self, mf_HC_tube, T, D_tube_in, tube_lengths, tube_roughness=45e-6):
		f_D = self.f_D(mf_HC_tube, T, D_tube_in, tube_roughness)
		rho = self.rho(T)
		V = mf_HC_tube/(N.pi*(D_tube_in/2.)**2.*rho)
		Dp = f_D*tube_lengths/D_tube_in*rho*V**2./2.
		return Dp

class Solar_salt():
	'''	
	H. Benoit et al. Renewable and Sustainable Energy Reviews 55 (2016) 298-315
	'''
	def __init__(self):
		self.Tmin = 533.
		self.Tmax = 873.

	def check_valid(self, T):
		if hasattr(T,'__len__'):
			Tlow = (T<self.Tmin).any()
			Thigh = (T>self.Tmax).any()
		else:
			Tlow = (T<self.Tmin)
			Thigh = (T>self.Tmax)
		if Tlow or Thigh:
			return False
		else:
			return True

	def rho(self, T):
		if self.check_valid(T):
			TC = T-273.15
			rho = 2090. - 0.636*TC
			return rho
		else:
			return N.nan

	def mu(self, T):
		if self.check_valid(T):
			TC = T-273.15
			mu = 2.2714e-2-1.2e-4*TC+2.281e-7*TC**2.-1.474e-10*TC**3.
			return mu
		else:
			return N.nan

	def k(self, T):
		if self.check_valid(T):
			TC = T-273.15
			k = 0.00019*TC + 0.443
			return k
		else:
			return N.nan

	def Cp(self, T):
		if self.check_valid(T):
			TC = T-273.15
			Cp = 1443. + 0.172*TC
			return Cp
		else:
			return N.nan

	def h(self, T):
		if self.check_valid(T):
			TC = T-273.15
			h = 1443.*TC+0.172/2.*TC**2.
			return h
		else:
			return N.nan

	def s(self, T):
		if self.check_valid(T):
			s = (1443.-0.172*273.15)*N.log(T) + 0.172*T
			return s
		else:
			return N.nan

	def V_tube(self, mf_HC_tube, T, D_tube_in):
		rho = self.rho(T)
		# Velocity in the tubes:
		V = mf_HC_tube/(N.pi*(D_tube_in/2.)**2.*rho)
		return V
		
	def h_conv_tube(self, mf_HC_tube, T, T_w, D_tube_in):
		rho = self.rho(T)
		mu = self.mu(T)
		Cp = self.Cp(T)
		k = self.k(T)
		# Velocity in the tubes:
		V = mf_HC_tube/(N.pi*(D_tube_in/2.)**2.*rho)
		# Reynolds:
		Re = D_tube_in*rho*V/mu
		# Prandtl:
		Pr = Cp*mu/k
		Nu = N.zeros(len(T))

		Wu_Pr = N.logical_and(Pr>=1.6,Pr<=23.9)
		Wu_trans_Re = N.logical_and(Re>=2300.,Re<=1e4)
		Wu_trans = N.logical_and(Wu_trans_Re,Wu_Pr)

		Wu_turb_Re = N.logical_and(Re>=1e4,Re<=4.3e4)
		Wu_turb = N.logical_and(Wu_turb_Re,Wu_Pr)

		Gnielinski_Pr = N.logical_and(Pr>=0.5,Pr<=2000.)
		Gnielinski_Re = N.logical_and(Re>=4e3,Re<=5e6)
		Gnielinski = N.logical_and(Gnielinski_Re,Gnielinski_Pr)

		valid = N.logical_or(N.logical_or(Wu_trans,Wu_turb),Gnielinski)

		if valid.all():

			Nu[Wu_trans] = 0.00154*Re[Wu_trans]**1.1*Pr[Wu_trans]**(1./3.) # Wu et al.

			Nu[Wu_turb] = 0.02948*Re[Wu_turb]**0.787*Pr[Wu_turb]**(1./3.) # Wu et al.

			f_D = self.f_D(mf_HC_tube, T, D_tube_in)
			Pr_w = self.Cp(T_w)*self.mu(T_w)/self.k(T_w)
			K = (Pr/Pr_w)**0.11
			Nu_3 = ((f_D/8.)*(Re-1000.)*Pr/(1.+12.7*N.sqrt(f_D/8.)*(Pr**(2./3.)-1.)))*K # Gnielinski
			Nu[Gnielinski] = Nu_3[Gnielinski]

		else:
			return N.nan
		# Heat transfer coefficient:
		u = Nu*k/(D_tube_in)
		return u

	def f_D(self, mf_HC_tube, T, D_tube_in, tube_roughness=45e-6):
		rho = self.rho(T)
		mu = self.mu(T)
		Cp = self.Cp(T)
		k = self.k(T)

		Pr = Cp*mu/k

		V = mf_HC_tube/(N.pi*(D_tube_in/2.)**2.*rho)
		Re = D_tube_in*rho*V/mu
		#f_D = (1.8*N.log10(Re)-1.5)**-2.
		
		S = N.log(Re/(1.816*N.log(1.1*Re/(N.log(1.+1.1*Re)))))
		f_D = (-2.*N.log10(tube_roughness/(3.71*D_tube_in)+2.18*S/Re))**(-2.) # Brkic using Lambert W-function approximation to solve Colebrook's implicit equation.
		return f_D

	def p_drop(self, mf_HC_tube, T, D_tube_in, tube_lengths):
		f_D = self.f_D(mf_HC_tube, T, D_tube_in)
		rho = self.rho(T)
		V = mf_HC_tube/(N.pi*(D_tube_in/2.)**2.*rho)
		Dp = f_D*tube_lengths/D_tube_in*rho*V**2./2.
		return Dp
	

if __name__=='__main__':
	import matplotlib.pyplot as plt
	plt.rc('font',size=8)
	sodium = Na()
	salt = Solar_salt()
	Tsalt = N.arange(salt.Tmin, salt.Tmax, 1.)
	Tsodium = N.arange(sodium.Tmin, sodium.Tmax, 1.)
	fig = plt.figure(figsize=(6./2.5,11.5/2.5))
	plt.subplots_adjust(top=0.82, right=0.95, hspace=0.15, bottom=0.08, left=0.25)
	alpha = 0.3

	plt.subplot(411)
	leg1, = plt.plot(Tsodium-273.15, sodium.rho(Tsodium), c='r')
	leg2, = plt.plot(Tsalt-273.15, salt.rho(Tsalt), c='b')
	ymin = 500
	ymax = 2000
	plt.fill_between([290,565],[ymax, ymax], [ymin, ymin], color='darkblue', alpha=alpha, lw=0)
	plt.fill_between([520,740],[ymax, ymax], [ymin, ymin], color='red', alpha=alpha, lw=0)
	#plt.xlabel('T [$^{\circ}$C]')
	plt.gca().set_xticklabels([])
	plt.ylabel(r'$\rho$ (kg${\cdot}$m${^{-3}}$)')
	plt.ylim(ymin, ymax)

	plt.subplot(412)
	plt.plot(Tsodium-273.15, sodium.Cp(Tsodium), c='r')
	plt.plot(Tsalt-273.15, salt.Cp(Tsalt), c='b')
	ymin = 1200
	ymax = 1600
	plt.fill_between([290,565],[ymax, ymax], [ymin, ymin], color='darkblue', alpha=alpha, lw=0)
	plt.fill_between([520,740],[ymax, ymax], [ymin, ymin], color='red', alpha=alpha, lw=0)
	#plt.xlabel('T [$^{\circ}$C]')
	plt.gca().set_xticklabels([])
	plt.ylabel(r'${Cp}$ (J${\cdot}$kg${^{-1}}$${\cdot}$K${^{-1}}$)')
	plt.ylim(ymin, ymax)

	plt.subplot(413)
	plt.plot(Tsodium-273.15, sodium.mu(Tsodium), c='r')
	plt.plot(Tsalt-273.15, salt.mu(Tsalt), c='b')
	ymin = 1e-4
	ymax = 5e-3
	plt.fill_between([290,565],[ymax, ymax], [ymin, ymin], color='darkblue', alpha=alpha, lw=0)
	plt.fill_between([520,740],[ymax, ymax], [ymin, ymin], color='red', alpha=alpha, lw=0)
	plt.yscale('log')
	#plt.xlabel('${T}$ ($^{\circ}$C)')
	plt.gca().set_xticklabels([])
	plt.ylabel(r'$\mu$ (kg${\cdot}$m${^{-1}}$${\cdot}$s${^{-1}}$)')
	plt.ylim(ymin, ymax)

	plt.subplot(414)
	plt.plot(Tsodium-273.15, sodium.k(Tsodium), c='r')
	plt.plot(Tsalt-273.15, salt.k(Tsalt), c='b')
	ymin = 0.1
	ymax = 1e2
	plt.fill_between([290,565],[ymax, ymax], [ymin, ymin], color='darkblue', alpha=alpha, lw=0)
	p1 = plt.Rectangle((0, 0), 0, 0, color='darkblue', alpha=alpha, lw=0)
	plt.fill_between([520,740],[ymax, ymax], [ymin, ymin], color='red', alpha=alpha, lw=0)
	p2 = plt.Rectangle((0, 0), 0, 0, color='red', alpha=alpha, lw=0)
	plt.yscale('log')
	plt.xlabel('${T}$ ($^{\circ}$C)')
	plt.ylabel(r'${k}$ (W${\cdot}$m${^{-1}}$${\cdot}$K${^{-1}}$)')
	plt.figlegend([leg2, leg1, p1, p2], ['Solar salt', 'Sodium', 'Low temperature system', 'High temperature system'], loc=1, ncol=1)
	plt.ylim(ymin, ymax)

	plt.savefig('/home/ael/Documents/Boulot/Material_properties/HC props.png', dpi=500)
