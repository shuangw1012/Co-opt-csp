import numpy as N
from CoolProp.CoolProp import PropsSI

def cyl_conv_loss_coeff_WSVH(height, diameter, velocity, T_wall, T_amb):
	# External cylinder loss from Winter, Sizmann and Vant Hull (the red book):
	g = 9.81 # gravity acceleration

	H = height # receiver height, m
	D = diameter # Receiver diameter, m
	v = velocity # velocity
	T_wall = T_wall # Average wall temperature K
	T_amb = T_amb # K

	T_m = (T_wall+T_amb)/2.

	rho = PropsSI('D', 'T', T_m, 'P', 101325., 'Air') # density
	mu = PropsSI('V', 'T', T_m, 'P', 101325., 'Air') # dynamic viscosity
	Cp = PropsSI('CP0MASS', 'T', T_m, 'P', 101325., 'Air') # specific enthalpy
	k = PropsSI('L', 'T', T_m, 'P', 101325., 'Air') # thermal conductivity
	B = PropsSI('isobaric_expansion_coefficient', 'T', T_m, 'P', 101325., 'Air') # bulk expansion coefficient

	nu = mu/rho # cinematic viscosity

	Re = rho*v*D/mu
	Pr = mu*Cp/k
	Gr = g*B*(T_wall-T_amb)*H**3./(nu**2.)

	# Forced convection:
	if Re<=1.:
		C = 0.
		m = 0.
	if 1.<=Re<4.:
		C = 0.891
		m = 0.330
	if 4.<=Re<40.:
		C = 0.821
		m = 0.385
	if 40.<=Re<4000.:
		C = 0.615
		m = 0.466
	if 4000.<=Re<40000.:
		C = 0.714
		m = 0.618
	if 40000.<=Re<4e6: # I increased the upper limit of the correlation from 4e5 to 4e6 to get something
		if Re>4e5:
			print "WARNING, the Reynolds (%s) number is %s beyond the validity range for existing forced convection correlations (4e5)."%(str(Re), str(Re-4e5))
		C = 0.0239
		m = 0.805

	Nu_f = C*Re**m*(0.785*T_wall/T_amb)**(m/4.)*1.167*Pr**0.45
	h_f = k*Nu_f/D

	# Natural convection:
	if 1e4<Gr*Pr<1e9:
		Nu_n = 0.59*(Gr*Pr)**0.25
	if 1e9<=Gr*Pr<1e14: # I increased the correlation limit from 1e12 to 1e13 to get something
		if Pr>1e12:
			print "WARNING, the Prandl number (%s) is %s beyond the validity range for existing natural convection correlations (1e12)."%(str(Pr), str(Pr-1e12))
		Nu_n = 0.13*(Gr*Pr)**0.333

	h_n = k*Nu_n/H

	h = N.pi/2.*(h_f**3.2+h_n**3.2)**(1./3.2)

	return h

def cyl_conv_loss_coeff_SK(height, diameter, pipe_radius, velocity, T_wall, T_amb):
	# External cylinder loss from Siebers and Kraabel (https://www.osti.gov/servlets/purl/6906848):
	g = 9.81 # gravity acceleration

	H = height # receiver height, m
	D = diameter # Receiver diameter, m
	v = velocity # velocity
	T_wall = T_wall # Average wall temperature K+6
	T_amb = T_amb # K
	
	T_m = (T_wall+T_amb)/2.

	rho = PropsSI('D', 'T', T_m, 'P', 101325., 'Air') # density
	mu = PropsSI('V', 'T', T_m, 'P', 101325., 'Air') # dynamic viscosity
	Cp = PropsSI('CP0MASS', 'T', T_m, 'P', 101325., 'Air') # specific enthalpy
	k = PropsSI('L', 'T', T_m, 'P', 101325., 'Air') # thermal conductivity
	B = PropsSI('isobaric_expansion_coefficient', 'T', T_m, 'P', 101325., 'Air') # bulk expansion coefficient

	nu = mu/rho # cinematic viscosity

	Re = rho*v*D/mu
	Pr = mu*Cp/k
	Gr = g*B*(T_wall-T_amb)*H**3./(nu**2.)

	# Forced convection:
	ks_D = pipe_radius/D

	def Nu_f_0(Re):
		return 0.3+0.488*Re**0.5*(1.+(Re/282000.)**0.625)**0.8

	def Nu_f_1(Re):
		return 0.0455*Re**0.81

	def interpolate(xs, ys, x):
		if ys[0] != ys[1]:
			a = (ys[1]-ys[0])/(xs[1]-xs[0])
			b = ys[1]/(a*xs[1])
			y = a*x+b
		else:
			y = ys[0]
		return y

	def Nu_f_interpolated(ks_D, Re):
		ks_D_data = [0., 75e-5, 300e-5, 900e-5]
		x0, x1 = N.amax(ks_D_data[ks_D_data<=ks_D]), N.amin(ks_D_data[ks_D_data>=ks_D])
		x_inter = [x0, x1]
		Nu_f_inter = []
		for x in x_inter:		
			if x == 0:
				Nu_f_inter.append(Nu_f_0(Re))
			if x == 75e-5:
				if Re<=7e5:
					Nu_f_inter.append(Nu_f_0(Re))
				if 7e5<Re<2.2e7:
					Nu_f_inter.append(2.57e-3*Re*0.98)
				if Re>=2.2e7:
					Nu_f_inter.append(Nu_f_1(Re))
			if x == 300e-5:
				if Re<=1.8e5:
					Nu_f_inter.append(Nu_f_0(Re))
				if 1.8e5<Re<4e6:
					Nu_f_inter.append(0.0135e-3*Re**0.89)
				if Re>=4e6:
					Nu_f_inter.append(Nu_f_1(Re))
			if x == 900e-5:
				if Re<=1e5:
					Nu_f_inter.append(Nu_f_0(Re))
				if Re>1e5:
					Nu_f_inter.append(Nu_f_1(Re))

		Nu_f = interpolate(x_inter, Nu_f_inter, ks_D)
		return Nu_f
	
	h_f = k*Nu_f_interpolated(ks_D, Re)/D

	# Natural convection:
	Nu_n = 0.098*Gr**(1./3.)*(T_wall/T_amb)**-0.14
	h_n = k*Nu_n/H
	#print height, diameter, pipe_radius, velocity, T_wall, T_amb,h_n,h_f
	h = (h_f**3.2+(N.pi/2.*h_n)**3.2)**(1./3.2)

	return h

if __name__ == '__main__':
	print cyl_conv_loss_coeff_SK(height=24., diameter=16., pipe_radius=30.15e-3, velocity=3., T_wall=520.+273.15, T_amb=25.+273.15)

