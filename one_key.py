import os
import re
import shutil
from cal_layout_r import *
from python_postprocessing import *
from cal_sun  import *
from matplotlib.patches import Circle
import matplotlib.pyplot as plt
import matplotlib as mpl
import random
import numpy as np
from sys import path
from scipy.optimize import curve_fit
from Deviation_aiming_new3 import *
from Open_CSPERB import *
from Open_CSPERB_plots import *
from HC import *
from Tube_materials import Inconel740H
from Flux_reader import *
from Loss_analysis import *


class one_key_start:
	def __init__(self, folder, subfolder, hst_w, hst_h, tower_h, num_hst, delta_r2,delta_r3,latitude,r_diameter,r_height,num_bundle):
		self.folder=folder # the main folder
		
		#self.subfolder='%s/SOLSTICE_opt'%self.folder
		self.subfolder=subfolder
		# for the heliostat field
		self.latitude=latitude # latitude of the field location
		self.hst_w=hst_w # heliostat width
		self.hst_h=hst_h # heliostat height
		self.delta_r2=delta_r2
		self.delta_r3=delta_r3
		dsep=0.
		self.num_hst=num_hst # number of heliostat
		self.DM=np.sqrt(hst_w**2+hst_h**2)+dsep 
		self.csv='%s/pos_and_aiming.csv'%(self.subfolder)
		self.csv_trimmed='%s/pos_and_aiming_trimmed.csv'%(self.subfolder)
		self.csv_aiming='%s/pos_and_aiming_new.csv'%(self.subfolder)
		
		# for the receiver
		self.tower_h=tower_h # tower height
		self.r_diameter=r_diameter # receiver diameter
		self.r_height=r_height # receiver height
		self.num_bundle=num_bundle # number of tube banks
		self.bins=50 # vertical binning of receiver surface
		self.num_rays=1000000

	def big_field_generation(self):
		# towerheight: the optical tower height
		pos_and_aim=radial_stagger(latitude=self.latitude, num_hst=self.num_hst, width=self.hst_w, height=self.hst_h, hst_z=7., 
		  towerheight=self.tower_h+0.5*self.r_height, R1=150., delta_r_g=[0.866,self.delta_r2,self.delta_r3,self.delta_r3], 
		  dsep=0., field='surround', savedir=self.subfolder, plot=False)
		aiming_cylinder(self.r_height,self.r_diameter, pos_and_aim, self.subfolder, c_aiming=0.)
		hst_info=np.loadtxt(self.csv,delimiter=',', skiprows=2)
		self.num_hst=len(hst_info)
		return self.subfolder
		
	def attenuation(self,csv):
		hst_info=np.loadtxt(csv,delimiter=',', skiprows=2)
		foc=hst_info[:,3]
		# to get the attenuation factor
		def func(x, b):
			return np.exp(-b * x)
		def fun_two(x):
			return 0.99321-0.0001176*x+1.97e-8*x**2
		xdata = np.linspace(0, np.max(foc), int(np.max(foc)*100))
		y = fun_two(xdata)
		ydata = y
		popt, pcov = curve_fit(func, xdata, ydata)
		y2 = [func(i, popt[0]) for i in xdata]
		att_factor =popt[0]
		return att_factor
		
	def equinox(self,csv_equinox):
		# run ray tracing simulation
		att_factor=self.attenuation(csv_equinox)
		dni=980.
		
		self.run_SOLSTICE(dni=dni,phi=0.,elevation=55.08,att_factor=att_factor,num_rays=self.num_rays,csv=csv_equinox)
		eta,q_results,eta_exc_intec=proces_raw_results('%s/vtk/simul'% self.subfolder,'%s/vtk'% self.subfolder)
		
		# calculate the efficiency of each heliostat
		hst_info=np.loadtxt(csv_equinox,delimiter=',', skiprows=2)
		num_hst=len(hst_info)
		Hst_data=np.arange(num_hst*10,dtype=float).reshape(num_hst,10)
		Hst_data[:,1:8]=hst_info
		q_hst_in,q_r_in=get_heliostat_to_receiver_data(simul='%s/vtk/simul'%self.subfolder, DNI=dni, receiver_name='cylinder')
		Hst_data[:,0]=np.arange(num_hst)
		Hst_data[:,8]=q_r_in
		Hst_data[:,9]=q_r_in[:]/q_hst_in[:] # heliostat optical efficiency
		
		return eta
		
	def get_I(self,elevation):
		I0=1363.
		zenith=90.-elevation
		AM=1./np.cos(zenith/180.*np.pi)
		I=I0*0.7**(AM**0.678)
		return I
	
	def run_SOLSTICE(self,dni,phi,elevation,att_factor,num_rays,csv): # the input is not a solstice style
		# transfer into SOLSTICE convention
		phi=270.-phi
		if phi > 360.:
			phi=phi-360.
		vtk_path='%s/vtk'%self.subfolder
		if os.path.exists(vtk_path):
			shutil.rmtree(vtk_path)
		file_path='%s/SOLSTICE.py' % self.subfolder
		old_file=file_path
		fopen=open(old_file,'r') 
		w_str=""
		for line in fopen:
			if re.search('dni_1=',line):
				line = 'dni_1=%s' % (dni) + '\n'
				w_str+=line
			elif re.search('azimuth_1=',line):
				line = 'azimuth_1=%s' % (phi) + '\n'
				w_str+=line
			elif re.search('zenith_1=',line):
				line = 'zenith_1=%s' % (elevation) + '\n'
				w_str+=line
			elif re.search('att_factor_1=',line):
				line = 'att_factor_1=%s' % (att_factor) + '\n'
				w_str+=line
			elif re.search('mainfolder_1=',line):
				line = "mainfolder_1='%s'" % (self.subfolder) + '\n'
				w_str+=line
			elif re.search('csv_1=',line):
				line = "csv_1='%s'" % (csv) + '\n'
				w_str+=line
			elif re.search('num_rays_1=',line):
				line = "num_rays_1=%s" % (num_rays) + '\n'
				w_str+=line
			elif re.search('r_cyl_1=',line):
				line = "r_cyl_1=%s" % (self.r_diameter/2.) + '\n'
				w_str+=line
			elif re.search('h_cyl_1=',line):
				line = "h_cyl_1=%s" % (self.r_height) + '\n'
				w_str+=line
			elif re.search('tower_h_1=',line):
				line = "tower_h_1=%s" % (self.tower_h) + '\n'
				w_str+=line
			elif re.search('num_bundle_1=',line):
				line = "num_bundle_1=%s" % (self.num_bundle) + '\n'
				w_str+=line
			else:
				w_str+=line
		wopen=open(old_file,'w')
		wopen.write(w_str)
		fopen.close()
		wopen.close()
		os.system('python %s/SOLSTICE.py ' % self.subfolder)
	
	def plot_hst_eff(self,DM,Hst_data,Vertex=[]):
		plt.rcParams['mathtext.default'] = 'rm'
		plt.rcParams['mathtext.fontset'] = 'stix'
		hst_info=Hst_data[:,1:8]
		Eff_hst=Hst_data[:,9]
		N_hst=len(Eff_hst)
		fig, ax = plt.subplots(figsize=(10,8))
		cmap = cm.rainbow
		norm = mpl.colors.Normalize(vmin = 0.3, vmax = 0.74)
		smap = cm.ScalarMappable(norm = norm, cmap = cmap)
		smap.set_array([])
		color_bar = fig.colorbar(mappable = smap, ax = ax, orientation = 'vertical',aspect=30)
		color_bar.set_label('$\it{\eta}_{\mathrm{a,hst}}$',fontsize=26)
		color_bar.ax.tick_params(labelsize=24) 
		for i in range(N_hst):
			cell = Circle(xy = (hst_info[i,0], hst_info[i,1]), radius=DM/2., edgecolor = 'black', 
			 linewidth=0.5,facecolor = smap.to_rgba(Eff_hst[i]))
			ax.add_patch(cell)
		plt.grid(linestyle='--')
		if Vertex!=[]:
			X_vertex=[]
			Y_vertex=[]
			for i in range(len(Vertex)):
				X_vertex.append(Vertex[i][0])
				Y_vertex.append(Vertex[i][1])
			plt.scatter(X_vertex,Y_vertex,s=10)
		plt.xlim(-1700.,1700.)
		plt.ylim(-1300.,1700.)
		ax.tick_params(axis='both', which='major', labelsize=24)
		ax.set_aspect(1)
		plt.savefig(self.subfolder+'/field.png', bbox_inches='tight',dpi=100)
		plt.close()
		plt.show()
	
if __name__=='__main__':
	
	folder=os.getcwd()
	latitude=34.85
	hst_w=12.2
	hst_h=12.2
	tower_h=175.
	num_hst=int(10500)
	num_bundle=int(16)
	
	
	
	r_diameter=16
	r_height=20
	delta_r2=0.91
	delta_r3=2.25
	
	subfolder='%s/SOLSTICE_run' % folder
	if not os.path.exists(subfolder):
		os.makedirs(subfolder)
	shutil.copy('%s/SOLSTICE.py' % folder, subfolder)
	
	Model=one_key_start(folder, subfolder, hst_w, hst_h, tower_h, num_hst, delta_r2,delta_r3,latitude,r_diameter,r_height,num_bundle)
	# generate the big field
	Model.big_field_generation()
	eta = Model.equinox('%s/pos_and_aiming.csv'%(subfolder))
	print (eta)
	#Model.annual_big_field()
	#Model.determine_field(target_aperture_power=620000000)
	#Model.flow_path(Qnet_demand=540.e6,velocity_limit=2.44)
	#Model.MDBA_aiming()
	#eta_tot,LCOE=Model.annual_trimmed_field()
	
	
