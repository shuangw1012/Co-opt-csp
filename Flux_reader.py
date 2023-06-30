import numpy as np
from sys import path
from matplotlib import *
import matplotlib.pyplot as plt
import matplotlib.cm as cm
plt.rcParams["font.family"] = "Times New Roman"

'''
This vtk is read from SOLSTICE-included cylinder, not from Blender.
Direction is: from East to N to W to S, from bottom to top.
'''


def read_data(folder,r_height,r_diameter,num_bundle,bins,flux_file=False,flux_map=False):
	elem_area=r_diameter*np.sin(np.pi/num_bundle)*r_height/bins
	simul='%s/vtk/simul' % folder
	#print simul
	with open(simul) as f:
		lines = f.readlines()

	# for cylinders
	receiver_name='cylinder'
	for line in lines:
		if receiver_name in line:
			k=lines.index(line)
	lines_loc=lines[k:-1]
	for line in lines_loc:
		if 'LOOKUP_TABLE' in line:
			b3=lines_loc.index(line) #(b3+1)is the index that energy starts
			break
	# store the flux into an array (triangle)
	num_cell=bins*num_bundle*2
	Flux_triangle=np.array([])
	for j in range(b3+1,b3+1+num_cell):
		Flux_triangle=np.append(Flux_triangle,(float(lines_loc[j].rstrip().split(' ', 1 )[0]))/1000.) #kW/m2	
	# store the flux into an array (element)
	num_elem=bins*num_bundle
	Flux_elem=np.zeros(num_elem,dtype=float)
	Energy=np.array([])
	for i in range(num_elem):
		Flux_elem[i]=0.5*(Flux_triangle[2*i]+Flux_triangle[2*i+1])
		Energy=np.append(Energy,(Flux_elem[i]*elem_area))
	
	# reverse the flux data
	Flux_rev=np.zeros(num_elem,dtype=float)
	for i in range(num_bundle):
		for j in range(bins):
			Flux_rev[bins*i+j]=Flux_elem[bins*(i+1)-j-1]
			#print i,j,bins*i+j,bins*(i+1)-j-1
	
	max_flux=np.max(Flux_rev)
	
	# fit to the style consistent with the receiver model
	Flux_blender=np.zeros(num_elem,dtype=float)
	Flux_blender[:int(0.25*num_elem)]=Flux_rev[int(0.75*num_elem):]
	Flux_blender[int(0.25*num_elem):]=Flux_rev[:int(0.75*num_elem)]
	
	
	Flux=np.arange(bins*(num_bundle+1),dtype=float).reshape(bins,num_bundle+1)
	#print Flux.shape
	for i in range(bins):
		Flux[i,0]=r_height-0.5*r_height/bins-r_height/bins*i 
		for j in range(1,num_bundle+1):
			Flux[i,j]=Flux_blender[bins*(j-1)+i]
			#print i,j,Flux[i,j],bins*(j-1)+i
	if flux_file==True:
		title=np.array(7*(num_bundle+1))
		title=np.array(['Flux intensity plot (units kW/m2)'])
		title=np.append(title,np.repeat('', num_bundle))
		title=np.append(title,['Max. flux',str(np.max(Flux_rev)),'kW/m2'])
		title=np.append(title,np.repeat('', num_bundle-2))
		title=np.append(title,['Min. flux',str(np.min(Flux_rev)),'kW/m2'])
		title=np.append(title,np.repeat('', num_bundle-2))
		title=np.append(title,['Ave. flux',str(np.sum(Flux_rev)/num_elem),'kW/m2'])
		title=np.append(title,np.repeat('', num_bundle-2))
		title=np.append(title,np.repeat('', num_bundle+1))
		title=np.append(title,['','Receiver_circumferential_position (0=North',' +East) [deg] -->'])
		title=np.append(title,np.repeat('', num_bundle-2))
		title=np.append(title,['Receiver_vertical_position_[m]'])
		title=np.append(title,np.repeat('', num_bundle))
						
		Flux=np.append(title, Flux)
		Flux=Flux.reshape(bins+7, num_bundle+1)
		for i in range(1,num_bundle+1):
			Flux[6,i]=180.-360./num_bundle*(i-1)-180./num_bundle
		flux_file='%s/flux-table.csv' % folder
		np.savetxt(flux_file, Flux, fmt='%s', delimiter=',')
	
	theta = np.linspace(-np.pi, np.pi, num_bundle+1)
	h = np.linspace(r_height, 0, bins)
	theta, h = np.meshgrid(theta, h)
	Flux=Flux[7:,1:].astype(np.float)
	norm2 = colors.Normalize(vmin=np.amin(Flux), vmax=np.amax(Flux))
	if flux_map==True:
		fig = plt.figure()
		ax1 = fig.add_subplot(111)
		plt.subplots_adjust(bottom=0.2)
		plt.xlabel('Augular angle (rad)',fontsize=20)
		plt.ylabel('Height (m)',fontsize=20)
		im = ax1.pcolormesh(theta, h, Flux,cmap=cm.plasma)
		cbar3=fig.colorbar(im, ax=ax1)
		cbar3.set_label('Flux (kW.m$^{-2}$)',fontsize=20)
		cbar3.ax.tick_params(labelsize=20) 
		ax1.tick_params(axis='both', which='major', labelsize=20)
		#plt.show()
		plt.savefig(open('%s/flux_map.png'%folder,'w'), dpi=400)
		plt.close(fig)
		
	return max_flux
	

if __name__=='__main__':
	folder=path[0]
	r_height=24.
	r_diameter=16.
	num_bundle=16
	bins=50 # the vertial binning number
	read_data(folder,r_height,r_diameter,num_bundle,bins,flux_file=True,flux_map=True)
