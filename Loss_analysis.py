from sys import path
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "Times New Roman"
from pandas import DataFrame
from sklearn import linear_model

#import statsmodels.api as sm
def receiver_correlation(r_diameter,r_height,folder):
	A_rec=np.pi*r_diameter*r_height
	R_data=np.loadtxt('%s/RELT.csv' % folder,delimiter=',', skiprows=1)
	abs_t=0.98
	ems_t=0.91
	eff_abs=abs_t/(2./np.pi*(1.-abs_t)+abs_t)
	eff_emi=ems_t/(2./np.pi*(1.-ems_t)+ems_t)
		
	# the reflecive loss: accurate
	Q_ref = (1-eff_abs)*R_data[:,0]
	
	# how to calculate T_ext_mean
	T_ext_mean=np.average(R_data[:,5])	
	
	# linear regression for T_ext
	T_ext={'X2': R_data[:,0], 'X3': (R_data[:,1]+273.15),'X4': R_data[:,2],'T_ext': R_data[:,5]}
	df = DataFrame(T_ext,columns=['X2','X3','X4','T_ext'])
	X = df[['X2','X3','X4']] 
	Y = df['T_ext']
	regr = linear_model.LinearRegression()
	regr.fit(X, Y)
	C0=regr.intercept_
	C=regr.coef_
	T_ext_linear=C0+C[0]*R_data[:,0]+C[1]*(R_data[:,1]+273.15)+C[2]*R_data[:,2]
	
	T_ext_4_mean=np.average(R_data[:,6])
	# linear regression for T_ext_4_mean
	T_ext_4={'X2': R_data[:,0], 'X3': (R_data[:,1]+273.15),'X4': R_data[:,2],'T_ext_4': R_data[:,6]}
	df = DataFrame(T_ext_4,columns=['X2','X3','X4','T_ext_4'])
	X = df[['X2','X3','X4']] 
	Y = df['T_ext_4']
	regr = linear_model.LinearRegression()
	regr.fit(X, Y)
	C1=regr.intercept_
	C1_1=regr.coef_
	T_ext_4_linear=C1+C1_1[0]*R_data[:,0]+C1_1[1]*(R_data[:,1]+273.15)+C1_1[2]*R_data[:,2]
	coefs_T=[C0,C,C1,C1_1]
	
	# how to calculate h_conv
	coefs=np.polyfit(R_data[:,2],R_data[:,7],4)
	h_conv=coefs[4]+coefs[3]*R_data[:,2]+coefs[2]*(R_data[:,2])**2+coefs[1]*(R_data[:,2])**3+coefs[0]*(R_data[:,2])**4 # h with polynominal fitting
	
	Q_emi=eff_emi*5.67e-8*A_rec*(T_ext_4_linear**4-(R_data[:,1]+273.15)**4)/1e6
	Q_conv=h_conv*A_rec*(T_ext_linear-(R_data[:,1]+273.15))/1e6
	Qnet=eff_abs*R_data[:,0]-Q_conv-Q_emi

	return coefs_T,coefs,eff_abs,eff_emi,A_rec
	
if __name__=='__main__':
	r_diameter=16.
	r_height=24.
	coefs_T,coefs,eff_abs,eff_emi,A_rec=receiver_correlation(r_diameter,r_height,folder=path[0])
	
	'''
	annual_info=np.loadtxt('%s/simulation_results_CA3.csv' % path[0] ,delimiter=',', skiprows=1)
	Rec_results=np.arange(len(annual_info)*12,dtype=float).reshape(len(annual_info),12)
	Rec_results[:,:-1]=annual_info
	
	Q_in = annual_info[:,3]*6764*12.2**2*annual_info[:,9]
	#Q_in = annual_info[:,7]
	
	Q_ref = (1-eff_abs)*Q_in/1e6
	
	T_ext_linear=coefs_T[0]+coefs_T[1][0]*Q_in/1e6+coefs_T[1][1]*(annual_info[:,4]+273.15)+coefs_T[1][2]*annual_info[:,5]
	T_ext_4_linear=coefs_T[2]+coefs_T[3][0]*Q_in/1e6+coefs_T[3][1]*(annual_info[:,4]+273.15)+coefs_T[3][2]*annual_info[:,5]
	
	#Q_emi=eff_emi*5.67e-8*A_rec*(T_ext_mean**4-(annual_info[:,4]+273.15)**4)/1e6
	Q_emi=eff_emi*5.67e-8*A_rec*(T_ext_4_linear**4-(annual_info[:,4]+273.15)**4)/1e6
	h_conv=coefs[4]+coefs[3]*annual_info[:,5]+coefs[2]*(annual_info[:,5])**2+coefs[1]*(annual_info[:,5])**3+coefs[0]*(annual_info[:,5])**4 # h with polynominal fitting
	Q_conv=h_conv*A_rec*(T_ext_linear-(annual_info[:,4]+273.15))/1e6
	Qnet=eff_abs*Q_in/1e6-Q_conv-Q_emi
	Eff_rec=np.zeros(len(annual_info))
	Diff=np.zeros(len(annual_info))
	for i in range(len(annual_info)):
		if Q_in[i]==0:
			continue
		#print int(annual_info[i,0]),int(annual_info[i,1]),int(annual_info[i,2]),annual_info[i,7]/1e6,Q_ref[i],Q_emi[i],Q_conv[i],Qnet[i],T_ext_linear[i]
		Eff_rec[i]=Qnet[i]/(Q_in[i]/1e6)
		Diff[i]=abs(annual_info[i,8]-Eff_rec[i])/(annual_info[i,8])
	#print (sum(annual_info[:,8]*Q_in/1e6)-sum(Eff_rec[:]*Q_in/1e6))/sum(annual_info[:,8]*Q_in/1e6)
	print sum(annual_info[:,8]*annual_info[:,7]/1e6),sum(Eff_rec[:]*Q_in/1e6)
	print (sum(annual_info[:,8]*annual_info[:,7]/1e6)-sum(Eff_rec[:]*Q_in/1e6))/sum(annual_info[:,8]*annual_info[:,7]/1e6)
	i = np.where(Diff==max(Diff))[0][0]
	#print annual_info[i,-2]/1e6,annual_info[i,-1],Eff_rec[i],Diff[i]
	
	Rec_results[:,-1]=Eff_rec
	
	title=np.array(['month', 'day', 'hour', 'DNI', 'T_abm', 'wind_speed', 'opt_eff', 'Qin', 'rec_eff','opt_eff_Cat','opt_eff_2D','rec_eff_Loss'])
	Rec_results=np.append(title, Rec_results)
	Rec_results=Rec_results.reshape(int(len(Rec_results)/12), 12)
	np.savetxt('%s/simulation_results_CA4.csv'%path[0], Rec_results, fmt='%s', delimiter=',')
	'''
