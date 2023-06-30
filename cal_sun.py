import numpy as N

class SunPosition:

    def __init__(self):
        pass


    def days(self, dd, mm):
        '''
        Arguement:
        dd - int, the day in the month
        mm - str, the month

        reference: J Duffie pp14, Table 1.6.1
        '''

        if mm=='Jan':
            days=dd

        elif mm=='Feb':
            days=31+dd

        elif mm=='Mar':
            days=59+dd

        if mm=='Apr':
            days=90+dd

        elif mm=='May':
            days=120+dd

        elif mm=='Jun':
            days=151+dd

        if mm=='Jul':
            days=181+dd

        elif mm=='Aug':
            days=212+dd

        elif mm=='Sep':
            days=243+dd

        if mm=='Oct':
            days=273+dd

        elif mm=='Nov':
            days=304+dd

        elif mm=='Dec':
            days=334+dd

        return days


    def declination(self, days, form='detail'):
        '''
        Reference: Solar Engineering of Thermal Processes, 4th edition, John A. Duffie and William A. Beckman, page 13

        declination angle: delta=23.45*sin(360*(284+day)/365)
 
        Arguement:
        day - int, day of the year (1-365)
        form - str, 'detail' or simple' model

        Return:
        delta - declination angle (deg)
        '''
        if form=='detail':
            B=float(days-1)*360./365.*N.pi/180.

            delta=(180./N.pi)*(0.006918 - 0.399912*N.cos(B) +0.070257*N.sin(B)- 0.006758*N.cos(2.*B) + 0.000907*N.sin(2.*B)- 0.002697*N.cos(3.*B) + 0.00148*N.sin(3.*B))

        else:
            delta=23.45*N.sin(360.*float(284+days)/365.*N.pi/180.) # deg

        return delta


    def solarhour(self, delta, latitude):

        '''        
        Reference: Solar Engineering of Thermal Processes, 4th edition, John A. Duffie and William A. Beckman, page 17

        Arguement:
        delta: declination angle- float, deg
        latitude: latitude angle: float, deg

        return: 
        hour: length of the daylight hour
        sunrise: the solar hour angle of the sunrise, deg

        '''

        sunset=N.arccos(-N.tan(latitude*N.pi/180.)*N.tan(delta*N.pi/180.))*180./N.pi # deg

        sunrise=-sunset

        hour=(sunset-sunrise)/15.
       
        return hour, sunrise


    def zenith(self, latitude, delta, omega):
        '''
        ref. eq.1.6.5
        Arguement:
        latitude: latitude angle, float, deg
        delta:  declination angle, float, deg
        omega: solar hour angle, float, deg

        return:
        theta: the zenith angle, float, deg

        '''        
        
        latitude*=N.pi/180.
        delta*=N.pi/180.
        omega*=N.pi/180.

        theta=N.arccos(N.cos(latitude)*N.cos(delta)*N.cos(omega)+N.sin(latitude)*N.sin(delta))*180./N.pi

        return theta
        
    def azimuth(self, latitude, theta, delta, omega):

        '''
        ref: eq. 1.6.6
        from South to West
        Arguement:
        latitude: latitude angle, deg
        delta: declination angle ,deg
        theta: zenith angle, deg
        omega: solar hour angle, deg

        return:
        phi: azimuth angle, deg, from South to west
        '''
        latitude*=N.pi/180.
        delta*=N.pi/180.
        theta*=N.pi/180.
    
        a1=N.cos(theta)*N.sin(latitude)-N.sin(delta)
        a2=N.sin(theta)*N.cos(latitude)
        b=a1/a2

        if abs(b+1.)<1e-10:
            phi=N.pi

        elif abs(b-1.)<1e-10:
            phi=0.
        else:
            phi=abs(N.arccos((N.cos(theta)*N.sin(latitude)-N.sin(delta))/(N.sin(theta)*N.cos(latitude)))) # unit radian

        if omega<0:
            phi=-phi

        phi*=180./N.pi

        return phi

    def convert_AZEL_to_declination_hour(self, theta, phi, latitude):

        '''
        Arguement:
        theta: zenith angle, deg
        phi: azimuth angle deg
        latitude: latitude latitude , deg

        return:
        delta: declination angle, deg
        omega: solar hour angle, deg
        '''      
        phi*=N.pi/180.
        theta*=N.pi/180.
        latitude*=N.pi/180.

        delta=N.arcsin(N.cos(theta)*N.sin(latitude)-N.cos(abs(phi))*N.sin(theta)*N.cos(latitude))

        omega=N.arccos((N.cos(theta)-N.sin(latitude)*N.sin(delta))/(N.cos(latitude)*N.cos(delta)))
        if phi<0:
            omega=-omega

        delta*=180./N.pi
        omega*=180./N.pi

        return delta, omega

    def convert_convention(self, tool):
        '''
        tool: str, 'Tracer', 'solstice' or else TBD
        '''
        #TODO develope the transfermation from conventions in different tools 
        pass

    def annual_angles(self, latitude, hemisphere, casefolder='NOTSAVE', nd=5, nh=5):
        '''
        Arguements:
        latitude: latitude of the location, deg
        hemisphere: str, 'North' or 'South'
        nd: number of points in the declination movement ( -- suggest nd>=5)
        nh: number of points in the solar hours (half of the day)
        '''
     
        # declination angle (deg)  
        # -23.45 ~ 23.45
        DELTA=N.linspace(-23.45, 23.45, nd)
        #print 'declination', delta

        # the maximum solar hour angle
        # summer solstice
        if hemisphere=='North':
            delta_summer=23.45
        elif hemisphere=='South':
            delta_summer=-23.45
        
        daymax, sunmax=self.solarhour(delta_summer, latitude)
        solartime=N.linspace(sunmax, 0, nh)
        #print 'solar', solar

        # solar time
        time=12.+solartime/15
        #print time

        table=N.zeros(((nh+3)*(nd+4)))
        table=table.astype(str)
        for i in xrange(len(table)):
            table[i]=' '

        table=table.reshape(nh+3,nd+4)
        table[0,6]='Declination (deg)'
        table[0,7]=''
        table[1,0]='Lookup table'
        table[1,4:]=N.arange(1,nd+1)
        table[2,2]='Solar time (h)'
        table[2,3]='Hour angle (deg)'
        table[2,4:]=DELTA
        table[3:,1]=N.arange(1,nh+1)
        table[3:,2]=time
        table[3:,3]=solartime

        c=1
        case_list=N.array(['Case','declination (deg)','solar hour angle (deg)', 'azimuth (deg) S-to-W ', 'zenith (deg)'])
        for i in xrange(len(DELTA)):
            delta=DELTA[i]
            hour, sunrise=self.solarhour(delta, latitude)
            sunset=-sunrise
            for j in xrange(len(solartime)):
                omega=solartime[j]
                if (omega>sunset or omega<sunrise):
                    table[3+j,4+i]='-' 
                else:
                    table[3+j, 4+i]=' case %s'%(c)
                    #zenith angle
                    theta=self.zenith(latitude, delta, omega)
                    # azimuth        
                    phi=self.azimuth(latitude, theta, delta, omega)

                    case_list=N.append(case_list, (c, delta, omega, phi, theta)) 
                    c+=1
                    
     
        case_list=case_list.reshape(len(case_list)/5,5)
        azimuth=case_list[1:,-1].astype(float)
        zenith=case_list[1:,-2].astype(float)

        if casefolder!='NOTSAVE':    
            N.savetxt(casefolder+'/table_view.csv', table, fmt='%s', delimiter=',')  
            N.savetxt(casefolder+'/annual_simulation_list.csv', case_list, fmt='%s', delimiter=',')          

        return azimuth, zenith

        


if __name__=='__main__':
	# example: PS10, summer solstice, solar noon
	latitude=34.85
	
	sun=SunPosition()
	'''
	dd=sun.days(21, 'Mar')
	delta=sun.declination(dd)
	print 'Declination angle', delta
	'''
	delta=-22.2927647419584

	daytime,sunrise=sun.solarhour(delta, latitude)
	
	
	omega=120. # solar noon
	theta=sun.zenith(latitude, delta, omega)
	phi=sun.azimuth(latitude, theta, delta, omega)
	
	sun.annual_angles(latitude, hemisphere='North', casefolder='.',nd=5, nh=9)
	
    
    



        
        
        
