#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  9 18:35:52 2021

@author: luciescharre
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 18 12:47:46 2021

@author: luciescharre
"""

import numpy as np 
import math
from collections import Counter
from astroquery.jplhorizons import Horizons
from astropy.coordinates import FK5
from astropy.coordinates import FK4
import astropy.units as u
from astropy.coordinates import SkyCoord
import time


long_UKST = ' 149°03 58.0 E'
lat_UKST = '31°16 36.5 S'
long_deg = 149.0661
long_dir= 'E'

prefix1 = []
prefix2 = []
emulsion = []

plate_no = []

RA_hh =[]
RA_mm =[]
RA_t =[]

Dec_sign = []
Dec_dd = []
Dec_mm = []

UT_y =[]
UT_m =[]
UT_d =[]

LST_hh = []
LST_mm = []

Exp_mmm =[]
Exp_t =[]

line_len = []
field_no=[]

# will read -0 as +0, read the sign as a character, if there is a minus sign, need to reflect 
# math functions will take rad as input
        

with open('UKST catalog.txt', 'r') as f:
    for line in f:
        line_len.append(len(line))
        
        
        # originally 17886 plates in catalog 
        # first filter (no dispersive element) brings number down to 15851
        # no T quality plates brings it down to 14435
       
        #only take non-dispersive element plates
        #if line[8] == '*':
        if line[7] == ' ':  #no dispersive element    
                # plate catalog length = n
                # line length = n+1
                # index starts at 0
                # the new line character \n is also part of length and length starts at 1, wheras index starts at 0, 
                # thus we need to take length -2 as the last index
                # at len 58 we have the abcde grading, no T's
                
                if (57 <= len(line) <= 62 and 'T' not in line[56:]) or \
                    (len(line)> 62 and 'T' not in line[56:61]):
                    #print(line)
                    #line_len.append(len(line))
                    
                    prefix1.append(line[0])
                    prefix2.append(line[1])
                    emulsion.append(line[40:46])
                    plate_no.append((line[2:7]))
                    
                    field_no.append(line[15:20])
                    
                    # RA hhmmt
                    RA_hh.append(float(line[20:22])) 
                    RA_mm.append(float(line[22:24])) 
                    RA_t.append(float(line[24:25])) 
                    
                    # Dec +/-ddmm
                    Dec_sign.append(line[25:26])  
                    Dec_dd.append(float(line[26:28]))  
                    Dec_mm.append(float(line[28:30]))  
                    
                    
                    # UT yymmdd
                    if float(line[30:32])== 00 or float(line[30:32])==1 or float(line[30:32])==2:
                        UT_y.append(2000+float(line[30:32])) 
                    else:
                        UT_y.append(1900+float(line[30:32]))
                        
                    UT_m.append(float(line[32:34])) 
                    UT_d.append(float(line[34:36])) 
                    
                    
                    # LST hhmm
                    if line[36:37]==' ' and line[37:38]!=' ':
                        LST_hh.append(float(line[37:38]))
                    
                    elif line[36:38]=='  ':
                        LST_hh.append(0.0)
                    
                    else: 
                        LST_hh.append(float(line[36:38]))
                    
                    LST_mm.append(float(line[38:40]))  
            
                    
                    # Exposure Time mmmt
                    if line[52:53]==' 'and line[53:54]!=' ':
                        Exp_mmm.append(float(line[53:55]))
                    
                    elif line[52:54]=='  ' and line[54:55]!=' ' :
                        Exp_mmm.append(float(line[54:55]))
                        
                    elif line[52:55]=='   ':
                        Exp_mmm.append(0.0)
                    
                    else: 
                        Exp_mmm.append(float(line[52:55]))
                    
                    Exp_t.append(float(line[55:56]))
                 

#print(len(plate_no))


#define time function(s)
        
def UT_date_to_JD(y,m,d,t):
    d += t
    if m==1 or m==2:
        y-=1
        m+=12
    
    A = int(y/100)
    B = 2-A+int(A/4)
    
    if y<0:
        C=int((365.25*y)-0.75)
    else:
        C=int(365.25*y)
        
    D = int(30.6001*(m+1))
    
    JD=B+C+D+d+1720994.5
    return JD
        
    
#test = UT_date_to_JD(1985,2,17,0.25)    
#print(test)


#have LST hhmm, want LST in decimal hours, (7)
#then GST in decimal hours using long_UKST, (15)
#convert to decimal hours in UT using JD date, (4,13)
#convert decimal hours to JD and add to the JD date (multiply by 1/24)
#want modified julian date at mid-exposure, add half of exposure
#might be better to have a function for this 

#want to find the LST in decimal hours including the mid exposure time

#def LST_exp_dec_hours(hh,mm,mmm,t):
#    mid_exp_t =(mmm+0.1*t)/2 #mid exposure time in minutes
#    LST_exp_t = hh + mm/60 + mid_exp_t/60  
#    return LST_exp_t       



# this gives time of beginning of exposure as opposed to at mid exposure
def LST_exp_dec_hours(hh,mm,mmm,t):
    #mid_exp_t =(mmm+0.1*t)/2 #mid exposure time in minutes
    LST_exp_t = hh + mm/60 
    return LST_exp_t      


def exp_t_h(mmm,t):
    exp_t_h =(mmm+0.1*t)/60 #mid exposure time in hours
    #LST_exp_t = hh + mm/60
    return exp_t_h  

#test = LST_exp_dec_hours(1,2,20,7)
#print(test)

#takes LST in decimal hours and the position to give GST in decimal hours
def LST_exp_to_GST(long_deg,long_dir,h):
    long_h = long_deg/15
    if long_dir == 'W':
        GST = h+long_h
        if GST > 24:
            GST -= 24
        elif GST < 0:
            GST += 24
        
    elif long_dir == 'E':
        GST = h-long_h
        if GST > 24:
            GST -= 24
        elif GST < 0:
            GST += 24  
    return GST   

#test = LST_exp_to_GST(64,'W',0.401453)
#print(test)


#take the GST in decimal hours and JD date to get UT in decimal hours 
def GST_to_UT(JD,GST):
    S = JD - 2451545.0
    T = S / 36525.0
    T0 = 6.697374558 + (2400.051336 * T)+ (0.000025862 * T**2)
    
    while T0 < 0:
        T0+=24
        
    while T0 > 24:
        T0-=24

    UT = GST-T0
    
    while UT < 0:
        UT+=24
        
    while UT > 24:
        UT-=24

    UT_h = UT * 0.9972695663 # UT in decimal hours
    #print(UT_h)
    UT_d =   UT_h/24                 #UT in decimal days, can be parsed in as t for the JD function
    return UT_d

#test = GST_to_UT(2444351.5,4.668119)
#print(test)

#perform subsequent time conversions
def cat_time_to_MJD(y,m,d,hh,mm,mmm,t,long_deg,long_dir):
    JD = UT_date_to_JD(y,m,d,0) # works
    LST_exp_t = LST_exp_dec_hours(hh,mm,mmm,t)  #works
    GST = LST_exp_to_GST(long_deg,long_dir,LST_exp_t)  #can't be sure
    UT = GST_to_UT(JD,GST)  #works
    JD = UT_date_to_JD(y,m,d,UT)  #works
    MJD = JD - 2400000.5
    return MJD


def cat_time_to_UT(y,m,d,hh,mm,mmm,t,long_deg,long_dir):
    JD = UT_date_to_JD(y,m,d,0) # works
    LST_exp_t = LST_exp_dec_hours(hh,mm,mmm,t)  #works
    GST = LST_exp_to_GST(long_deg,long_dir,LST_exp_t)  #can't be sure
    UT_frac = GST_to_UT(JD,GST)  #works
    UT = str(y)+str(m)+str(d+UT_frac)
    #JD = UT_date_to_JD(y,m,d,UT)  #works
    #MJD = JD - 2400000.5
    return UT

    
                
#define position functions  
#plate catalog, hours and minutes to decimal degrees
def RA_to_deg(hh,mm,t):
    RA_deg = (hh+(mm+0.1*t)/60)*15
    return RA_deg

#test=RA_to_deg(2,30,45)
#print(test)

##plate catalog, degrees and minutes to decimal degrees, pay attention to sign
def Dec_to_deg(sign,dd,mm):
    if sign == '-':
        Dec_deg = - (dd+mm/60)
    
    else:
        Dec_deg = dd+mm/60 
    return Dec_deg

#test=Dec_to_deg('-',100,20)
#print(test)
        

def deg_to_sexagesimal(deg):
    d = int(deg)
    m = int((deg - d) * 60)
    s = round((deg - d - m/60) * 3600,2)
    sexa_str = str(d) + ' ' + str(abs(m)) + ' ' + str(abs(s))
    return sexa_str

def deg_to_h(deg):
    h = int(deg/15)
    m = int((deg/15 - h) * 60)
    s = round((deg/15 - h - m/60) * 3600,2)
    sexa_str = str(h) + ' ' + str(m) + ' ' + str(s)
    return sexa_str

#define the tangent plate projection, takes RA and Dec values of centre of the plate (A,D) and other RA Dec values (alpha,delta)
#xi and eta are computed relative to the plate centre
def tangent_plane_projection(A,D,alpha,delta):
    cos_theta = math.sin(delta)*math.sin(D)+math.cos(delta)*math.cos(D)*math.cos(alpha-A)
    
    xi = math.cos(delta)*math.sin(alpha-A)/cos_theta     #make sure it takes them in radians
    
    eta = (math.sin(delta)*math.cos(D)-math.cos(delta)*math.sin(D)*math.cos(alpha-A))/cos_theta
    return xi, eta
    

#test = tangent_plane_projection(0,0,math.radians(-3),math.radians(-3))
#print(test)


#transform xi and eta from radians to mm using the plate scale,
#offset from the centre to one corner of the plate 
#by adding 177.5mm (call this transform [xi,eta] -> [x, y]) 
def transform_to_mm(xi,eta):   #
    plate_scale = 67.14  #arcsec/mm
    offset = 177.5 #mm
    xi_arcsec = xi*(3600 * 180)/math.pi
    eta_arcsec = eta*(3600 * 180)/math.pi
    
    # x and y positions relative to lower left corner, plate is 350x350mm
    x = xi_arcsec/plate_scale + offset
    y = eta_arcsec/plate_scale + offset
    return x,y


def error_projection(A,D,alpha,delta,RA_sigma,Dec_sigma):
    d_xi_d_alpha = (math.cos(delta) * (math.cos(D) * math.cos(delta) + math.cos(A - alpha) * math.sin(D) * math.sin(delta)))/(math.cos(D) * math.cos(A - alpha) * math.cos(delta) + math.sin(D) * math.sin(delta))**2
    d_xi_d_delta = (math.sin(D) * math.sin(A - alpha))/(math.cos(D) * math.cos(delta) * math.cos(A - alpha) + math.sin(D) * math.sin(delta))**2
    
    #print(RA_sigma)
    #print(d_xi_d_alpha)
    #print(d_xi_d_delta)
    
    xi_sigma = math.sqrt(( d_xi_d_alpha * RA_sigma )**2 + ( d_xi_d_delta * Dec_sigma )**2)
    #print(xi_sigma)
    
    d_eta_d_alpha = -(math.sin(A - alpha) * math.sin(2*delta))/(2 * (math.cos(D) * math.cos(A - alpha) * math.cos(delta) + math.sin(D) * math.sin(delta))**2)
    d_eta_d_delta = math.cos(A - alpha)/(math.cos(D) * math.cos(delta) * math.cos(A - alpha) + math.sin(D) * math.sin(delta))**2
       
    eta_sigma = math.sqrt(( d_eta_d_alpha * RA_sigma )**2 + ( d_eta_d_delta * Dec_sigma )**2)
    
    return xi_sigma, eta_sigma

def error_mm(xi_sigma,eta_sigma):
    plate_scale = 67.14  #arcsec/mm
     
    x_sigma = xi_sigma * 1/plate_scale * (3600 * 180)/math.pi
    y_sigma = eta_sigma   * 1/plate_scale * (3600 * 180)/math.pi 
    
    return x_sigma, y_sigma


#find limiting magnitude values    
def exp_time(mmm,t):
    exp_t =mmm+0.1*t
    return exp_t


#find the prediction direction and trail length
# JPL gives ra_rate*cosD and DEC_rate in arcsec/h
#delta should be given in radians, rough declination as given by JPL (degrees)
#delta in big loop already converted to rad
#exposure time should be in h

def location_prediction(RA_rate,DEC_rate,delta, exp_t):
    #RA_rate = RA_rate / math.cos(delta)
    
    #print(RA_rate)
    #print(DEC_rate)
    #print(exp_t)
    RA_length = RA_rate * exp_t  * math.pi/(180 * 3600)  #arcsec to radians
    DEC_length = DEC_rate * exp_t  * math.pi/(180 * 3600)  #arcsec to radians
    
    orient_angle = math.atan(RA_length/DEC_length)     #radians, angle from north to east
    arc_length = math.sqrt(RA_length ** 2 + DEC_length ** 2) * math.pi/(180 * 3600)  #arcsec to radians
    
    return orient_angle, RA_length, DEC_length,arc_length

# expect for R letter all filter and emulsions uniquely defined, 
# all single prefixes are always at the second index, 
# remember to include the condition that the first one has to be empty


# can compute the mag limit in same loop in 319 and only check against object in object loop
# i have exposure time and SNR

def mag_limit(i):
    # find the index in prefix list
    
    exp_t = Exp_mmm[i]+0.1*Exp_t[i] # find exposure time for the given plate index

    #print(prefix1[i])
    if prefix1[i]==' ' and prefix2[i] == 'B':  
        hb_exp_t = 60
        hb_mag_limit = 21
      
    elif prefix1[i]==' ' and prefix2[i] == 'J':  
        hb_exp_t = 60
        hb_mag_limit = 22.5
                 
    elif prefix1[i]==' ' and prefix2[i] == 'V':  
        hb_exp_t = 60
        hb_mag_limit = 21      
  
    elif prefix1[i]=='O' and prefix2[i] == 'R':  
        hb_exp_t = 60
        hb_mag_limit = 21     
        
    elif prefix1[i]==' ' and prefix2[i] == 'R':    #line[0:1]
        if emulsion[i] == '098-04':
            hb_exp_t = 60
            hb_mag_limit = 21   
            
        if emulsion[i] == 'IIIaF ':  #line[40:46]
            hb_exp_t = 90
            hb_mag_limit = 22
            
        else:
            hb_exp_t = 0.5
            hb_mag_limit = 1000
            
             
    elif prefix1[i]==' ' and prefix2[i] == 'I':  
        hb_exp_t = 90
        hb_mag_limit = 19.5    
    
    else: #for the filters I don't have data for create an artificially low limit value 
        hb_exp_t = 0.5
        hb_mag_limit = 1000
    
    if exp_t == 0:  #for the filters with no exposure time get a very high value, as no object will be found
        mag_limit = -1000

    else:
        flux_scaling_factor = (exp_t/hb_exp_t)**(1/2)
        safety = 1
        mag_limit = hb_mag_limit + 2.5 * math.log10(flux_scaling_factor) + safety
    
    return mag_limit


#computing the time, positions and magnitude limit for each plate
plates_MJD = []
plates_RA = []
plates_Dec = []
plates_mag_limit = []
plates_UT = []
plates_exp_t_h = []


for i in range(len(plate_no)):
    MJD = cat_time_to_MJD(UT_y[i] , UT_m[i] , UT_d[i] , LST_hh[i] , LST_mm[i] , Exp_mmm[i] , Exp_t[i] , long_deg , long_dir)
    

    UT_date = cat_time_to_UT(UT_y[i] , UT_m[i] , UT_d[i] , LST_hh[i] , LST_mm[i] , Exp_mmm[i] , Exp_t[i] , long_deg , long_dir)
    

    plates_MJD.append(MJD)
    plates_UT.append(UT_date)
    
    RA_deg = RA_to_deg(RA_hh[i] , RA_mm[i] , RA_t[i])
    Dec_deg = Dec_to_deg(Dec_sign[i] , Dec_dd[i] , Dec_mm[i])
        
    plates_RA.append(RA_deg)
    plates_Dec.append(Dec_deg)

    mag_lim = mag_limit(i)

    plates_mag_limit.append(mag_lim)

    
    plates_exp_t_h.append(exp_t_h(Exp_mmm[i] , Exp_t[i]))


plates_MJD=np.array(plates_MJD)
plates_RA=np.array(plates_RA)
plates_Dec=np.array(plates_Dec)



 # the following read in allows to isolate the sentry, esa and mpc objects    

"""object_ids = ['Ceres','Vesta','Pallas','Hygiea','Interamnia','Davida']"""

object_ids = []
sentry_objects = []
ESA_objects = []
MPC_objects = []    
        
with open('sentry 2.txt', 'r') as f:
    for line in f:
        tokens = line.split('   ')
        object_ids.append(str(tokens[0]).rstrip())
        sentry_objects.append(str(tokens[0]).rstrip())
        
        

with open('ESA.txt', 'r') as f:
    for line in f:
        object_ids.append(str(line[0:4]).rstrip()+' '+str(line[4:9]).rstrip())
        #ESA_objects.append(str(line[0:4]).rstrip()+' '+str(line[4:8]).rstrip())
        

          
ESA_sentry = list(dict.fromkeys(object_ids))

ESA_objects = []
for fruit in ESA_sentry:
    if fruit not in sentry_objects:
        ESA_objects.append(fruit)

object_ids = ESA_objects
#object_ids = sentry_objects

        
with open('MPC Emoid sorted 2.txt', 'r') as f:
    for line in f:
        tokens = line.split('\t')
        object_ids.append(str(tokens[0]).rstrip())
        #MPC_objects.append(line[27:37].rstrip())

ESA_MPC = list(dict.fromkeys(object_ids))
  
MPC_objects = []
for fruit in ESA_MPC:
    if fruit not in ESA_sentry:
        MPC_objects.append(fruit)
   
object_ids = MPC_objects    
#object_ids = ['1950 DA', '1999 RQ36', '1979 XB', '2007 FT3', '2004 MN4']



#split the plate in batches to parse into programmatic interface, 300 was too large
#270 before
step = 100
batches = np.arange(0,(len(plate_no)),step)
batches = np.append(batches,(len(plate_no)))

count=0
#for every object call out to interface in batches 
for celest_object in object_ids:
    #print(celest_object)
    object_MJD = []
    object_RA = []
    object_Dec = []
    object_mag = []
    object_RA_sigma = []
    object_Dec_sigma = []
    object_RA_rate  = []
    object_Dec_rate = []
    object_POS_3sigma = []
    
    for j in range(len(batches)-1):
        
        if celest_object == "499" or "Ceres" or 'Vesta':  #mars needs id type specified
            id_type = 'majorbody'
        else:
            id_type = "smallbody"
        
        eph = Horizons(id=celest_object, id_type=id_type, location='260', epochs=plates_MJD[batches[j]:batches[j+1]]).ephemerides(quantities='1,3,9,38',refsystem="B1950")
        
        
        object_Dec = np.append(object_Dec,eph['DEC'].data)
        object_RA = np.append(object_RA,eph['RA'].data)
        object_MJD = np.append(object_MJD,eph['datetime_jd'].data)
        
        #object_RA_sigma = np.append(object_RA_sigma,eph['RA_3sigma'].data)      #arcseconds
        #object_Dec_sigma = np.append(object_Dec_sigma,eph['DEC_3sigma'].data)   #arcseconds
        object_POS_3sigma = np.append(object_POS_3sigma,eph['RSS_3sigma'].data)
        object_mag = np.append(object_mag,eph['V'].data)
        object_RA_rate = np.append(object_RA_rate,eph['RA_rate'].data)
        object_Dec_rate = np.append(object_Dec_rate,eph['DEC_rate'].data)


#collecting arrays first is neater, then I can work with a singular index

    for i in range(len(plate_no)):
       
        if object_Dec[i] - plates_Dec[i] >= -4 and object_Dec[i] - plates_Dec[i] <= 4:
            
            A = math.radians(plates_RA[i])
            D = math.radians(plates_Dec[i])
            alpha = math.radians(object_RA[i])
            delta = math.radians(object_Dec[i])
            
            theta, RA_move, DEC_move,arclen = location_prediction(object_RA_rate[i] , object_Dec_rate[i] , delta , plates_exp_t_h[i])

            RA_fin = alpha + RA_move 
            DEC_fin = delta + DEC_move
            #initial position
            xi, eta = tangent_plane_projection(A,D,alpha,delta) #ephemeris in plate coordinates, centre is the zero point
            # final position
            xi_fin, eta_fin = tangent_plane_projection(A,D,RA_fin,DEC_fin) #ephemeris in plate coordinates, centre is the zero point
            

            
            if xi>= -0.05240777928304121 and xi <= 0.05240777928304121 and eta>= -0.05240777928304121 and eta <= 0.05240777928304121:

                if object_mag[i] <= plates_mag_limit[i] and  object_mag[i] <=24:                    
                    
                    alpha_sig = object_POS_3sigma[i] * math.pi/(180 * 3600)  #conversion arcsec to radians
                    delta_sig = object_POS_3sigma[i] * math.pi/(180 * 3600) 
                    
                    xi_sigma, eta_sigma = error_projection(A,D,alpha,delta,alpha_sig,delta_sig)
                    x_sigma,y_sigma = error_mm(xi_sigma,eta_sigma)


                    x,y = transform_to_mm(xi,eta) 
                
                    print(str(celest_object), end =' ')
                    print(','+str(prefix1[i])+str(prefix2[i]),',',str(plate_no[i]),',',str(field_no[i]),',',str(deg_to_h(object_RA[i]))+' '+str(deg_to_sexagesimal(object_Dec[i])),',',deg_to_h(object_POS_3sigma[i]/3600),' ',deg_to_sexagesimal(object_POS_3sigma[i]/3600),',',object_mag[i],',',plates_mag_limit[i])
                    
                    print("x: "+str(round(x,4))+" $\pm$ "+str(round(x_sigma,4)))
                    print("y: "+str(round(y,4))+" $\pm$ "+str(round(y_sigma,4)))
                    
                    
                    print(' ')
                    print('FK4 B1950')
                    print('START: ',deg_to_h(object_RA[i]),deg_to_sexagesimal(object_Dec[i]))
                    print('END: ' ,deg_to_h(math.degrees(RA_fin)), deg_to_sexagesimal(math.degrees(DEC_fin)))
                    print('RA error :',deg_to_h(object_POS_3sigma[i]/3600))
                    print('DEC error :',deg_to_sexagesimal(object_POS_3sigma[i]/3600))
                    

                    print('Orientation Angle: ',math.degrees(theta))
                    print('Arclength:',math.degrees(arclen))
                    
                    print('FK5 J2000')

                    c1 = SkyCoord(str(deg_to_h(object_RA[i])),str(deg_to_sexagesimal(object_Dec[i])), unit=(u.hourangle, u.deg), frame=FK4(equinox = 'B1950'))
                    c1=c1.transform_to(FK5(equinox='J2000'))  
                    print('START: '+c1.to_string('hmsdms'))
                    c2 = SkyCoord(deg_to_h(math.degrees(RA_fin)), deg_to_sexagesimal(math.degrees(DEC_fin)), unit=(u.hourangle, u.deg), frame=FK4)
                    
                    c2=c2.transform_to(FK5(equinox='J2000'))  
                    
                    print('END: '+c2.to_string('hmsdms'))
                    print(' ')

                 
                    
