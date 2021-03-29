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
long_deg = 116.865
long_dir= 'W'


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
        

with open('POSS-II catalog.txt', 'r') as f:
    for line in f:
        line_len.append(len(line))
        
        
        # originally 17886 plates in catalog 
        # first filter (no dispersive element) brings number down to 15851
        # no T quality plates brings it down to 14435
       
        #only take non-dispersive element plates
        if line[8:15] == 'POSS-II' and line[57] != 'R':
            if float(line[26:28]) >= 2.0:
            #if line[7] == ' ':      
                # plate catalog length = n
                # line length = n+1
                # index starts at 0
                # the new line character \n is also part of length and length starts at 1, wheras index starts at 0, 
                # thus we need to take length -2 as the last index
                # at len 58 we have the abcde grading, no T's
                
             #   if (57 <= len(line) <= 62 and 'T' not in line[56:]) or \
             #       (len(line)> 62 and 'T' not in line[56:61]):
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
                    
                    
                    if line[55:56] == ' ':
                        Exp_t.append(0.0)
                        
                    else:
                        Exp_t.append(float(line[55:56]))
                    
                 

print(len(plate_no))


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

#define the tangent plate projection, takes RA and Dec values of centre of the plate and other RA Dec values
#xi and eta are computed relative to the plate centre
def tangent_plane_projection(A,D,alpha,delta):
    cos_theta = math.sin(delta)*math.sin(D)+math.cos(delta)*math.cos(D)*math.cos(alpha-A)
    
    xi = math.cos(delta)*math.sin(alpha-A)/cos_theta     #make sure it takes them in radians
    
    eta = (math.sin(delta)*math.cos(D)-math.cos(delta)*math.sin(D)*math.cos(alpha-A))/cos_theta
    return xi, eta
    

#test = tangent_plane_projection(0,0,math.radians(-3),math.radians(-3))
#print(test)

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
    arc_length = math.sqrt(RA_length ** 2 + DEC_length ** 2) 
    
    return orient_angle, RA_length, DEC_length,arc_length

def arclen_angle_errors(x,y,dx,dy):
    d_arclen = math.sqrt((dx**2*x**2+dy**2*x**2)/(x**2+y**2))
    d_angle= math.sqrt((dx**2*y**2 + dy**2*dx**2)*(1/((x/y)**2+1)**2)/y**4)
    return d_arclen,d_angle

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


# expect for R letter all filter and emulsions uniquely defined, 
# all single prefixes are always at the second index, 
# remember to include the condition that the first one has to be empty


# can compute the mag limit in same loop in 319 and only check against object in object loop
# i have exposure time and SNR

def mag_limit(i):
    # find the index in prefix list
    
    exp_t = Exp_mmm[i]+0.1*Exp_t[i] # find exposure time for the given plate index

    
    if prefix1[i]==0 and prefix2[i] == 'B':  
        hb_exp_t = 60
        hb_mag_limit = 21
      
    elif prefix1[i]==0 and prefix2[i] == 'J':  
        hb_exp_t = 60
        hb_mag_limit = 22.5
                 
    elif prefix1[i]==0 and prefix2[i] == 'V':  
        hb_exp_t = 60
        hb_mag_limit = 21      
  
    elif prefix1[i]=='O' and prefix2[i] == 'R':  
        hb_exp_t = 60
        hb_mag_limit = 21     
        
    elif prefix1[i]==0 and prefix2[i] == 'R':    #line[0:1]
        if emulsion[i] == '098-04':
            hb_exp_t = 60
            hb_mag_limit = 21   
            
        if emulsion[i] == 'IIIaF ':  #line[40:46]
            hb_exp_t = 90
            hb_mag_limit = 22
             
    elif prefix1[i]==0 and prefix2[i] == 'I':  
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
    Dec_deg = Dec_to_deg(Dec_sign[i] , Dec_dd[i] , Dec_mm[i])
   
    MJD = cat_time_to_MJD(UT_y[i] , UT_m[i] , UT_d[i] , LST_hh[i] , LST_mm[i] , Exp_mmm[i] , Exp_t[i] , long_deg , long_dir)
    
    plates_MJD.append(MJD)
    
    
    
    UT_date = cat_time_to_UT(UT_y[i] , UT_m[i] , UT_d[i] , LST_hh[i] , LST_mm[i] , Exp_mmm[i] , Exp_t[i] , long_deg , long_dir)
    plates_UT.append(UT_date)   
    
    RA_deg = RA_to_deg(RA_hh[i] , RA_mm[i] , RA_t[i])
    
        
    plates_RA.append(RA_deg)
    
    plates_Dec.append(Dec_deg)
    
    plates_exp_t_h.append(exp_t_h(Exp_mmm[i] , Exp_t[i]))

plates_MJD=np.array(plates_MJD)
plates_RA=np.array(plates_RA)
plates_Dec=np.array(plates_Dec)


 # the following read in allows to isolate the sentry, esa and mpc objects    

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

        
with open('MPC Emoid sorted.txt', 'r') as f:
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
object_ids = ['1950 DA', '2007 FT3', '2001 CA21', '2007 DX40', '2017 SH33', '2020 TJ3', '2007 WP3', '2017 PY26', '2014 MV67', '2014 UD57', '2020 TY1', '2005 WG57', '2018 PZ21', '2010 TW149', '2012 WS3', '2013 RS43', '2009 HC', '2017 DA120', '2016 AZ193', '2014 CH13', '2017 UA45', '2014 QN266', '2020 YB5', '2004 FY3', '2018 BP3']
object_ids = ['2021 EC4', '2021 EX2', '2010 MA113', '2017 SM33', '2014 MG68', '2019 TJ5']
object_ids  = ['1999 AN10', '2000 EK26', '2009 WM1', '2011 AG5', '2018 LK', '2017 YZ1', '2006 SU49', '2004 MX2', '2002 NT7', '2000 QK130', '1998 SC15', '2011 JA', '2015 RE36', '2005 YU55', '2003 CR20', '2013 ED28', '1997 XF11', '2013 QR1', '2012 SW20', '2001 VK5', '2009 KK', '1999 JU3', '2011 SR5', '2004 GU9', '1993 KH', '2018 XG5', '2014 DM22', '2002 EZ11', '2012 HZ33', '2015 WA2', '2004 VC17', '2008 ER7', '2003 KO2', '2006 KV86', '2012 LK11', '2016 CB194', '2003 YK118', '2014 KZ45', '2004 VD17', '2001 XU30', '2006 FX', '2001 WN5', '2000 YG29', '2018 FO5', '2010 KR10', '2011 SD173', '2013 BP73', '2008 KZ5', '1999 DB7', '2007 CN26', '2009 XT6', '1998 HH49', '2002 PZ39', '2000 GV147', '2020 TF', '2003 MH4', '2011 WN15', '1999 RM45', '2018 JA', '2000 WO107', '2011 LT17', '2013 LM31', '1947 XC', '2018 BP', '2006 HQ30', '1998 VD35', '1998 QA1', '2002 QW47', '2011 YV62', '2018 DA1', '1989 FC', '2008 TC4', '2004 UE', '2018 JD2', '2018 TY2', '2004 KH17', '2007 LF', '1973 EA', '2002 JZ8', '2009 BE58', '2001 PM9', '2016 LN1', '1999 XL136', '2012 HG8', '2003 WH166', '2003 QO104', '2018 BM3', '2017 RV', '2008 XM', '2021 DM2', '2002 AW', '2010 JU39', '2019 OR1', '2012 GV17', '2001 SG286', '2002 JE9', '2011 UG20', '2003 RM10', '2005 SQ', '2010 FR', '2001 GQ2', '1998 KJ9', '1993 EA', '2005 UL5', '2017 EN', '2002 AJ129', '1988 TA', '2004 DV24', '2002 QF15', '2008 LV16', '2011 GN44', '2007 TB23', '2017 NR6', '1997 US2', '2005 GC120', '2016 EZ157', '2012 UZ33', '2012 CA21', '2008 DE', '2016 EU85', '2006 BC10', '2007 AE12', '2000 QS7', '2008 YU32', '1989 UP', '2003 YM137', '1989 AC', '2011 GA', '2011 CT4', '1999 YR14', '2012 HG31', '2020 FM6', '2017 MC4', '2017 CH1', '2005 YS8', '2002 LV', '2017 GM7', '2002 DU3', '2005 NB7', '2012 FQ1', '2010 DW1', '2019 PG1', '2000 EJ26', '1998 FW4', '2008 PF1', '2014 QK434', '2013 WT45', '2003 RN10', '2008 WQ63', '2017 VR12', '2012 TS78', '2015 RP82', '2006 GY2', '2016 HL', '2016 GD135', '2017 HT2', '2000 RS11', '2017 SV20', '2006 CU', '2003 YG118', '2014 HK129', '2007 LQ19', '2015 CV13', '2009 KC3', '2004 MD6', '2017 WV13', '1990 OS', '2018 HL2', '2004 HK33', '2016 WF9', '2020 DR2', '2019 BC1', '2004 HW', '2004 DC', '2019 TN', '2014 QL433', '1977 HA', '1998 ST27', '2006 BQ6', '2007 SQ6', '2009 ST19', '2005 TU50', '2006 GB', '2018 ED4', '2015 EO61', '2015 NU13', '2013 LE16', '2019 CE4', '2014 KP84', '2010 WZ8', '2000 MU1', '2002 AP3', '2019 XQ3', '2008 WZ13', '2010 RA147', '2011 SV71', '2000 SL10', '2000 GE2', '2008 BT18', '2005 UK1', '2002 XR14', '2020 BW14', '2009 CB3', '1998 HJ3', '2007 YQ56', '2020 BP13', '2010 CL19', '2003 YX1', '2007 HE15', '2009 DL46', '2019 WW4', '2000 EV70', '2002 UK11', '2004 TB10', '2019 WL4', '1990 UA', '2011 UH20', '1996 JG', '2016 CU193', '2014 ER49', '2000 CT101', '1994 CN2', '2018 EJ4', '2009 YG', '2013 FW13', '2011 TX8', '2016 WW9', '2001 WS1', '1999 NB5', '2018 FJ4', '2012 TV78', '2011 EM51', '2003 MA3', '1998 SF36', '2011 YH28', '1984 QA', '2001 ME1', '2005 GU', '2000 DO1', '2001 OY13', '2003 XB22', '2008 RG1', '1982 HR', '2009 EO2', '2010 TK54', '2005 WC1', '2000 PN9', '2006 OC5', '2005 GP21', '2018 GJ1', '2006 UK', '2012 VF37', '2010 XP69', '1986 JK', '2001 FC58', '2014 YL14', '2011 GS60', '2016 NV', '2008 OO', '2006 XG1', '2014 QO296', '2000 DP107', '2012 XD112', '2001 BO61', '2015 XE352', '2013 QU1', '2003 UX34', '1996 JA1', '2002 TB70', '2018 GM', '2020 FX4', '1998 OR2', '2013 BC18', '2014 GG17', '2002 KJ4', '2006 VW2', '2001 XP1', '2002 GT', '2017 VV', '2011 LC19', '2006 HC2', '2011 BN24', '2014 PS59', '1994 CC', '1998 QS52', '2000 GK137', '2008 TZ3', '1991 AQ', '2019 NY2', '1989 VB', '1990 MF', '2001 SN289', '2002 PE130', '2007 TL23', '1998 ML14', '2005 EJ225', '2006 SV19', '2001 SG10', '2005 AD13', '2005 PO', '1993 VD', '2011 DU', '2000 OL8', '2015 AR45', '2004 OT11', '2015 LG2', '2007 WV4', '2003 BD44', '2013 VO5', '2006 XD2', '2001 DF47', '2014 VL6', '2004 GA1', '2003 UW29', '2010 UQ7', '2010 LN14', '2005 EE', '2006 VC', '2008 DG5', '1997 QK1', '2011 UL21', '2000 YV137', '2020 UQ4', '1999 CF9', '1998 SH2', '1994 NE', '2004 UR1', '2005 JS108', '2013 US3', '2000 SP43', '2011 CP4', '2016 AH193', '1986 PA', '2007 CA19', '2005 VC', '1983 TB', '2017 OP19', '2019 VG6', '2011 BO24', '2002 WQ4', '2016 TJ18', '1998 SU4', '2010 UT7', '2006 WQ29', '2014 XJ3', '2006 WJ3', '2007 UL12', '2003 BK47', '1984 KD', '2008 VB1', '2002 YP2', '2001 KF54', '1991 BN', '2005 WA1', '2013 JM14', '2011 CC22', '2008 UW91', '1999 GK4', '2003 TK2', '2015 BN311', '2004 TB18', '2011 WO41', '2009 VZ', '2012 BN11', '2015 TE323', '2004 FJ11', '2006 UL217', '2015 PU228', '2008 UU1', '2003 AY2', '2015 DP155', '2001 XR30', '2001 RA12', '1989 JA', '2014 CD13', '2021 EW2', '2014 YB35', '1991 EE', '1996 EN', '2001 JV1', '2013 BK18', '2004 KB', '2004 TG10', '2011 BE38', '2017 RL', '2017 FE65', '2013 GS66', '2016 BZ14', '2005 LX36', '2011 OV18', '2014 PG51', '2013 TE135', '2004 AS1', '2004 CL', '2003 GY', '2014 GY48', '2004 PJ2', '2007 VZ137', '2008 HH', '2002 SM', '2007 ML24', '2000 AF6', '2013 WF108', '2015 BY310', '2004 RZ164', '2009 UN3', '2011 WV134', '1994 UG', '2016 XM', '2016 GT220', '2002 HK12', '1988 EG', '2017 BP31', '2008 SE85', '2016 CS29', '2019 EN', '2011 BX10', '2011 MU', '2006 CF', '1991 JX', '2011 YG6', '2016 CL32', '1994 PM', '2003 MU', '2011 JK', '2016 PF8', '2007 EF', '2018 WR1', '2000 QW69', '2010 CN44', '2001 GN2', '2013 BJ18', '2013 GH84', '1992 UY4', '2014 YW34', '1999 SL5', '2009 SG18', '2014 LJ21', '2000 GJ147', '1999 JM8', '1932 HA', '2001 UY4', '2007 FE', '2015 RR150', '1999 GL4', '2003 KU2', '1989 QF', '2013 PX6', '2002 OD20', '1991 VH', '2020 KO3', '2002 VR85', '2009 CC3', '2012 MJ6', '2015 UR51', '2007 MB24', '2016 FC15', '2003 RX7', '2009 AE16', '2001 YJ4', '2012 VO6', '2016 QJ44', '2017 UR4', '1998 SJ70', '2017 EQ13', '2010 DJ56', '2009 VW', '2008 UZ94', '2000 OH', '2017 HP3', '2004 LC2', '2004 LJ', '1993 DQ1', '2007 FA', '2002 TW55', '2002 NV16', '2000 JG5', '1959 LM', '2017 YR1', '2013 BO76', '2002 CD14', '2011 GQ61', '2011 BX18', '2018 MV4', '2013 RZ73', '2005 GJ8', '2008 CR118', '2016 WG9', '2001 QC34', '1987 SY', '2003 BH', '2017 UE45', '2001 PT9', '2018 KB1', '2008 SR1', '1998 OH', '2002 AC9', '2015 LK24', '2017 AQ19', '2014 SM260', '2005 JE46', '2016 PO39', '2005 YY128', '2000 QW7', '2018 AT', '2001 UZ16', '2016 KR3', '2015 FP118', '1951 RA', '2011 CG2', '2001 XQ', '2003 WG', '2003 WR21', '2014 MK55', '2004 WS2', '2015 DC200', '2006 UO', '2001 XU10', '2019 VL4', '2000 YN29', '2010 XC11', '2019 MF1', '2008 QY', '2003 YS17', '2009 EV', '2011 WO4', '2002 VU94', '2010 JL33', '2002 JB9', '2013 QC11', '2020 XD5', '2001 SG262', '2006 VD13', '2018 YH', '2003 WP21', '1990 SM', '2012 XY6', '2018 LZ13', '1994 CJ1', '2019 OM1', '2017 BK92', '2012 UR136', '2016 DP', '2013 TR132', '2002 AM31', '2007 JX2', '2012 VC82', '2001 UA5', '2011 UT20', '2010 CF19', '2006 SU131', '2009 EK1', '2008 XN', '2002 BM26', '2015 GZ13', '2002 VX94', '2015 HX176', '2014 OT111', '2018 XV5', '2004 JA27', '2005 EA', '2018 YJ2', '2016 CC30', '2001 SZ269', '1995 SA', '2013 WT44', '2003 AD23', '2009 HV2', '2012 MS4', '2001 YV3', '2017 AN19', '2019 XA3', '1998 WT', '2018 SX2', '1999 GU3', '2009 VA26', '2006 HZ51', '2011 XA3', '1998 QE2', '1989 UR', '2012 VO76', '2020 RH7', '2011 UV63', '1999 VF22', '2019 TP1', '2007 FF1', '2016 RN1', '2017 VX1', '2012 HJ1', '2015 TY237', '2020 QC6', '2000 XK47', '2010 JE87', '2006 SU19', '2004 FE5', '2015 TK178', '2013 CW32', '2009 CN5', '2013 JR28', '2001 HY7', '2016 XH2', '2011 AH37', '2002 FV5', '2001 LD', '2006 KD40', '2010 FH81', '2009 SG2', '2014 NE52', '2009 DZ42', '2003 CC', '2014 SR339', '2005 AY28', '1998 UT18', '2002 SR41', '1997 BQ', '1994 XL1', '2009 FU4', '2008 NO3', '2001 TA2', '2008 HS3', '2012 VL6', '1998 VF32', '1982 XB', '2013 UP8', '2003 YE45', '2010 TH19', '2006 BZ7', '2007 LD', '2014 UQ56', '2012 OO', '2017 OA', '2000 EU70', '2014 GE35', '2011 AK5', '2012 TP139', '2002 PM6', '2018 EM2', '2016 AL8', '1994 RC', '2015 FT118', '2002 EL6', '2008 EC69', '2007 AV2', '2003 NC', '2004 BW58', '2013 NF19', '2005 EK94', '2002 UQ3', '2018 NF15', '2011 OR15', '2007 VW137', '2005 FH', '2019 UX13', '2002 JV15', '2018 NY14', '2014 WD201', '1998 OK1', '2008 XW2', '2001 MG1', '2016 AC65', '2020 KT2', '1996 AW1', '2014 TA36', '2009 ES', '2017 TV4', '2013 LN31', '2017 NH', '1999 XA143', '2005 BY2', '2017 SY32', '2009 LW2', '2008 WZ94', '2019 LX4', '2013 CU83', '2001 XT1', '2011 CY46', '2017 HY50', '2002 FG7', '1998 HL1', '2002 GZ8', '2017 YT5', '2002 LX64', '2010 US7', '2018 BY2', '2001 GT2', '2003 EF54', '1996 FO3', '2009 XZ1', '1999 YG3', '2017 WC15', '2009 MS', '2013 VY13', '2017 ED4', '2020 MZ3', '2001 VB76', '2012 YO1', '2015 XP129', '2005 XL80', '2005 PY16', '2015 DE198', '1979 VA', '2005 MO13', '2010 VZ', '2011 KQ12', '2016 AA10', '2007 SJ', '2003 FG', '2010 AF30', '1989 DA', '2003 VE1', '2018 PQ23', '2015 VO66', '2016 BR80', '2015 XA352', '2004 RY10', '2011 EG17', '2003 RB', '2006 KY86', '2015 SY16', '1981 ET3', '2009 DZ', '2017 QW17', '2015 DB1', '1990 UQ', '2019 UK', '2003 DX10', '2012 QQ10', '2001 XN254', '2002 EV11', '2010 ON101', '2007 LE', '2016 OP5', '1992 SK', '2001 SO73', '2013 FX7', '2012 PS4', '2000 BF19', '2007 DL41', '2003 HB', '2001 JM1', '2005 WB1', '2012 TF53', '2007 YV29', '1987 SB', '1998 YM4', '2002 BK25', '2003 WD158', '2018 FS3', '2013 EW27', '2004 PS42', '2001 KY66', '2016 UW40', '2005 JF21', '2020 KN3', '2011 WU95', '2002 AC5', '2010 XY72', '2015 HY116', '2000 BO28', '1993 BX3', '2020 SL1', '2005 NE7', '2004 RY109', '2013 ER89', '2012 XS111', '2016 BT39', '1991 VK', '2012 LK9', '2000 ED104', '1999 JD6', '2015 UM67', '2017 MU8', '2011 SO32', '1999 VT25', '2009 KD5', '2019 UO9', '2006 TB', '2015 RA36', '2013 KN6', '2014 TW32', '2009 ST103', '2008 WN2', '2012 BT23', '2012 KC6', '2000 DK79', '2001 SN263', '2002 XP90', '2020 EN', '1999 VR6', '2008 JG', '2016 VS5', '2003 CA', '2003 UX5', '2015 XD130', '2000 WO67', '1999 SM5', '2014 SV141', '2015 NU2', '1999 RD32', '2009 QJ9', '2012 VE82']

object_ids = ['2000 EK26','2006 SU49','2011 AG5']
object_ids = ['1999 TU228','2000 SO238']
object_ids= ['2006 SU49','2001 QT202','4799 P-L', '2014 NF61','2000 AT147','2000 KM35']

#object_ids = ['Ceres','Vesta','Pallas','Hygiea','Interamnia','Davida']



#split the plate in batches to parse into programmatic interface, 300 was too large
step = 50
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

        eph = Horizons(id=celest_object, id_type=id_type, location='261', epochs=plates_MJD[batches[j]:batches[j+1]]).ephemerides(quantities='1,3,9,38',refsystem="B1950")

        
        object_Dec = np.append(object_Dec,eph['DEC'].data)
        object_RA = np.append(object_RA,eph['RA'].data)
        object_MJD = np.append(object_MJD,eph['datetime_jd'].data)
        object_mag = np.append(object_mag,eph['V'].data)

        object_POS_3sigma = np.append(object_POS_3sigma,eph['RSS_3sigma'].data)
        
        object_RA_rate = np.append(object_RA_rate,eph['RA_rate'].data)
        object_Dec_rate = np.append(object_Dec_rate,eph['DEC_rate'].data)

#collecting arrays first is neater, then I can work with a singular index

    for i in range(len(plate_no)):
       
        if object_Dec[i] - plates_Dec[i] >= -4 and object_Dec[i] - plates_Dec[i] <= 4:
            
            A = math.radians(plates_RA[i])
            D = math.radians(plates_Dec[i])
            alpha = math.radians(object_RA[i])
            delta = math.radians(object_Dec[i])
            
            xi, eta = tangent_plane_projection(A,D,alpha,delta) #ephemeris in plate coordinates, centre is the zero point
            # xi and eta in radians
            
            theta, RA_move, DEC_move,arclen = location_prediction(object_RA_rate[i] , object_Dec_rate[i] , delta , plates_exp_t_h[i])
            d_arclen,d_angle = arclen_angle_errors(RA_move,DEC_move,object_POS_3sigma[i]/3600,object_POS_3sigma[i]/3600)
            
            RA_fin = alpha + RA_move 
            DEC_fin = delta + DEC_move
            
            if xi>= -0.05240777928304121 and xi <= 0.05240777928304121 and eta>= -0.05240777928304121 and eta <= 0.05240777928304121:
               
                
                # if a certain filter was used and the exposure time is at a certain value, 
                # check limiting magnitude value against magnitude estimation
                     
                if object_mag[i] <= 24: #24 is max limit
                   
                    if object_POS_3sigma[i]/3600 < 3:

                        print(str(celest_object), end =' ')
                        print(','+str(prefix1[i])+str(prefix2[i]),',',str(plate_no[i]),',',str(field_no[i]),',',str(deg_to_h(object_RA[i]))+' '+str(deg_to_sexagesimal(object_Dec[i])),',',deg_to_h(object_POS_3sigma[i]/3600),' ',deg_to_sexagesimal(object_POS_3sigma[i]/3600),',',object_mag[i])

                        
                        # x,y positions on plate
                        #alpha_sig = object_POS_3sigma[i] * math.pi/(180 * 3600)  #conversion arcsec to radians
                        #delta_sig = object_POS_3sigma[i] * math.pi/(180 * 3600) 
                        #x,y = transform_to_mm(xi,eta)
                        #xi_sigma, eta_sigma = error_projection(A,D,alpha,delta,alpha_sig,delta_sig)
                        #x_sigma,y_sigma = error_mm(xi_sigma,eta_sigma)
                        #print("x: "+str(round(x,4))+" $\pm$ "+str(round(x_sigma,4)))
                        #print("y: "+str(round(y,4))+" $\pm$ "+str(round(y_sigma,4)))
                        #print(' ')

                 
                        print(' ')
                        print('FK4 B1950')
                        print('START: ',deg_to_h(object_RA[i]),deg_to_sexagesimal(object_Dec[i]))
                        print('END: ' ,deg_to_h(math.degrees(RA_fin)), deg_to_sexagesimal(math.degrees(DEC_fin)))
                        print('RA error :',deg_to_h(object_POS_3sigma[i]/3600))
                        print('DEC error :',deg_to_sexagesimal(object_POS_3sigma[i]/3600))
                        
    
                       
                        print(plates_UT[i])
                        print('Orientation Angle: ',str(round(math.degrees(theta),2))+'('+str(round(d_angle,2))+')')

                        
                        print('Arclength:',math.degrees(arclen))
                        print(d_arclen)

                        print('FK5 J2000')
    
                        c1 = SkyCoord(str(deg_to_h(object_RA[i])),str(deg_to_sexagesimal(object_Dec[i])), unit=(u.hourangle, u.deg), frame=FK4(equinox = 'B1950'))
                        c1=c1.transform_to(FK5(equinox='J2000'))  
                        print('START: '+c1.to_string('hmsdms'))

                        c2 = SkyCoord(deg_to_h(math.degrees(RA_fin)), deg_to_sexagesimal(math.degrees(DEC_fin)), unit=(u.hourangle, u.deg), frame=FK4)
                        
                        c2=c2.transform_to(FK5(equinox='J2000'))  
                        
                        print('END: '+c2.to_string('hmsdms'))
                        print(' ')
                        print(' ')