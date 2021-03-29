#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb  6 11:59:51 2021

@author: luciescharre
"""

#import urllib 
#import selenium


from selenium import webdriver
from selenium.webdriver.common.keys import Keys

from selenium.webdriver.support.ui import Select

import time
from selenium.webdriver.firefox.options import Options

options = Options()
options.headless = True


driver = webdriver.Firefox()


url = 'http://www-wfau.roe.ac.uk/sss/pixel.html'
#urllib.request.urlopen(url, data=None, cafile=None, capath=None, cadefault=False, context=None)


#filein = open("UKST output.txt", "r")    

#filein = open("POSS-II output.txt", "r")      #sf,sj
#filein = open("POSS-I output.txt", "r")     
#filein = open("UKST output bright all.txt", "r")       


filein = open("UKST MPC.txt", "r")    #0-10

#filein = open("POSS-II MPC.txt", "r")      #0-25
#filein = open("POSS-I MPC.txt", "r")     #0-31
#filein = open("ESO MPC.txt", "r")       #0-1

name =[]
waveb = []
plate_no = []           

field_no = [] 
coords=[]          

for line in filein.readlines():                      #iterates through the lines in the file
    #print(line)
    tokens = line.split(',')                   # breaks line into tokens seperated by ,
                
    name.append(tokens[0])           #appends data in their respective lists
    waveb.append(tokens[1])
    
    plate_no.append(tokens[2])            
    field_no.append(tokens[3])                     #appends array to position list
    coords.append(tokens[4])         
#coordinates =  "18 46 34.21 0 10 23.63"
#wb = 'R'        


#driver.get("http://www-wfau.roe.ac.uk/sss/")

driver.get(url)

# in loop i need to put in J, R, I and coordinates 

for i in range(0,11):

    
    #UKST
    if waveb[i] == 'OR ':
        wb = 'R'
        name_wb = 'OR'
    
    elif waveb[i] == ' J ':
        wb = 'J'
        name_wb = 'J'
        
    elif waveb[i] == ' I ':
        wb = 'I'
        name_wb = 'I'
     
        
    #POSS-II
    elif waveb[i] == 'XJ ':
        wb = 'B'
        name_wb = 'J'
    
    elif waveb[i] == 'XF ':
        wb = 'F'
        name_wb = 'F'
    
    elif waveb[i] == 'XN ':
        wb = 'N'
        name_wb = 'N'
     
    
    #POSS-I
    elif waveb[i] == 'AE ':
        wb = 'A'
        name_wb = 'A'
    
    elif waveb[i] == 'BE ':
        wb = 'P'
        name_wb = 'P'
        
    #ESO    
    elif waveb[i] == ' R ':
        wb = 'E'
        name_wb = 'E'
    coordinates = coords[i]
        
    print(str(name[i])+str(name_wb)+str(plate_no[i]))
    #assert "SuperCOSMOS Sky Surveys" in driver.title
    elem = driver.find_element_by_name("coords")
    elem.clear()
    elem.send_keys(coordinates)
    
    
    select_eq = Select(driver.find_element_by_name('equinox'))
    select_eq.select_by_value('1')
    
    elem_size = driver.find_element_by_name("size")
    elem_size.clear()
    elem_size.send_keys('15')
    
    select_wb = Select(driver.find_element_by_name('waveband'))
    select_wb.select_by_value(wb)
    
    elem_size.send_keys(Keys.RETURN)
    
    time.sleep(6)
    
    driver.find_element_by_link_text('here').click()
    
    driver.find_element_by_link_text('Another extraction?').click()
    

driver.close()







