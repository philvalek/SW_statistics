# -*- coding: utf-8 -*-
"""
Created on Tue Jul 17 12:05:14 2018

@author: Phil

The purpose of this code is to open a webbrowser, navigate to the NASA
CDAWeb (Coordinated Data Analysis Web) site, and read in the OMNI Solar Wind
data.Once the data is read in, this code will perform some analysis to determine
what the average conditions were during a given time period. To make this code
a bit more spiffy, it will take the time range either from the command line or
a default range
"""


"""I honestly probably could have knocked this task out in under an 
hour if I downloaded the data files by hand and used IDL, but I really want to
better understand python. I am going to use In this code I am going to use 
the Selenium package to navigate to the CDAweb site, and use the matplotlib
package for the analysis, plus anything else I find useful
"""

#import the packages
import requests, csv, time, datetime
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC


def NavigateToInput():
    #navigate to the correct page
    #use Firefox as the browser
    browser = webdriver.Firefox()
    
    #goto the CDAWeb homepage and select the OMNI data

    #-->goto the home page and select the OMNI data set
    browser.get('https://cdaweb.sci.gsfc.nasa.gov/index.html/')
    omni_box = browser.find_element_by_id('OMNI__Combined_1AU_IP_Data;_Magnetic_and_Solar_Indices_')
    omni_box.click()
        
    omni_box.submit()

#-->goto the next page and select the data time resolution
#put in a wait to give the page time to load
    element = WebDriverWait(browser, 15).until(
              EC.presence_of_element_located((By.ID, "OMNI_HRO_1MIN"))
              )

#only want 1min, so deselect the other three options
    omni_hro_no = browser.find_element_by_id('OMNI_HRO_5MIN')
    omni_hro_no.click() 
   
    omni_hro_no = browser.find_element_by_id('OMNI2_H0_MRG1HR')
    omni_hro_no.click() 
    
    omni_hro_no = browser.find_element_by_id('OMNI_COHO1HR_MERGED_MAG_PLASMA')
    omni_hro_no.click() 
    
    omni_hro_no.submit()

#-->goto the page where we select the data elements
#now load in the page where we pick time range and data elements to get
#put in a wait to give the page time to load
    element = WebDriverWait(browser, 15).until(
              EC.presence_of_element_located((By.ID, "flow_speed"))
              ) 
    return browser


def SelectValues(browser):
#The default for this page is to plot the results, to get the data set 
#list_data = 1 (i.e. true), other for plots

#The routine ReadData() expects the formating from this call. If the variables
#are changed here, then ReadData needs to be updated.

    list_data = 1
    if list_data:
        dataElem = browser.find_element_by_id('list')
        dataElem.click()
    

#now pick the values to check
         #speed
    dataElem = browser.find_element_by_id('flow_speed')
    dataElem.click()
    
        #IMF
    dataElem = browser.find_element_by_id('BX_GSE')
    dataElem.click()
    dataElem = browser.find_element_by_id('BY_GSE')
    dataElem.click()
    dataElem = browser.find_element_by_id('BZ_GSE')
    dataElem.click()
    
        #velocity
    dataElem = browser.find_element_by_id('Vx')
    dataElem.click()
    dataElem = browser.find_element_by_id('Vy')
    dataElem.click()
    dataElem = browser.find_element_by_id('Vz')
    dataElem.click()

        #density    
    dataElem = browser.find_element_by_id('proton_density')
    dataElem.click()
    
        #dynamic pressure    
    dataElem = browser.find_element_by_id('Pressure')
    dataElem.click()


    
    
def ReadData(curURL):

#This routine will read through the lines of the CDAWeb text file. This code is 
#set up to read a file that is requesting an OMNI data file that has 
#flowspeed, Vx, Vy, Vz, density and dynamic pressure. If the call in SelectValues()
#is changed, then this routine will also need to be updated.

    swdata = requests.get(curURL)
#use the 
        
#find the line where the data headings start
    strt  = swdata.text.find('dd-mm-yyyy')
    stend = swdata.text.find('#', strt)
    
    lineStart = swdata.text.find('nPa',strt) + 3
#split up the text to pull out the data
    noHeader = swdata.text[lineStart:stend].split()
        
#now pull the sub sets and place into arrays, and also convert to float  
#and put into a nparray
    bx_gse    = np.asarray([float(y) for y in noHeader[2::11]])
    by_gse    = np.asarray([float(y) for y in noHeader[3::11]])
    bz_gse    = np.asarray([float(y) for y in noHeader[4::11]])
    flowSpeed = np.asarray([float(y) for y in noHeader[5::11]])
    vx        = np.asarray([float(y) for y in noHeader[6::11]])
    vy        = np.asarray([float(y) for y in noHeader[7::11]])
    vz        = np.asarray([float(y) for y in noHeader[8::11]])
    den       = np.asarray([float(y) for y in noHeader[9::11]])
    dyP       = np.asarray([float(y) for y in noHeader[10::11]])

    
#now we need to remove the bad elements. if the flow speed is 99999.9, then it
#is a bad meausrement and will be rejected, the IMF (B) measurements are from
#a different instrument, so they may not have good values at the same time
#we will only include values where both the mag and particle values are good
    gd = (flowSpeed < 99999) & (bx_gse < 9999)
    bx_gse    = np.extract(gd, bx_gse)
    by_gse    = np.extract(gd, by_gse)
    bz_gse    = np.extract(gd, bz_gse)
    flowSpeed = np.extract(gd, flowSpeed)
    vx        = np.extract(gd, vx)
    vy        = np.extract(gd, vy)
    vz        = np.extract(gd, vz)
    den       = np.extract(gd, den)
    dyP       = np.extract(gd, dyP)

    
    return (bx_gse, by_gse, bz_gse, flowSpeed,vx, vy, vz, den, dyP)
    
    
    
           
def enterStartStop(browser, month, year):
#this routine will take a start and stop month and year, and submit the request    
   if month <= 11:    
      start_time = str(year) + '/' + str(month).zfill(2) + '/01 00:00:00'
      stop_time  = str(year) + '/' + str(month + 1).zfill(2) + '/01 00:00:00'
   else:
      start_time = str(year) + '/' + str(month).zfill(2) + '/01 00:00:00'
      stop_time  = str(year + 1) + '/01/01 00:00:00'
            
   print(start_time + ' to ' + stop_time)
#enter start and stop time
   dataElem = browser.find_element_by_id('start')
   dataElem.clear()
   dataElem.send_keys(start_time)

   dataElem = browser.find_element_by_id('stop')
   dataElem.clear()
   dataElem.send_keys(stop_time)
        
#submit the form
   dataElem.submit()

#switch to the next window
   browser.switch_to.window(browser.window_handles[1])

#now read in the data and start working with it
#   element = WebDriverWait(browser, 30).until(
#             EC.presence_of_element_located((By.LINK_TEXT, 'OMNI_HRO_1MIN_'))
#             )
   element = WebDriverWait(browser, 30).until(
             EC.presence_of_element_located((By.PARTIAL_LINK_TEXT, 'OMNI_HRO_1MIN_'))
             )

#navigate to the data text file, read in the URL and download it using requests
   linkElem = browser.find_element_by_partial_link_text('OMNI_HRO_1MIN_')
   linkElem.click()   
       
   #put a 1 second wait here, I had a few times where the current URL was returned
   #as "about:blank", what I think means a new tab was opened but the page wasn't loaded
   element = WebDriverWait(browser,2)
   
   browser.switch_to.window(browser.window_handles[2])
   curURL = browser.current_url

   return curURL
    
def closeWindows(browser):
#this routine will close the two windows opened from the calls in enterStartStop()
    
    browser.switch_to.window(browser.window_handles[2])
    browser.close()
    
    browser.switch_to.window(browser.window_handles[1])
    browser.close()
    
    browser.switch_to.window(browser.window_handles[0])
    
   
#!!!!!!!!!!!!!!!!!here starts the main code
    
#save_file_name is the file where thesw data will be stored
save_file_name = 'C:/py/SW_statistics/swdata.csv'

st_year = 2013
sp_year = 2014 #end year inclusive

reRead = bool(0) #set to non-zero to re-read in the data from the web, 0 otherwise
if reRead:
    #navigate to the page where values are selected     
    browser = NavigateToInput()   

    #select the variabels
    SelectValues(browser)
    
    #loop through year and month. When we submit the form, the site will open a 
    #new window. The code will then jump back to the first window and update the 
    #start and stop times and repeat
    
    firstFlag = bool(1)
    for year in range(st_year, sp_year+1):
        for month in range(1, 13):
#loop over data range, one month at a time, CDAWeb doesn't like requests larger than this

            curURL = enterStartStop(browser, month, year)
            print(curURL)
        
                 
#read in the data, and reject any bad measurements
            tbx_gse, tby_gse, tbz_gse, tflowSpeed,tvx, tvy, tvz, tden, tdyP = ReadData(curURL)
            
            if firstFlag:
                bx_gse    = tbx_gse
                by_gse    = tby_gse
                bz_gse    = tbz_gse
                flowSpeed = tflowSpeed
                vx        = tvx
                vy        = tvy
                vz        = tvz
                den       = tden
                dyP       = tdyP
                firstFlag = bool(0)
            else:
                bx_gse    = np.concatenate((bx_gse   , tbx_gse))
                by_gse    = np.concatenate((by_gse   , tby_gse))
                bz_gse    = np.concatenate((bz_gse   , tbz_gse))
                flowSpeed = np.concatenate((flowSpeed, tflowSpeed))
                vx        = np.concatenate((vx       , tvx))
                vy        = np.concatenate((vy       , tvy))
                vz        = np.concatenate((vz       , tvz))
                den       = np.concatenate((den      , tden))
                dyP       = np.concatenate((dyP      , tdyP))
            
                        
#close daughter windows
            closeWindows(browser)
    #when done reading in data, close the last window
    browser.close()
#save data to a file so we won't have to read it in again
#write a csv file
    with open(save_file_name,'w') as save_file:
        print('write data to ' + save_file_name)
        header = '#SW data from Omni from' + str(st_year) 
        header += ' to ' + str(sp_year) + "\n"
        header += '#bx_gse, by_gse, bz_gse, flowSpeed, vx, vy, vz, den, dyP\n'
        save_file.write(header)
        
        for i in range(len(vx)):
            out_string  = " "  + str(bx_gse[i])            
            out_string += ", " + str(by_gse[i])            
            out_string += ", " + str(bz_gse[i])            
            out_string += ", " + str(flowSpeed[i])
            out_string += ", " + str(vx[i])
            out_string += ", " + str(vy[i])
            out_string += ", " + str(vz[i])
            out_string += ", " + str(den[i])
            out_string += ", " + str(dyP[i])
            out_string += '\n'
            save_file.write(out_string)
            
            
            
            
else: 
    #the code commented below work, but is slow, I am now going to try
    #to use the csv functions
    firstFlag = bool(1)
    print('read in data from ' + save_file_name)
    cnt = 0
    with open(save_file_name, newline='') as csvfile:
        csvreader = csv.reader(csvfile, delimiter = ',' )
        t0 = time.time()    
#appending values takes for ever. predefine the arrays to see the speed effect
#as an asside, reading in ~ 1e6 lines by appending the array would have taken about
#and hour and a half, pre-setting large arrays and then triming takes a minute(ish)
        elements = 2000000  #this is a number the should be larger than the file
        bx_gse    = np.empty(elements)
        by_gse    = np.empty(elements)
        bz_gse    = np.empty(elements)
        flowSpeed = np.empty(elements)
        vx        = np.empty(elements)
        vy        = np.empty(elements)
        vz        = np.empty(elements)
        den       = np.empty(elements)
        dyP       = np.empty(elements)
        
        
        for row in csvreader:
            if((row[0][0]) == "#"):continue
            
            bx_gse[cnt]    = float(row[0])
            by_gse[cnt]    = float(row[1])
            bz_gse[cnt]    = float(row[2])
            flowSpeed[cnt] = float(row[3])
            vx[cnt]        = float(row[4])
            vy[cnt]        = float(row[5])
            vz[cnt]        = float(row[6])
            den[cnt]       = float(row[7])
            dyP[cnt]       = float(row[8])
            cnt += 1    
               

            if (cnt%100000)==0:
                print(cnt)        
         #now trim the tail off the arrays
        bx_gse    = np.resize(bx_gse, cnt)
        by_gse    = np.resize(by_gse, cnt)
        bz_gse    = np.resize(bz_gse, cnt)
        flowSpeed = np.resize(flowSpeed, cnt)
        vx        = np.resize(vx, cnt)
        vy        = np.resize(vy, cnt)
        vz        = np.resize(vz, cnt)
        den       = np.resize(den, cnt)
        dyP       = np.resize(dyP, cnt)
            
        
        



#calculate mean, standard deviation, and plot a histogram for now.


#make a grid of plots including a histogram of flow speed, den and dyP, and a 
#dial plot which shows flow direction
pdfFileName = 'C:/py/SW_statistics/swdata.pdf'
with PdfPages(pdfFileName) as pdf:
                
    fig = plt.figure(figsize=(7.5, 10))
    plt.subplot(4,1,1)
    title_str = 'Data from ' + str(st_year) + ' to ' + str(sp_year)
    plt.title(title_str)
    
#    plt.tight_layout()
#plot histogram of the flowspeed
    plt.xlabel('FlowSpeed (km/s)')
    plt.ylabel('Minutes')
    plt.hist(flowSpeed, range = (200, 1000), bins = 100, log=1)
    
#plot histogram of the density
    plt.subplot(4,1,2)
    plt.hist(den, range = (0, 40), bins = 100, log=1)
    plt.xlabel('Density (cc)')
    plt.ylabel('Minutes')

#plot histogram of dynamic pressure
    plt.subplot(4,2,5)
    plt.hist(dyP, range = (0, 10), bins = 100, log=1)
    plt.xlabel('Pdyn (nPa)')
    plt.ylabel('Minutes')

    
#make a plot of dynamic pressure as a function of 
    plt.subplot(4,2,6)
#    vd_hist = np.histogram2d(flowSpeed, den, 
#                             bins = [100, 100], range = [[0, 1000],[0, 20]])
#    im_plot = plt.imshow(vd_hist[0], cmap = 'nipy_spectral', 
#                         extent = (min(vd_hist[1]), max(vd_hist[1]), 
#                                   max(vd_hist[2]), min(vd_hist[2])), 
#                         aspect='auto')
#    vd_hist = np.histogram2d(flowSpeed, den,
#                             bins = [100, 100], range = [[0, 500],[0, 40]])
    vd_hist = np.histogram2d(den, flowSpeed,
                             bins = [100, 100])#, range = [[0, 500],[0, 40]])
#NOTE: the histrogram does not follow the Cartesian convention where x
#values are on the abscissa and y valyues on the ordinate axis. Rather x is 
#histogrammed along the first dimension of the array (verical), this ensures 
#compatibility with histogramdd.
    plt.imshow(np.log10(vd_hist[0])    , cmap = "nipy_spectral", 
               extent = (min(vd_hist[1]), max(vd_hist[1]), 
                         min(vd_hist[2]), max(vd_hist[2])),
               origin = 'lower', aspect = 'auto') 

    #make some contours of pressure (isobars)
    #dyP = 2e-6*n*v^2
    vels = np.arange(1,1000)
    for p in range(1,11):
        n = p/(2e-6 * vels**2)
        plt.plot(vels, n)
        plt.axis([min(vd_hist[1]), max(vd_hist[1]), 
                  min(vd_hist[2]), max(vd_hist[2])])
    plt.colorbar().set_label('Minutes')
    plt.xlabel('FlowSpeed(km/s)')
    plt.ylabel('Density(cc)')
     
    theta = np.arctan2(vz, vy)
    phi   = np.degrees(np.arctan2(np.sqrt(vz**2 + vy**2), -vx))
    
    angx = phi*np.cos(theta)
    angy = phi*np.sin(theta)
    plt.subplot(4,2,7)
    rng = 15
    vel_hist = np.histogram2d(angx, angy, bins = 200, 
                              range= ([-rng, rng], [-rng, rng]))
    plt.imshow(np.log10(vel_hist[0])    , cmap = "nipy_spectral", 
               extent = (min(vel_hist[1]), max(vel_hist[1]), 
                         min(vel_hist[2]), max(vel_hist[2])),
               origin = 'lower') 
    plt.colorbar().set_label('log(events)')
    plt.xlabel('Vy(deg)')
    plt.ylabel('Vz(deg)')
    
    
    
    plt.subplot(4,2,8)
    rng = 20
    b_hist = np.histogram2d(by_gse, bz_gse, bins = 200, 
                            range= ([-rng, rng], [-rng, rng]))
    plt.imshow(np.log10(b_hist[0])    , cmap = "nipy_spectral", 
               extent = (min(b_hist[1]), max(b_hist[1]), 
                         min(b_hist[2]), max(b_hist[2])),
               origin = 'lower')
    plt.colorbar().set_label('log(events)')
    plt.xlabel('By_gse')
    plt.ylabel('Bz_gse')
    plt.tight_layout()   
    pdf.savefig(orientation='portrait')
    plt.tight_layout().set_label('log(events)')
    plt.show()
    plt.close()
    
    # Set the file's metadata via the PdfPages object:
    d = pdf.infodict()
    d['Title'] = 'SW statistics'
    d['Author'] = u'Philip Valek'
    d['CreationDate'] = datetime.datetime(2018, 8, 7)
    d['ModDate'] = datetime.datetime.today()
 
                    
                    
                    