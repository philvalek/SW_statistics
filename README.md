# SW_statistics
This python code uses the selenium package to read in data from NASA's CDAWeb and make plots of the solar wind.

Data from downloaded from the NASA CDAWeb site (https://cdaweb.gsfc.nasa.gov/index.html/) is used to make plots of the solar wind properties for the years 2013 and 2014. This is used in the design of a plasma instrument so we know the measurement range we need to have.

The code first uses the selenium package to navigate CDAWeb and download the solar wind data from OMNI. (OMNI is not an instrument, but rather a group that propagates the data taken at L1 to Earth.) I wouldn't be surprised if there are tools out there to read data from either CDAWeb or OMNI specifically into python, but I didn't look too hard. Anyways, this was a good exercise to write some code to do web scraping.

After the data is read in, I do some very simple histograms to look at the data.
