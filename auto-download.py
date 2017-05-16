#Ya know. For when you ruin your data, and just need to redownload all ~100 tiles.
#Cause I'm dumb!
import urllib.request
import numpy as np

good_ims = ["014","015","016","017","018","019","020","021","026","027","028","029","030","031","032","033",
           "038","039","040","041","042","043","044","045","050","051","052","053","054","055","056","057",
           "062","063","064","065","066","067","068","069","074","075","076","077","078","079","080","081",
           "086","087","088","089","090","091","092","093","098","099","100","101","102","103","104","105",
           "110","111","112","113","114","115","116","117"]

url = "http://irsa.ipac.caltech.edu/data/COSMOS/images/subaru/original_psf/"

for i in good_ims:
    urllib.request.urlretrieve(url + "B/subaru_B_"+i+"_sci_20.fits","subaru_B_"+i+"_sci_20.fits")
    urllib.request.urlretrieve(url + "IB427/subaru_IB427_"+i+"_sci_10.fits","subaru_IB427_"+i+"_sci_10.fits")
    

