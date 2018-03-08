from urllib.request import urlopen
#import json
#import pickle


def anu():
    url = 'https://qrng.anu.edu.au/API/jsonI.php?length=10&type=uint16'
    page = urlopen(url, timeout=10)
    print(page)
    #data = page.read()
    #ty, le, da, st = pickle.load(page, encoding = 'latin1')
    # print(da)
