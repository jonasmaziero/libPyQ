import urllib2
#from urllib.request import urlopen

def getRandomBin():
    url = 'http://150.203.48.55/RawBin.php'
    page = urlopen(url, timeout=5)
    data = page.read()
    num = data.split('"rng"')[1].split('<td>\n')[1].split('</td>')[0]
    return num
