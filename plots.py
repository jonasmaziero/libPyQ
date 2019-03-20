import matplotlib.pyplot as plt
import numpy as np


def error_bar():
    x = [1,2,3,4,5]
    y = [1,2,3,4,5]
    yerr = [0.1,0.2,0.2,0.3,0.1]
    plt.errorbar(x,y,yerr)
    #plt.plot(x,y)
    plt.show()


def plot2d(x, y):
    # x = [1,2,3]
    # y1 = [1,2,3]
    # y2 = [2,3,4]
    plt.plot(x, y, label='')
    # plt.plot(x,y2,label='y2')
    plt.xlabel('x')
    plt.ylabel('y')
    # plt.title('title\nsubtitle')
    plt.legend()
    plt.show()


def plot2ds(x, y1, y2, y3, y4, y5, y6):
    plt.plot(x, y1, label='di')
    plt.plot(x, y2, label='cc')
    plt.plot(x, y3, label='im')
    plt.plot(x, y4, label='diE')
    plt.plot(x, y5, label='ccE')
    plt.plot(x, y6, label='imE')
    plt.xlabel('p')
    plt.legend()
    plt.show()


def plotBar():
    x = [2, 4, 6, 8, 10]
    y = [6, 7, 8, 2, 4]
    x2 = [1, 3, 5, 7, 9]
    y2 = [7, 8, 2, 9, 1]
    plt.bar(x, y, label='barPlot1', color='blue')
    plt.bar(x2, y2, label='barPlot2', color='red')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.legend()
    plt.show()


def plotHistogram():
    populationAges = [22, 55, 62, 45, 21, 22, 34, 42, 42, 4, 99, 102, 55,
                      44, 66, 77, 33, 22, 99, 88, 77, 66, 55, 44, 33, 11, 23,
                      45, 67, 89]
#    ids = [x for x in range(len(populationAges))]
#    plt.bar(ids,populationAges)
    bins = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110]
    plt.hist(populationAges, bins, histtype='bar', rwidth=0.8)
    plt.xlabel('x')
    plt.ylabel('y')
    plt.legend()
    plt.show()


def plotScatter(x, y):
    # x = [1,2,3,4,5,6,7,8]
    # y = [3,5,6,7,5,4,3,4]
    plt.scatter(x, y, label='scatter', color='blue', s=5, marker='*')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.legend()
    plt.show()


def plotStack():
    days = [1, 2, 3, 4, 5]
    sleeping = [7, 6, 8, 11, 7]
    eating = [2, 3, 4, 3, 2]
    working = [7, 8, 7, 2, 2]
    playing = [8, 5, 7, 8, 13]
    plt.plot([], [], color='magenta', label='sleeping', linewidth=5)
    plt.plot([], [], color='cyan', label='eating', linewidth=5)
    plt.plot([], [], color='gray', label='working', linewidth=5)
    plt.plot([], [], color='green', label='playing', linewidth=5)
    plt.stackplot(days, sleeping, eating, working, playing,
                  colors=['magenta', 'cyan', 'gray', 'green'])
    plt.xlabel('x')
    plt.ylabel('y')
    plt.legend()
    plt.show()


def plotPie():
    days = [1, 2, 3, 4, 5]
    sleeping = [7, 6, 8, 11, 7]
    eating = [2, 3, 4, 3, 2]
    working = [7, 8, 7, 2, 2]
    playing = [8, 5, 7, 8, 13]
    slices = [7, 2, 2, 13]
    activities = ['sleeping', 'eating', 'working', 'playing']
    cols = ['magenta', 'cyan', 'gray', 'green']
    plt.pie(slices, labels=activities, colors=cols, startangle=90,
            shadow=True, explode=[0, 0.2, 0.1, 0], autopct='%1.1f%%')
    plt.title('pie chart')
    plt.show()


def plotFiles():
    '''
    # 1st way
    import csv
    x = []
    y = []
    with open('data.txt','r') as csvfile:
        plots = csv.reader(csvfile, delimiter=',')
        for row in plots:
            x.append(int(row[0]))
            y.append(int(row[1]))
    plt.plot(x,y, label='datafile')
    '''
    # 2nd way
    x, y = np.loadtxt('data.txt', delimiter=',', unpack=True)
    plt.plot(x, y, label='datafile')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.legend()
    plt.show()


import urllib
import matplotlib.dates as mdates


def bytespdate2num(fmt, encoding='utf-8'):
    strconverter = mdates.strpdate2num(fmt)

    def bytesconverter(b):
        s = b.decode(encoding)
        return strconverter(s)
    return bytesconverter


def plotInternet(stock):
    stock_price_url = 'http://chartapi.finance.yahoo.com/instrument/1.0/' + \
        stock+'/chartdata;type=quote;range=10y/csv'
    source_code = urllib.request.urlopen(stock_price_url).read().decode()
    stock_data = []
    split_source = source_code.split('\n')
    for line in split_source:
        split_line = line.split(',')
        if len(split_line) == 6:
            if 'values' not in line and 'labels' not in line:
                stock_data.append(line)
    date, closep, highp, lowp, openp, volume = np.loadtxt(stock_data,
                                                          delimiter=',',
                                                          unpack=True,
                                                          # %Y = full year. 2015
                                                          # %y = partial year. 15
                                                          # %m = number month
                                                          # %D = number day
                                                          # %H = hours
                                                          # %M = minutes
                                                          # %S = seconds
                                                          # 12-26-2014 -> %m-%D-%Y
                                                          converters={0: bytespdates2num('%Y%m%D')})
    plt.plot_date(date, closep, '-', label='prices from web')
    plt.xlabel('date')
    plt.ylabel('price')
    plt.legend()
    plt.show()


def plotCustomize():
    x = [1, 2, 3, 4, 5, 6, 7, 8]
    y = [3, 5, 6, 7, 5, 4, 3, 4]

    fig = plt.figure()
    # The 1st tuple is for how many plots. The 2nd is for the origin.
    ax1 = plt.subplot2grid((1, 1), (0, 0))

    ax1.scatter(x, y, label='scatter', color='blue', s=300, marker='*')
    for label in ax1.xaxis.get_ticklabels():
        label.set_rotation(45)  # this rotates the label
    ax1.grid(True, color='g', linestyle='-', linewidth=0.5)

    plt.xlabel('x')
    plt.ylabel('y')
    plt.legend()
    plt.subplots_adjust(left=0.09, bottom=0.18, right=0.94, top=0.80, wspace=0.2, hspace=0)
    plt.show()


def level_curves():
    import scipy.interpolate
    N = 1000  # number of points for plotting/interpolation
    x, y, z = np.genfromtxt(
        r'/home/jonasmaziero/Dropbox/Research/qnesses/interplay/dipolar/dipolarCalc/EErhotpi2.dat',
        unpack=True)
    xll = x.min()
    xul = x.max()
    yll = y.min()
    yul = y.max()
    xi = np.linspace(xll, xul, N)
    yi = np.linspace(yll, yul, N)
    zi = scipy.interpolate.griddata(
        (x, y), z, (xi[None, :], yi[:, None]), method='linear', rescale=False)
    contours = plt.contour(xi, yi, zi, 6, colors='black')
    plt.clabel(contours, inline=True, fontsize=7)
    plt.imshow(zi, extent=[xll, xul, yll, yul], origin='lower', cmap=plt.cm.jet, alpha=0.9)
    plt.xlabel(r'$w$')
    plt.ylabel(r'$p$')
    plt.clim(0, 1)
    # plt.colorbar()
    plt.savefig('/home/jonasmaziero/Dropbox/Research/qnesses/interplay/dipolar/dipolarCalc/EErhotpi2.eps',
                format='eps', dpi=100)
    plt.show()
