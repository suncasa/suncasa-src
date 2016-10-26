#= Smoothing of a 1D signal =

#This method is based on the convolution of a scaled window with the signal.
#The signal is prepared by introducing reflected window-length copies 
#of the signal at both ends so that boundary effect are minimized
#in the beginning and end part of the output signal.

#== Code ==

#{{{#!python
import numpy

def smooth(x,xshift=4,window_len=11,window='hanning'):
    """smooth the data using a window with requested size.
    
    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal 
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.
    
    input:
        x: the input signal 
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal
        
    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)
    
    see also: 
    
    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter
 
    TODO: the window parameter could be the window itself if an array instead of a string
    NOTE: length(output) != length(input), to correct this: return y[(window_len/2-1):-(window_len/2)] instead of just y.
    """ 
     
    if numpy.ndim(x) != 1:
        raise ValueError, "smooth only accepts 1 dimension arrays."

    if numpy.size(x) < window_len:
        raise ValueError, "Input vector needs to be bigger than window size."
        

    if window_len<3:
        return x
    
    
    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"

    if not (xshift and (type(xshift) == int)):
        xshift=int(4)
        
    s=numpy.r_[x[window_len-1:0:-1],x,x[-1:-window_len:-1]]
    #print(len(s))
    if window == 'flat': #moving average
        w=numpy.ones(window_len,'d')
    else:
        w=eval('numpy.'+window+'(window_len)')

    ylen=len(x)
    #print 'number of original data points: ',ylen
    y=numpy.convolve(w/w.sum(),s,mode='valid')
    #print 'number of convolved data points: ',len(y)
    temp=numpy.zeros(ylen)
    if xshift > len(y)-ylen:
        temp[0:len(y)-xshift]=y[xshift:len(y)]
        #print len(y)-xshift
        #print ylen
        #temp[len(y)-xshift+1:ylen]=0
        #temp[len(y)-xshift+1:ylen]=x[len(y)-xshift+1:ylen]
        temp[len(y)-xshift+1:ylen]=y[len(y)-xshift+1:ylen]
        y=temp
    else:
        y=y[xshift:xshift+ylen]
    # print 'number of shifted data points',len(y)
    return y

    
    
    
from numpy import *
from pylab import *
import task_smooth
def smooth_demo():    
    
    t=linspace(-4,4,100)
    x=sin(t)
    xn=x+randn(len(t))*0.1
    y=task_smooth.smooth(x)
    
    ws=31
    
    subplot(211)
    plot(ones(ws))
    
    windows=['flat', 'hanning', 'hamming', 'bartlett', 'blackman']
    
    hold(True)
    for w in windows[1:]:
        eval('plot('+w+'(ws) )')
    
    axis([0,30,0,1.1])
    
    legend(windows)
    title("The smoothing windows")
    subplot(212)
    plot(x)
    plot(xn)
    for w in windows:
        plot(task_smooth.smooth(xn,10,w))
    l=['original signal', 'signal with noise']
    l.extend(windows)
    
    #legend(l)
    title("Smoothing a noisy signal")
    show()

    
if __name__=='__main__':
    smooth_demo()
#}}}

#== Figure ==

#inline:smoothsignal.jpg

#= Smoothing of a 2D signal =

#Convolving a noisy image with a gaussian kernel (or any bell-shaped
#curve) blurs the noise out and leaves the low-frequency details of the
#image standing out.

#== Functions used ==

#{{{#!python
def gauss_kern(size, sizey=None):
    """ Returns a normalized 2D gauss kernel array for convolutions """
    size = int(size)
    if not sizey:
        sizey = size
    else:
        sizey = int(sizey)
    x, y = mgrid[-size:size+1, -sizey:sizey+1]
    g = exp(-(x**2/float(size)+y**2/float(sizey)))
    return g / g.sum()

def blur_image(im, n, ny=None) :
    """ blurs the image by convolving with a gaussian kernel of typical
        size n. The optional keyword argument ny allows for a different
        size in the y direction.
    """
    g = gauss_kern(n, sizey=ny)
    improc = signal.convolve(im,g, mode='valid')
    return(improc)
#}}}

#== Example ==

#{{{#!python numbers=disable
#from scipy import *

#X, Y = mgrid[-70:70, -70:70]
#Z = cos((X**2+Y**2)/200.)+ random.normal(size=X.shape)
#}}}

#inline:noisy.png

#{{{#!python numbers=disable
#blur_image(Z, 3)
#}}}

#inline:convolved.png


#The attachment attachment:cookb_signalsmooth.py contains a version of this script with some stylistic cleanup.

#= See Also =
#["Cookbook/FiltFilt"] which can be used to smooth the data by low-pass filtering and does not delay the signal (as this smoother does).  