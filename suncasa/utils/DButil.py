import glob
import json
import os
import pickle
from functools import wraps

import astropy.units as u
import numpy as np

__author__ = ["Sijie Yu"]
__email__ = "sijie.yu@njit.edu"


def img2html_movie(imgprefix, outname='movie', img_fmt='png'):
    imgfiles = glob.glob(imgprefix + '*.' + img_fmt)
    imgfiles = sorted(imgfiles)

    # try:
    import matplotlib.image as mpimg
    img = mpimg.imread(imgfiles[0])
    height, width, dummy = img.shape
    # except:
    #     from PIL import Image
    #     img = Image.open(imgfiles[0])
    #     width, height = img.size
    #     img.close()

    nfiles = len(imgfiles)
    htmlfile = os.path.join(os.path.dirname(imgprefix), '{}.html'.format(outname))
    fi = open(htmlfile, 'w')
    fi.write('<HTML>\n')
    fi.write('\n')
    fi.write('<HEAD>\n')
    fi.write('<TITLE> Javascript Movie Player </TITLE>\n')
    fi.write('\n')
    fi.write('</HEAD>\n')
    fi.write('\n')
    fi.write('<BODY FONT=3 BGCOLOR=FFFFFF LINK=#0000CC VLINK=#0000CC TEXT=#000000> \n')
    fi.write('\n')
    fi.write('<SCRIPT language="Javascript">\n')
    fi.write('<!--\n')
    fi.write('//\n')
    fi.write('// DYNAMICALLY CREATED HTML - DO NOT EDIT\n')
    fi.write('//\n')
    fi.write('// Javascript program to produce fast animations by reading from cache\n')
    fi.write('// Written: 7 July 1997, Zarro, NASA/GSFC\n')
    fi.write('// Modified: 27 July 1997, Freeland, LSAL - added nice SWING option\n')
    fi.write('// Modified: 9 Oct 1998, Zarro, SMA/GSFC - added context window button\n')
    fi.write('// Modified: 14 Oct 1998, Zarro, SMA/GSFC - load images asynchronously\n')
    fi.write('// Modified: 23 Aug 1999, William Thompson, GSFC - make button text black\n')
    fi.write('\n')
    fi.write('// Following variables control movie speed, labels, etc.\n')
    fi.write('\n')
    fi.write('var imax = {}, inc = 1.50, delay = 250;\n'.format(nfiles))
    fi.write('var num_loaded_images = 0;\n')
    fi.write('var frame=-1, speed=10;\n')
    fi.write('var timeout_id=null;\n')
    fi.write('var dir=1, playing=0, swingon=0, run=0;\n')
    fi.write('var bname = "Reverse";\n')
    fi.write('var url_path = ".";\n')
    fi.write('var context = "frame00.gif";\n')
    fi.write('var ctitle ="test";\n')
    fi.write('var url_context=url_path+"/"+context;\n')
    fi.write('var iwidth = {},iheight = {};\n'.format(width, height))
    fi.write('var index=0;\n')
    fi.write('\n')
    fi.write('// -->\n')
    fi.write('</SCRIPT>\n')
    fi.write('\n')
    fi.write('<CENTER>\n')
    fi.write('\n')
    fi.write('<H1> Javascript Movie Player </H1>\n')
    fi.write('<P>\n')
    fi.write('\n')
    fi.write('<TABLE BORDER="10" CELLPADDING="8">\n')
    fi.write('<TR>\n')
    fi.write('<TD align="center"> \n')
    fi.write('<img NAME=animation ALT="FRAME" width={} height={}>\n'.format(width, height))
    fi.write('</TR>\n')
    fi.write('</TABLE>\n')
    fi.write('<P>\n')
    fi.write('\n')
    fi.write('<FORM NAME="form">\n')
    fi.write(' <FONT COLOR="Black">\n')
    fi.write(' <INPUT TYPE=button VALUE="Start" onClick="start_play();">\n')
    fi.write(' <INPUT TYPE=button VALUE="Pause" onClick="stop_play();">\n')
    fi.write(' <INPUT TYPE=button VALUE="Faster" onClick="speed*=inc; show_speed();">\n')
    fi.write(' <INPUT TYPE=button VALUE="Slower" onClick="speed/=inc; show_speed();">\n')
    fi.write(' <INPUT TYPE=button VALUE="Step" onClick="oneStep();">\n')
    fi.write(' <INPUT TYPE=button NAME="direction" VALUE="Reverse" onClick="reverse();">\n')
    fi.write(' <INPUT TYPE=button VALUE="Swing Mode:" onClick="swing_mode();">\n')
    fi.write(' <INPUT TYPE=text VALUE="OFF" NAME="swing" SIZE=3>\n')
    fi.write('\n')
    fi.write(' </FONT>\n')
    fi.write(' <BR>\n')
    fi.write(' Frame: <INPUT TYPE=text VALUE="" NAME="frame" SIZE=22> \n')
    fi.write(' &nbsp; Speed:<INPUT TYPE=text VALUE="" NAME="rate" SIZE=4> (frames/sec)\n')
    fi.write('</FORM>\n')
    fi.write(' \n')
    fi.write('</FORM>\n')
    fi.write('</CENTER>\n')
    fi.write('\n')
    fi.write('<P>\n')
    fi.write('<HR>\n')
    fi.write('<B>Document</B>: <I><SCRIPT>document.write(document.title);</SCRIPT></I><BR>\n')
    fi.write('<B>URL</B>: <I><SCRIPT>document.write(document.URL);</SCRIPT></I><BR>\n')
    fi.write('<B>Last Update</B>: <I><SCRIPT>document.write(document.lastModified);</SCRIPT></I><BR>\n')
    fi.write('\n')
    fi.write('\n')
    fi.write('<SCRIPT LANGUAGE="JavaScript">\n')
    fi.write('<!--\n')
    fi.write('\n')
    fi.write('function get_window_width(fac,def) {  // Return window width\n')
    fi.write('\n')
    fi.write(' \n')
    fi.write('var width;\n')
    fi.write('if (!fac) {fac=.75;}\n')
    fi.write('if (!def) {def=512;}\n')
    fi.write('\n')
    fi.write('if (window.screen) {\n')
    fi.write(' width=parseInt(fac*parseFloat(window.screen.width));\n')
    fi.write('} else {\n')
    fi.write(' width=def;\n')
    fi.write('}\n')
    fi.write('\n')
    fi.write('return width;\n')
    fi.write('\n')
    fi.write('}\n')
    fi.write('\n')
    fi.write('/////////////////////////////////////////////////////////////////////////////\n')
    fi.write('\n')
    fi.write('function get_window_height(fac,def) {   // Return window height\n')
    fi.write('\n')
    fi.write(' \n')
    fi.write('var height;\n')
    fi.write('if (!fac) {fac=.75;}\n')
    fi.write('if (!def) {def=512;}\n')
    fi.write('\n')
    fi.write('if (window.screen) {\n')
    fi.write(' height=parseInt(fac*parseFloat(window.screen.height));\n')
    fi.write('} else {\n')
    fi.write(' height=def;\n')
    fi.write('}\n')
    fi.write('\n')
    fi.write('return height;\n')
    fi.write('\n')
    fi.write('}\n')
    fi.write('\n')
    fi.write('/////////////////////////////////////////////////////////////////////////////\n')
    fi.write('// Javascript Pop-up image Window\n')
    fi.write('// Written: Zarro (SAC/GSFC), October 1998, (dzarro@solar.stanford.edu)\n')
    fi.write('\n')
    fi.write('var win;\n')
    fi.write('\n')
    fi.write('function pop_img(url,title,width,height) {\n')
    fi.write('\n')
    fi.write('if (!url) {\n')
    fi.write(" alert('Image URL not entered');\n")
    fi.write(' return false;\n')
    fi.write('}\n')
    fi.write('\n')
    fi.write('// default to fit within 75% of browser window, or 512x512\n')
    fi.write('\n')
    fi.write('if (!width) {width=get_window_width(.75,512);}\n')
    fi.write('if (!height) {height=get_window_height(.75,512);}\n')
    fi.write('\n')
    fi.write('if (!win || win.closed) {\n')
    fi.write(
        ' win = open("","img","width="+width.toString()+",height="+height.toString()+",scrollbars=yes,resizable=yes"); \n')
    fi.write('}\n')
    fi.write('\n')
    fi.write('// dynamically create HTML, adding a CLOSE button\n')
    fi.write('\n')
    fi.write('if (win) {\n')
    fi.write(' d=win.document;\n')
    fi.write(' d.write(\'<html><head><title>Image Window</title></head><body bgcolor="white"><center>\');\n')
    fi.write(' if (title) {\n')
    fi.write("  d.write('<h1>'+title+'</h1>');\n")
    fi.write(' }\n')
    fi.write(" d.write('<img src='+url+'>');\n")
    fi.write(" d.write('<br><br><br>');\n")
    fi.write(
        ' d.write(\'<form><b><input type="button" value="CLOSE" onClick="self.close();"></b></form></center>\'); \n')
    fi.write(" d.write('</html></body>');\n")
    fi.write(' d.close();\n')
    fi.write(' win.focus();\n')
    fi.write('}\n')
    fi.write('\n')
    fi.write('return true;     \n')
    fi.write('}\n')
    fi.write('\n')
    fi.write('///////////////////////////////////////////////////////////////////////////\n')
    fi.write('\n')
    fi.write('function load_img() {        // asynchronously load all images into cache\n')
    fi.write(' for (i=0; i < imax ; i++) {\n')
    fi.write('  id[i]=setTimeout("load_src()",0);\n')
    fi.write(' }\n')
    fi.write(' return;\n')
    fi.write('}\n')
    fi.write('/////////////////////////////////////////////////////////////////////////////\n')
    fi.write('\n')
    fi.write('function load_src() {      // load individual images into cache\n')
    fi.write('\n')
    fi.write(' if (index < imax) {\n')
    fi.write('  if (iwidth && iheight) {\n')
    fi.write('   images[index] = new Image(iwidth,iheight);\n')
    fi.write('  } else {\n')
    fi.write('   images[index] = new Image();\n')
    fi.write('  }\n')
    fi.write('  images[index].onload=count_images;\n')
    fi.write('  images[index].src = urls[index];\n')
    fi.write('  index++;\n')
    fi.write(' }\n')
    fi.write(' return;\n')
    fi.write('}\n')
    fi.write('\n')
    fi.write('/////////////////////////////////////////////////////////////////////////////\n')
    fi.write('\n')
    fi.write("function clear_ids() {         // clear asynchronous id's\n")
    fi.write(' for (i=0; i < imax ; i++) {clearTimeout(id[i]);}\n')
    fi.write(' return;\n')
    fi.write('}\n')
    fi.write('\n')
    fi.write('\n')
    fi.write('/////////////////////////////////////////////////////////////////////////////\n')
    fi.write('\n')
    fi.write('function count_images() // count images as they are loaded into cache\n')
    fi.write('{ \n')
    fi.write(' if (++num_loaded_images == imax) {\n')
    fi.write('  show_speed();\n')
    fi.write('  clear_ids();\n')
    fi.write('  animate();\n')
    fi.write(' } else {\n')
    fi.write('  document.animation.src=images[num_loaded_images-1].src;\n')
    fi.write('  document.form.frame.value="Loading "+num_loaded_images+" of "+imax; \n')
    fi.write(' }\n')
    fi.write('}\n')
    fi.write('\n')
    fi.write('///////////////////////////////////////////////////////////////////////////\n')
    fi.write('\n')
    fi.write('function image_abort() //  abort loading images\n')
    fi.write('{ \n')
    fi.write(' imax=num_loaded_images;\n')
    fi.write(' if (!images[num_loaded_images].complete) imax=imax-1;\n')
    fi.write(' alert("Aborting");\n')
    fi.write(' if (imax > -1) animate(); \n')
    fi.write('}\n')
    fi.write('\n')
    fi.write('///////////////////////////////////////////////////////////////////////////\n')
    fi.write('\n')
    fi.write('function image_error(message) //  abort loading images\n')
    fi.write('{ \n')
    fi.write(' alert(message);\n')
    fi.write('}\n')
    fi.write('\n')
    fi.write('///////////////////////////////////////////////////////////////////////////\n')
    fi.write('\n')
    fi.write('function start_play()  // start movie\n')
    fi.write('{\n')
    fi.write(' if (playing == 0) {\n')
    fi.write('  if (timeout_id == null && num_loaded_images==imax) animate();\n')
    fi.write(' }\n')
    fi.write('} \n')
    fi.write('\n')
    fi.write('///////////////////////////////////////////////////////////////////////////\n')
    fi.write('\n')
    fi.write('function stop_play() // stop movie\n')
    fi.write('{ \n')
    fi.write(' if (timeout_id) clearTimeout(timeout_id); \n')
    fi.write('  timeout_id=null;\n')
    fi.write('  playing = 0;\n')
    fi.write('} \n')
    fi.write('\n')
    fi.write('///////////////////////////////////////////////////////////////////////////\n')
    fi.write('\n')
    fi.write('function swing_mode()    // set swing mode\n')
    fi.write('{\n')
    fi.write(' if (swingon) {\n')
    fi.write('  swingon=0;\n')
    fi.write('  document.form.swing.value="OFF";\n')
    fi.write(' }\n')
    fi.write('  else {\n')
    fi.write('  swingon=1;\n')
    fi.write('  document.form.swing.value="ON";\n')
    fi.write(' }\n')
    fi.write('}\n')
    fi.write('\n')
    fi.write('///////////////////////////////////////////////////////////////////////////\n')
    fi.write('  \n')
    fi.write('function animate()  // control movie loop\n')
    fi.write('{\n')
    fi.write(' var j;\n')
    fi.write(' frame=(frame+dir+imax)%imax;\n')
    fi.write(' j=frame+1;\n')
    fi.write(' if (images[frame].complete) {\n')
    fi.write('  document.animation.src=images[frame].src;\n')
    fi.write('  document.form.frame.value="Displaying "+j+" of "+imax;\n')
    fi.write('  if (swingon && (j == imax || frame == 0)) reverse();\n')
    fi.write('  timeout_id=setTimeout("animate()",delay);\n')
    fi.write('  playing=1;\n')
    fi.write(' }\n')
    fi.write('}\n')
    fi.write('\n')
    fi.write('///////////////////////////////////////////////////////////////////////////\n')
    fi.write('\n')
    fi.write('function oneStep() // step frames\n')
    fi.write('{\n')
    fi.write(' var j;\n')
    fi.write(' if (timeout_id) clearTimeout(timeout_id); timeout_id=null;\n')
    fi.write(' frame=(frame+dir+imax)%imax;\n')
    fi.write(' j=frame+1;\n')
    fi.write(' if (images[frame].complete) {\n')
    fi.write('  document.animation.src=images[frame].src;\n')
    fi.write('  document.form.frame.value="Displaying "+j+" of "+imax;\n')
    fi.write('  playing=0;\n')
    fi.write(' }\n')
    fi.write('}\n')
    fi.write('\n')
    fi.write('///////////////////////////////////////////////////////////////////////////\n')
    fi.write('\n')
    fi.write('function reverse()  // reverse direction\n')
    fi.write('{\n')
    fi.write(' dir=-dir;\n')
    fi.write(' if (dir > 0) document.form.direction.value="Reverse"; bname="Reverse";\n')
    fi.write(' if (dir < 0) document.form.direction.value="Forward"; bname="Forward";\n')
    fi.write('}\n')
    fi.write('\n')
    fi.write('///////////////////////////////////////////////////////////////////////////\n')
    fi.write('\n')
    fi.write('function show_speed()      // show speed\n')
    fi.write('{\n')
    fi.write('  document.form.rate.value=Math.round(speed);\n')
    fi.write('  delay = 1000.0/speed; \n')
    fi.write('}\n')
    fi.write('\n')
    fi.write('///////////////////////////////////////////////////////////////////////////\n')
    fi.write('// actual image loading is done here\n')
    fi.write('\n')
    fi.write('show_speed();\n')
    fi.write('images = new Array(imax);          \n')
    fi.write('urls= new Array(imax);\n')
    fi.write('id= new Array(imax);\n')
    fi.write('\n')
    fi.write('\n')
    for i in range(nfiles):
        fi.write('urls[{}]=url_path+"/{}";\n'.format(i, os.path.basename(imgfiles[i])))
    fi.write('\n')
    fi.write('load_img();\n')
    fi.write('\n')
    fi.write('// -->\n')
    fi.write('</SCRIPT>\n')
    fi.write('\n')
    fi.write('</BODY>\n')
    fi.write('</HTML>\n')
    fi.write('\n')
    fi.close()  # pdb.set_trace()
    print('\nWrite movie to {}'.format(htmlfile))
    # print('To open -----> firefox {} &'.format(os.path.abspath(htmlfile)))


def my_timer(orig_func):
    import time

    @wraps(orig_func)
    def wrapper(*args, **kwargs):
        t1 = time.time()
        result = orig_func(*args, **kwargs)
        t2 = time.time() - t1
        print('{} ran in: {} sec'.format(orig_func.__name__, t2))
        return result

    return wrapper


def spectrogram2wav(spec, threshld=None, gama=1, fs=1.0, t_length=None, w=1, wavfile='output.wav'):
    '''
    Convert spectrogram to audio in WAV format
    :param spec: spec.shape (nfreq, ntime)
    :param threshld: below which is set to be threshold
    :param gama:
    :param fs:
    :param t_length: time duration of output WAV file
    :param w: width of the smooth window, if apply
    :param wavfile:
    :return:
    '''
    from scipy import signal
    from scipy.io.wavfile import write
    import numpy as np

    if np.sum(np.isnan(spec)) > 0:
        for idx in range(spec.shape[1]):
            p_slice = spec[:, idx]
            mask_nan = np.isnan(p_slice)
            if np.sum(~mask_nan) > 0:
                p_slice[mask_nan] = np.interp(np.flatnonzero(mask_nan), np.flatnonzero(~mask_nan), p_slice[~mask_nan])
                spec[:, idx] = p_slice
        for idx in range(spec.shape[0]):
            p_slice = spec[idx, :]
            mask_nan = np.isnan(p_slice)
            if np.sum(~mask_nan) > 0:
                p_slice[mask_nan] = np.interp(np.flatnonzero(mask_nan), np.flatnonzero(~mask_nan), p_slice[~mask_nan])
                spec[idx, :] = p_slice

    # smooth the dynamic spectrum
    if w > 1:
        k = 1 - np.abs(np.linspace(-1, 1, w))
        kernel = k.reshape(w, 1) * k.reshape(1, w)
        kernel /= kernel.sum()  # kernel should sum to 1!  :)
        Zxx = signal.convolve2d(spec, kernel, mode='same')
    else:
        Zxx = spec
    if threshld:
        Zxx = np.where(Zxx >= threshld, Zxx, threshld)
    _, xrec = signal.istft(Zxx ** gama, fs)

    scaled = np.int16(xrec / np.max(np.abs(xrec)) * 32767)
    if t_length:
        sample_rate = int(np.round(len(scaled) / t_length))
    else:
        sample_rate = 44100
    write(wavfile, sample_rate, scaled)
    return


def smooth(x, window_len=11, window='hanning'):
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


    NOTE: length(output) != length(input), to correct this: return y[(window_len/2-1):-(window_len/2)] instead of just y.
    """

    if x.ndim != 1:
        raise ValueError("smooth only accepts 1 dimension arrays.")

    if x.size < window_len:
        raise ValueError("Input vector needs to be bigger than window size.")

    if window_len < 3:
        return x

    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError("Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'")

    s = np.r_[x[window_len - 1:0:-1], x, x[-2:-window_len - 1:-1]]
    # print(len(s))
    if window == 'flat':  # moving average
        w = np.ones(window_len, 'd')
    else:
        w = eval('np.' + window + '(window_len)')

    y = np.convolve(w / w.sum(), s, mode='same')
    if len(x) == len(y):
        return y
    else:
        return y[window_len - 1:-(window_len - 1)]


# def smooth(x, window_len=11, window='hanning', mode='same'):
#     """smooth the data using a window with requested size.
#
#     This method is based on the convolution of a scaled window with the signal.
#     The signal is prepared by introducing reflected copies of the signal
#     (with the window size) in both ends so that transient parts are minimized
#     in the begining and end part of the output signal.
#
#     input:
#         x: the input signal
#         window_len: the dimension of the smoothing window; should be an odd integer
#         window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
#             flat window will produce a moving average smoothing.
#
#     output:
#         the smoothed signal
#
#     example:
#
#     t=linspace(-2,2,0.1)
#     x=sin(t)+randn(len(t))*0.1
#     y=smooth(x)
#
#     see also:
#
#     numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
#     scipy.signal.lfilter
#
#     TODO: the window parameter could be the window itself if an array instead of a string
#     NOTE: length(output) != length(input), to correct this: return y[(window_len/2-1):-(window_len/2)] instead of just y.
#     """
#     import numpy
#     if x.ndim != 1:
#         raise ValueError("smooth only accepts 1 dimension arrays.")
#
#     if x.size < window_len:
#         raise ValueError("Input vector needs to be bigger than window size.")
#
#     if window_len < 3:
#         return x
#
#     if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
#         raise ValueError("Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'")
#
#     s = numpy.r_[x[window_len - 1:0:-1], x, x[-2:-window_len - 1:-1]]
#     # print(len(s))
#     if window == 'flat':  # moving average
#         w = numpy.ones(window_len, 'd')
#     else:
#         w = eval('numpy.' + window + '(window_len)')
#
#     y = numpy.convolve(w / w.sum(), s, mode=mode)
#     if mode == 'same':
#         return y[np.int(window_len) - 1:-np.int(window_len) + 1]
#     else:
#         return y[np.int(window_len / 2 - 1):-np.int(window_len / 2)]

def img2movie(imgprefix='', img_ext='png', outname='movie', size=None, start_num=0, crf=15, fps=10, overwrite=False,
              crop=[], title=[], dpi=200, keeptmp=False, usetmp=False, autorotate=True):
    '''

    :param imgprefix:
    :param img_ext:
    :param outname:
    :param size:
    :param start_num:
    :param crf:
    :param fps:
    :param overwrite:
    :param title: the timestamp on each frame
    :param crop: 4-tuple of integer specifies the cropping pixels [x0, x1, y0, y1]
    :param dpi:
    :param keeptmp:
    :param usetmp: use the image in the default tmp folder
    :return:
    '''
    import subprocess, os
    from tqdm import tqdm
    import matplotlib.pyplot as plt
    if type(imgprefix) is list:
        imgs = imgprefix
        imgprefix = os.path.dirname(imgprefix[0])
    else:
        imgs = glob.glob(imgprefix + '*.' + img_ext)
    if imgs:
        imgs = sorted(imgs)
        tmpdir = os.path.join(os.path.dirname(imgprefix), 'img_tmp') + '/'
        if not os.path.exists(tmpdir):
            os.makedirs(tmpdir)
        if usetmp:
            pass
        else:
            if crop != [] or title != []:
                plt.ioff()
                fig, ax = plt.subplots(figsize=[float(ll) / dpi for ll in size.split('x')])
                for idx, ll in enumerate(tqdm(imgs)):
                    data = plt.imread(ll)
                    if crop != []:
                        x0, x1, y0, y1 = crop
                        data = data[y0:y1 + 1, x0:x1 + 1, :]
                    if idx == 0:
                        im = ax.imshow(data)
                        ax.set_axis_off()
                    else:
                        im.set_array(data)
                    if title != []:
                        ax.set_title(title[idx])
                    if idx == 0:
                        fig.tight_layout()
                    fig.savefig('{}/{:04d}.{}'.format(tmpdir, idx, img_ext), dpi=dpi, format=img_ext)
                plt.close(fig)
            else:
                for idx, ll in enumerate(imgs):
                    os.system('cp {} {}/{:04d}.{}'.format(ll, tmpdir, idx, img_ext))
        outd = {'r': fps, 's': size, 'start_number': start_num, 'crf': crf}
        if size is None:
            outd.pop('s')
        if overwrite:
            ow = '-y'
        else:
            ow = ''
        if autorotate:
            noautorotate = ''
        else:
            noautorotate = '-noautorotate'
        try:
            outdstr = ' '.join(['-{} {}'.format(k, v) for k, v in outd.iteritems()])
        except:
            outdstr = ' '.join(['-{} {}'.format(k, outd[k]) for k in outd])
        try:
            cmd = 'ffmpeg {4} -r {3} -f image2 -i {0}%04d.{1} -vcodec libx264 -pix_fmt yuv420p {2} '.format(tmpdir,
                                                                                                            img_ext,
                                                                                                            outdstr,
                                                                                                            fps,
                                                                                                            noautorotate) + '{0} {1}.mp4'.format(
                ow, os.path.join(os.path.dirname(imgprefix), outname))
            subprocess.check_output(['bash', '-c', cmd])
        except:
            cmd = 'ffmpeg {4} -r {3} -f image2 -i {0}%04d.{1} -vcodec libx264 -pix_fmt yuv420p {2} -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2"'.format(
                tmpdir, img_ext,
                outdstr,
                fps, noautorotate) + '{0} {1}.mp4'.format(
                ow, os.path.join(os.path.dirname(imgprefix), outname))
            subprocess.check_output(['bash', '-c', cmd])
        print(cmd)
        if not keeptmp:
            os.system('rm -rf {}'.format(tmpdir))
    else:
        print('Images not found!')


def image_fill_gap(image):
    for idx in range(image.shape[0]):
        image_slice = image[idx, :]
        mask_nan = np.isnan(image_slice)
        if np.sum(~mask_nan) > 0 and np.sum(mask_nan) > 0:
            image_slice[mask_nan] = np.interp(np.flatnonzero(mask_nan), np.flatnonzero(~mask_nan),
                                              image_slice[~mask_nan])
            image[idx, :] = image_slice
    for idx in range(image.shape[1]):
        image_slice = image[:, idx]
        mask_nan = np.isnan(image_slice)
        if np.sum(~mask_nan) > 0 and np.sum(mask_nan) > 0:
            image_slice[mask_nan] = np.interp(np.flatnonzero(mask_nan), np.flatnonzero(~mask_nan),
                                              image_slice[~mask_nan])
            image[:, idx] = image_slice
    return image


def getspwfromfreq(vis, freqrange):
    from taskinit import ms
    ms.open(vis)
    axisInfo = ms.getdata(["axis_info"], ifraxis=True)
    spwInfo = ms.getspectralwindowinfo()
    freqInfo = axisInfo["axis_info"]["freq_axis"]["chan_freq"].swapaxes(0, 1) / 1e9
    freqInfo_ravel = freqInfo.ravel()
    timeInfo = axisInfo["axis_info"]["time_axis"]['MJDseconds']
    mstimran = ms.range(["time"])
    ms.close()
    freq0, freq1 = freqrange.split(' ')[0].split('~')
    freq0, freq1 = float(freq0), float(freq1)
    for ll in [freq0, freq1]:
        if not freqInfo_ravel[0] <= ll <= freqInfo_ravel[-1]:
            raise ValueError('Selected frequency out of range!!!')
    freqIdx0 = np.where(freqInfo == freq0)
    freqIdx1 = np.where(freqInfo == freq1)
    sz_freqInfo = freqInfo.shape
    ms_spw = ['{}'.format(ll) for ll in range(freqIdx0[0], freqIdx1[0] + 1)]
    if len(ms_spw) == 1:
        ms_chan = ['{}~{}'.format(freqIdx0[1][0], freqIdx1[1][0])]
    else:
        ms_chan = ['{}~{}'.format(freqIdx0[1][0], sz_freqInfo[1] - 1)] + ['0~{}'.format(sz_freqInfo[1] - 1) for ll in
                                                                          range(freqIdx0[0] + 1, freqIdx1[0])]
        ms_chan.append('0~{}'.format(freqIdx1[1][0]))
    spw = ','.join('{}:{}'.format(t[0], t[1]) for t in zip(ms_spw, ms_chan))
    return spw


def initconfig(suncasa_dir):
    if not os.path.exists(suncasa_dir + 'DataBrowser/config.json'):
        os.system(
            'cp {} {}'.format(suncasa_dir + 'DataBrowser/config_init.json', suncasa_dir + 'DataBrowser/config.json'))
        return True
    else:
        return False


# def mkunidir(dirlist, isdir=True):
#     '''
#     get the list of all unique directories from the filelist, make these directories
#     :param dirlist:
#     :param isdir: if dirlist is a list of directories.
#                 if not, dirlist is a list of full path of files,
#                 the base name will be removed from the paths.
#     :return:
#     '''
#   DEll19432017  import os
#     if not isdir:
#         dirlist = [os.path.dirname(ff) for ff in dirlist]
#     dirs = list(set(dirlist))
#     for ll in dirs:
#         if not os.path.exists(ll):
#             os.makedirs(ll)

def ProgressBar(iteration, total, prefix='', suffix='', decimals=1, length=100, empfill=' ', fill='#'):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
        empfill     - Optional  : empty bar fill character (Str)
    """
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + empfill * (length - filledLength)
    # return '%s |%s| %s%% %s' % (prefix, bar, percent, suffix)
    return '{} |{}| {}% {}'.format(prefix, bar, percent, suffix)


def getcurtimstr(prefix='CleanID_', suffix=''):
    import time
    return prefix + time.strftime("%Y%m%d_%H%M%S") + suffix


def getlatestfile(directory='./', prefix='CleanID_', suffix=''):
    filelist = glob.glob('{}/{}*{}'.format(directory, prefix, suffix))
    if len(filelist) > 0:
        latest_file = max(filelist, key=os.path.getctime)
        print(latest_file)
        return {'items': filelist, 'latest': latest_file, 'idx': filelist.index(latest_file)}
    else:
        print('No file found!')
        return None


def loadjsonfile(jsonfile, mustexist=True):
    if os.path.exists(jsonfile):
        with open(jsonfile, 'r') as fp:
            data = json.load(fp)
        return data
    else:
        if mustexist:
            raise SystemExit('{} not found!'.format(jsonfile))
        else:
            return None


def updatejsonfile(jsonfile, data):
    with open(jsonfile, 'w') as fp:
        json.dump(data, fp)


def getSDOdir(config, database_dir, suncasa_dir):
    try:
        if config['datadir']['SDOdir']:
            SDOdir = config['datadir']['SDOdir']
            if SDOdir.startwith('$'):
                SDOdir = os.path.expandvars(SDOdir)
            if not os.path.exists(SDOdir):
                os.makedirs(SDOdir)
        else:
            raise ValueError
    except:
        SDOdir = database_dir + 'Download/'
        config['datadir']['SDOdir'] = SDOdir
        fout = suncasa_dir + 'DataBrowser/config.json'
        updatejsonfile(fout, config)
    return SDOdir


def getsdodir(filename, unique=True):
    '''
    return a list of the data path relative to the SDOdir
    :param filename:
    :return:
    '''
    if type(filename) == str:
        filename = [filename]
    dirlist = []
    timstrs = []
    if unique:
        for ll in filename:
            timstr = os.path.basename(ll).split('.')[2].split('T')[0]
            ymd = timstr.replace('-', '/')
            dirlist.append('/{}/_{}'.format(ymd, timstr))
        dirlist = list(set(dirlist))
        for ll, dd in enumerate(dirlist):
            dirlist[ll], timstr = dd.split('_')
            timstrs.append(timstr)
    else:
        for ll in filename:
            timstr = os.path.basename(ll).split('.')[2].split('T')[0]
            ymd = timstr.replace('-', '/')
            dirlist.append('/{}/'.format(ymd))
            timstrs.append(timstr)
    return {'dir': dirlist, 'timstr': timstrs}


def FileNotInList(file2chk, filelist):
    '''
    return the index of files not in the list
    :param file2chk: files to be check
    :param filelist: the list
    :return:
    '''
    idxs = []
    if len(file2chk) > 0:
        filelist = [os.path.basename(ll) for ll in filelist]
        for idx, items in enumerate(file2chk):
            if items not in filelist:
                idxs.append(idx)
    return idxs


def getfreeport():
    import socket
    s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    s.bind(('localhost', 0))
    port = s.getsockname()[1]
    s.close()
    return port


def normalize_aiamap(aiamap):
    '''
    do expisure normalization of an aia map
    :param aia map made from sunpy.map:
    :return: normalised aia map
    '''
    import sunpy.map as smap
    try:
        if aiamap.observatory == 'SDO' and aiamap.instrument[0:3] == 'AIA':
            data = aiamap.data.copy().astype(np.float)
            idxpix = ~np.isnan(data)
            data[idxpix] = data[idxpix] / aiamap.exposure_time.value
            data[data < 0] = 0
            aiamap.meta['exptime'] = 1.0
            aiamap = smap.Map(data, aiamap.meta)
            return aiamap
        else:
            raise ValueError('input sunpy map is not from aia.')
    except:
        raise ValueError('check your input map. There are some errors in it.')


def tplt(mapcube):
    from astropy.time import Time
    t = []
    for idx, mp in enumerate(mapcube):
        if mp.meta.has_key('t_obs'):
            tstr = mp.meta['t_obs']
        else:
            tstr = mp.meta['date-obs']
        t.append(tstr)
    return Time(t)


def sdo_aia_scale_hdr(amap, sigma=None):
    import sunpy.map as smap
    from astropy.coordinates import SkyCoord
    wavelnth = '{:0.0f}'.format(amap.wavelength.value)
    pxscale_x, pxscal_y = amap.scale
    pxscale = (pxscale_x + pxscal_y) / 2
    r_sun = amap.rsun_obs / pxscale
    x = np.arange(amap.dimensions.x.value)
    y = np.arange(amap.dimensions.y.value)
    xx, yy = np.meshgrid(x, y) * u.pix
    crpix1, crpix2 = amap.world_to_pixel(SkyCoord(0 * u.arcsec, 0 * u.arcsec, frame=amap.coordinate_frame))
    rdist = np.sqrt((xx - (crpix1 - 1.0 * u.pix)) ** 2 + (yy - (crpix2 - 1.0 * u.pix)) ** 2)
    ind_disk = np.where(rdist <= r_sun)
    # ind_limb = np.where(rdist < r_sun)
    rdist[ind_disk] = r_sun
    rfilter = rdist / r_sun - 1
    rfilter = rfilter.value
    if sigma:
        mapdata = amap.data * np.exp(rfilter * sigma)
    else:
        if wavelnth == '94':
            mapdata = amap.data * np.exp(rfilter * 4)
        elif wavelnth == '131':
            mapdata = amap.data * (np.sqrt(rfilter * 5) + 1)
        elif wavelnth == '171':
            mapdata = amap.data * np.exp(rfilter * 5)
        elif wavelnth == '193':
            mapdata = amap.data * np.exp(rfilter * 3)
        elif wavelnth == '211':
            mapdata = amap.data * np.exp(rfilter * 3)
        elif wavelnth == '304':
            mapdata = amap.data * np.exp(rfilter * 5)
        elif wavelnth == '335':
            mapdata = amap.data * np.exp(rfilter * 3)
        elif wavelnth == '6173':
            pass
        elif wavelnth == '1':
            pass
        else:
            sigma = 5.0
            mapdata = amap.data * np.exp(rfilter * sigma)
    return smap.Map(mapdata, amap.meta)


def sdo_aia_scale_dict(wavelength=None, imagetype='image'):
    '''
    rescale the aia image
    :param image: normalised aia image data
    :param wavelength:
    :return: byte scaled image data
    '''
    if type(wavelength) is not str:
        wavelength = '{:0.0f}'.format(wavelength)
    if wavelength == '94':
        if imagetype == 'image':
            return {'low': 0.1, 'high': 3000, 'log': True}
        elif imagetype == 'RDimage':
            return {'low': -30, 'high': 30, 'log': False}
        elif imagetype == 'BDimage':
            return {'low': -30, 'high': 30, 'log': False}
        elif imagetype == 'RDRimage':
            return {'low': -1.5, 'high': 1.5, 'log': False}
        elif imagetype == 'BDRimage':
            return {'low': -1.5, 'high': 1.5, 'log': False}
    elif wavelength == '131':
        if imagetype == 'image':
            return {'low': 0.5, 'high': 10000, 'log': True}
        elif imagetype == 'RDimage':
            return {'low': -100, 'high': 100, 'log': False}
        elif imagetype == 'BDimage':
            return {'low': -100, 'high': 100, 'log': False}
        elif imagetype == 'RDRimage':
            return {'low': -1.5, 'high': 1.5, 'log': False}
        elif imagetype == 'BDRimage':
            return {'low': -1.5, 'high': 1.5, 'log': False}
    elif wavelength == '171':
        if imagetype == 'image':
            return {'low': 20, 'high': 5000, 'log': True}
        elif imagetype == 'RDimage':
            return {'low': -400, 'high': 400, 'log': False}
        elif imagetype == 'BDimage':
            return {'low': -400, 'high': 400, 'log': False}
        elif imagetype == 'RDRimage':
            return {'low': -1.5, 'high': 1.5, 'log': False}
        elif imagetype == 'BDRimage':
            return {'low': -1.5, 'high': 1.5, 'log': False}
    elif wavelength == '193':
        if imagetype == 'image':
            return {'low': 30, 'high': 5000, 'log': True}
        elif imagetype == 'RDimage':
            return {'low': -1500, 'high': 1500, 'log': False}
        elif imagetype == 'BDimage':
            return {'low': -1500, 'high': 1500, 'log': False}
        elif imagetype == 'RDRimage':
            return {'low': -1.5, 'high': 1.5, 'log': False}
        elif imagetype == 'BDRimage':
            return {'low': -1.5, 'high': 1.5, 'log': False}
    elif wavelength == '211':
        if imagetype == 'image':
            return {'low': 10, 'high': 8000, 'log': True}
        elif imagetype == 'RDimage':
            return {'low': -300, 'high': 300, 'log': False}
        elif imagetype == 'BDimage':
            return {'low': -300, 'high': 300, 'log': False}
        elif imagetype == 'RDRimage':
            return {'low': -1.5, 'high': 1.5, 'log': False}
        elif imagetype == 'BDRimage':
            return {'low': -1.5, 'high': 1.5, 'log': False}
    elif wavelength == '304':
        if imagetype == 'image':
            return {'low': 1, 'high': 10000, 'log': True}  # return {'low': 1, 'high': 500, 'log': True}
        elif imagetype == 'RDimage':
            return {'low': -300, 'high': 300, 'log': False}
        elif imagetype == 'BDimage':
            return {'low': -300, 'high': 300, 'log': False}
        elif imagetype == 'RDRimage':
            return {'low': -1.5, 'high': 1.5, 'log': False}
        elif imagetype == 'BDRimage':
            return {'low': -1.5, 'high': 1.5, 'log': False}
    elif wavelength == '335':
        if imagetype == 'image':
            return {'low': 0.1, 'high': 800, 'log': True}
        elif imagetype == 'RDimage':
            return {'low': -15, 'high': 15, 'log': False}
        elif imagetype == 'BDimage':
            return {'low': -15, 'high': 15, 'log': False}
        elif imagetype == 'RDRimage':
            return {'low': -1.5, 'high': 1.5, 'log': False}
        elif imagetype == 'BDRimage':
            return {'low': -1.5, 'high': 1.5, 'log': False}
    elif wavelength == '1600':
        if imagetype == 'image':
            return {'low': 20, 'high': 2500, 'log': True}
        elif imagetype == 'RDimage':
            return {'low': -800, 'high': 800, 'log': False}
        elif imagetype == 'BDimage':
            return {'low': -800, 'high': 800, 'log': False}
        elif imagetype == 'RDRimage':
            return {'low': -1.5, 'high': 1.5, 'log': False}
        elif imagetype == 'BDRimage':
            return {'low': -1.5, 'high': 1.5, 'log': False}
    elif wavelength == '1700':
        if imagetype == 'image':
            return {'low': 300, 'high': 10000, 'log': True}
        elif imagetype == 'RDimage':
            return {'low': -1500, 'high': 1500, 'log': False}
        elif imagetype == 'BDimage':
            return {'low': -1500, 'high': 1500, 'log': False}
        elif imagetype == 'RDRimage':
            return {'low': -1.5, 'high': 1.5, 'log': False}
        elif imagetype == 'BDRimage':
            return {'low': -1.5, 'high': 1.5, 'log': False}
    else:
        return {'high': None, 'log': False, 'low': None}


def sdo_aia_scale(image=None, wavelength=None):
    '''
    rescale the aia image
    :param image: normalised aia image data
    :param wavelength:
    :return: byte scaled image data
    '''
    from scipy.misc import bytescale
    clrange = sdo_aia_scale_dict(wavelength)
    image[image > clrange['high']] = clrange['high']
    image[image < clrange['low']] = clrange['low']
    image = np.log10(image)
    return bytescale(image)


def insertchar(source_str, insert_str, pos):
    return source_str[:pos] + insert_str + source_str[pos:]


# def get_trange_files(trange):
#     '''
#     Given a timerange, this routine will take all relevant SDO files from that time range,
#     put them in a list, and return that list.
#     :param trange: two elements list of timestamps in Julian days
#     :return:
#     '''


def readsdofile(datadir=None, wavelength=None, trange=None, isexists=False, timtol=1):
    '''
    read sdo file from local database
    :param datadir:
    :param wavelength:
    :param trange: the timestamp or timerange in Julian days. if is timerange, return a list of files in the timerange
    :param isexists: check if file exist. if files exist, return file name
    :param timtol: time difference tolerance in days for considering data as the same timestamp
    :return:
    '''
    from astropy.time import Time
    import sunpy.map
    from datetime import date
    from datetime import timedelta as td

    trange = Time(trange)
    wavelength = str(wavelength)
    wavelength = wavelength.lower()
    if timtol < 12. / 3600 / 24:
        timtol = 12. / 3600 / 24

    if len(trange.iso) == 2:
        if trange[1] < trange[0]:
            raise ValueError('start time must be occur earlier than end time!')
        else:
            sdofitspath = []
            jdtimestr = [Time(ll, format='jd').iso for ll in trange]
            ymd = [ll.split(' ')[0].split('-') for ll in jdtimestr]
            d1 = date(int(ymd[0][0]), int(ymd[0][1]), int(ymd[0][2]))
            d2 = date(int(ymd[1][0]), int(ymd[1][1]), int(ymd[1][2]))
            delta = d2 - d1
            for i in range(delta.days + 1):
                ymd = d1 + td(days=i)
                sdofitspathtmp = glob.glob(
                    datadir + '/{:04d}/{:02d}/{:02d}/aia.lev1_*Z.{}.image_lev1.fits'.format(ymd.year, ymd.month,
                                                                                            ymd.day, wavelength))
                if len(sdofitspathtmp) > 0:
                    sdofitspath = sdofitspath + sdofitspathtmp
            if len(sdofitspath) == 0:
                if isexists:
                    return sdofitspath
                else:
                    raise ValueError(
                        'No SDO file found under {} at the time range of {} to {}. Download the data with EvtBrowser first.'.format(
                            datadir,
                            jdtimestr[0],
                            jdtimestr[1]))
            sdofits = [os.path.basename(ll) for ll in sdofitspath]
            sdotimeline = Time(
                [insertchar(insertchar(ll.split('.')[2].replace('T', ' ').replace('Z', ''), ':', -4), ':', -2) for ll in
                 sdofits],
                format='iso', scale='utc')
            sdofitspathnew = [x for (y, x) in sorted(zip(sdotimeline.jd, sdofitspath))]
            sdofitsnew = [os.path.basename(ll) for ll in sdofitspathnew]
            sdotimelinenew = Time(
                [insertchar(insertchar(ll.split('.')[2].replace('T', ' ').replace('Z', ''), ':', -4), ':', -2) for ll in
                 sdofitsnew], format='iso',
                scale='utc')
            sdofile = list(
                np.array(sdofitspathnew)[
                    np.where(np.logical_and(trange[0].jd <= sdotimelinenew.jd, sdotimelinenew.jd <= trange[1].jd))[0]])
            return sdofile
    else:
        jdtimstr = trange.iso
        ymd = jdtimstr.split(' ')[0].split('-')
        sdofitspath = glob.glob(
            datadir + '/{}/{}/{}/aia.lev1_*Z.{}.image_lev1.fits'.format(ymd[0], ymd[1], ymd[2], wavelength))
        if len(sdofitspath) == 0:
            return []  # raise ValueError('No SDO file found under {}.'.format(datadir))
        sdofits = [os.path.basename(ll) for ll in sdofitspath]
        sdotimeline = Time(
            [insertchar(insertchar(ll.split('.')[2].replace('T', ' ').replace('Z', ''), ':', -4), ':', -2) for ll in
             sdofits],
            format='iso', scale='utc')
        if timtol < np.min(np.abs(sdotimeline.jd - trange.jd)):
            return []  # raise ValueError('No SDO file found at the select timestamp. Download the data with EvtBrowser first.')
        idxaia = np.argmin(np.abs(sdotimeline.jd - trange.jd))
        sdofile = sdofitspath[idxaia]
        if isexists:
            return sdofile
        else:
            try:
                sdomap = sunpy.map.Map(sdofile)
                return sdomap
            except:
                raise ValueError('File not found or invalid input')


def readsdofileX(datadir=None, filelist=None, wavelength=None, trange=None, isexists=False, timtol=1):
    '''
    read sdo file from local database
    :param datadir:
    :param wavelength:
    :param trange: time object. the timestamp or timerange. if is timerange, return a list of files in the timerange
    :param isexists: check if file exist. if files exist, return file name
    :param timtol: time difference tolerance in days for considering data as the same timestamp
    :return:
    '''
    from astropy.time import Time
    import sunpy.map
    from datetime import date
    from datetime import timedelta as td

    trange = Time(trange)
    wavelength = str(wavelength)
    wavelength = wavelength.lower()
    if timtol < 12. / 3600 / 24:
        timtol = 12. / 3600 / 24
    # import pdb
    # pdb.set_trace()
    if len(trange.iso) == 2:
        if trange.jd[1] < trange.jd[0]:
            raise ValueError('start time must be occur earlier than end time!')
        else:
            if filelist:
                sdofitspath = filelist
            else:
                sdofitspath = []
                jdtimestr = trange.iso
                ymd = [ll.split(' ')[0].split('-') for ll in jdtimestr]
                d1 = date(int(ymd[0][0]), int(ymd[0][1]), int(ymd[0][2]))
                d2 = date(int(ymd[1][0]), int(ymd[1][1]), int(ymd[1][2]))
                delta = d2 - d1
                for i in range(delta.days + 1):
                    ymd = d1 + td(days=i)
                    sdofitspathtmp = glob.glob(
                        datadir + '/aia.lev1_*{0}*{1}*{2}*Z.{3}.image*.fits'.format(ymd.year, ymd.month, ymd.day,
                                                                                    wavelength))
                    if len(sdofitspathtmp) > 0:
                        sdofitspath = sdofitspath + sdofitspathtmp
            if len(sdofitspath) == 0:
                if isexists:
                    return sdofitspath
                else:
                    raise ValueError(
                        'No SDO file found under {} at the time range of {} to {}. Download the data with EvtBrowser first.'.format(
                            datadir,
                            jdtimestr[0],
                            jdtimestr[1]))
            sdofits = [os.path.basename(ll) for ll in sdofitspath]
            if 'hmi' in sdofits[0]:
                sdotimeline = Time(
                    [insertchar(insertchar(
                        insertchar(insertchar(ll.split('.')[2].replace('_TAI', '').replace('_', ' '), ':', -4), ':',
                                   -2), '-',
                        6), '-', 4) for ll in
                        sdofits],
                    format='iso', scale='utc')
            else:
                sdotimeline = Time(
                    [insertchar(insertchar(ll.split('.')[2].replace('T', ' ').replace('Z', ''), ':', -4), ':', -2) for
                     ll in
                     sdofits],
                    format='iso', scale='utc')
            sdofitspathnew = [x for (y, x) in sorted(zip(sdotimeline.jd, sdofitspath))]
            sdofitsnew = [os.path.basename(ll) for ll in sdofitspathnew]
            if 'hmi' in sdofits[0]:
                sdotimelinenew = Time(
                    [insertchar(insertchar(
                        insertchar(insertchar(ll.split('.')[2].replace('_TAI', '').replace('_', ' '), ':', -4), ':',
                                   -2), '-',
                        6), '-', 4) for ll in
                        sdofitsnew],
                    format='iso', scale='utc')
            else:
                sdotimelinenew = Time(
                    [insertchar(insertchar(ll.split('.')[2].replace('T', ' ').replace('Z', ''), ':', -4), ':', -2) for
                     ll in
                     sdofitsnew],
                    format='iso', scale='utc')
            sdofile = list(
                np.array(sdofitspathnew)[
                    np.where(np.logical_and(trange.jd[0] <= sdotimelinenew.jd, sdotimelinenew.jd <= trange.jd[1]))[0]])
            return sdofile
    else:
        jdtimstr = trange.iso
        ymd = jdtimstr.split(' ')[0].split('-')
        if filelist:
            sdofitspath = filelist
        else:
            sdofitspath = glob.glob(
                datadir + '/aia.lev1_*{0}*{1}*{2}*Z.{3}.image*.fits'.format(ymd[0], ymd[1], ymd[2], wavelength))
        if len(sdofitspath) == 0:
            return []  # raise ValueError('No SDO file found under {}.'.format(datadir))
        sdofits = [os.path.basename(ll) for ll in sdofitspath]
        if 'hmi' in sdofits[0]:
            sdotimeline = Time(
                [insertchar(insertchar(
                    insertchar(insertchar(ll.split('.')[2].replace('_TAI', '').replace('_', ' '), ':', -4), ':', -2),
                    '-',
                    6), '-', 4) for ll in
                    sdofits],
                format='iso', scale='utc')
        else:
            sdotimeline = Time(
                [insertchar(insertchar(ll.split('.')[2].replace('T', ' ').replace('Z', ''), ':', -4), ':', -2) for ll in
                 sdofits],
                format='iso', scale='utc')
        if timtol <= np.min(np.abs(sdotimeline.jd - trange.jd)):
            return []  # raise ValueError('No SDO file found at the select timestamp. Download the data with EvtBrowser first.')
        idxaia = np.argmin(np.abs(sdotimeline.jd - trange.jd))
        sdofile = sdofitspath[idxaia]
        if isexists:
            return sdofile
        else:
            try:
                sdomap = sunpy.map.Map(sdofile)
                return sdomap
            except:
                raise ValueError('File not found or invalid input')


def findDist(x, y):
    dx = np.diff(x)
    dy = np.diff(y)
    dist = np.hypot(dx, dy)
    return np.insert(dist, 0, 0.0)


def paramspline(x, y, length, s=0):
    from scipy.interpolate import splev, splprep
    tck, u = splprep([x, y], s=s)
    unew = np.linspace(0, u[-1], length)
    out = splev(unew, tck)
    xs, ys = out[0], out[1]
    grads = get_curve_grad(xs, ys)
    return {'xs': xs, 'ys': ys, 'grads': grads['grad'], 'posangs': grads['posang']}


def polyfit(x, y, length, deg, keepxorder=False):
    if keepxorder:
        xs = np.linspace(x[0], x[-1], length)
    else:
        xs = np.linspace(np.nanmin(x), np.nanmax(x), length)
    z = np.polyfit(x=x, y=y, deg=deg)
    p = np.poly1d(z)
    ys = p(xs)
    if keepxorder:
        return {'xs': xs, 'ys': ys}
    else:
        grads = get_curve_grad(xs, ys)
        return {'xs': xs, 'ys': ys, 'grads': grads['grad'], 'posangs': grads['posang']}


def htfit_warren2011(x, y, cutlength):
    from scipy.optimize import curve_fit

    def fit_func(t, h0, vt, a0, tau):
        return h0 + vt * t + a0 * tau ** 2 * (np.exp(-t / tau) - 1)

    x0 = x[0]
    params = curve_fit(fit_func, x - x0, y)

    [h0, vt, a0, tau] = params[0]
    xs = np.linspace(np.nanmin(x), np.nanmax(x), cutlength)
    ys = fit_func(xs - x0, h0, vt, a0, tau)
    grads = get_curve_grad(xs, ys)
    return {'xs': xs, 'ys': ys, 'grads': grads['grad'], 'posangs': grads['posang']}


def spline(x, y, length, s=0):
    from scipy.interpolate import splev, splrep
    xs = np.linspace(x.min(), x.max(), length)
    tck = splrep(x, y, s=s)
    ys = splev(xs, tck)
    grads = get_curve_grad(xs, ys)
    return {'xs': xs, 'ys': ys, 'grads': grads['grad'], 'posangs': grads['posang']}


def get_curve_grad(x, y):
    '''
    get the grad of at data point
    :param x:
    :param y:
    :return: grad,posang
    '''
    deltay = np.roll(y, -1) - np.roll(y, 1)
    deltay[0] = y[1] - y[0]
    deltay[-1] = y[-1] - y[-2]
    deltax = np.roll(x, -1) - np.roll(x, 1)
    deltax[0] = x[1] - x[0]
    deltax[-1] = x[-1] - x[-2]
    grad = deltay / deltax
    posang = np.arctan2(deltay, deltax)
    return {'grad': grad, 'posang': posang}


def improfile(z, xi, yi, interp='cubic'):
    '''
    Pixel-value cross-section along line segment in an image
    :param z: an image array
    :param xi and yi: equal-length vectors specifying the pixel coordinates of the endpoints of the line segment
    :param interp: interpolation type to sampling, 'nearest' or 'cubic'
    :return: the intensity values of pixels along the line
    '''
    import scipy.ndimage
    imgshape = z.shape
    if len(xi) != len(yi):
        raise ValueError('xi and yi must be equal-length!')
    if len(xi) < 2:
        raise ValueError('xi or yi must contain at least two elements!')
    for idx, ll in enumerate(xi):
        if not 0 < ll < imgshape[1]:
            return np.nan
        if not 0 < yi[idx] < imgshape[0]:
            return np.nan
    if len(xi) == 2:
        length = np.hypot(np.diff(xi), np.diff(yi))[0]
        x, y = np.linspace(xi[0], xi[1], length), np.linspace(yi[0], yi[1], length)
    else:
        x, y = xi, yi
    if interp == 'cubic':
        zi = scipy.ndimage.map_coordinates(z, np.vstack((x, y)))
    else:
        zi = z[np.floor(y).astype(np.int), np.floor(x).astype(np.int)]

    return zi


def canvaspix_to_data(smap, x, y):
    import astropy.units as u
    '''
    Convert canvas pixel coordinates in MkPlot to data (world) coordinates by using
    `~astropy.wcs.WCS.wcs_pix2world`.

    :param smap: sunpy map
    :param x: canvas Pixel coordinates of the CTYPE1 axis. (Normally solar-x)
    :param y: canvas Pixel coordinates of the CTYPE2 axis. (Normally solar-y)
    :return: world coordinates
    '''
    # xynew = smap.pixel_to_data(x * u.pix, y * u.pix)
    mesh = smap.pixel_to_world(x * u.pix, y * u.pix)
    mapx, mapy = mesh.Tx, mesh.Tx
    xnew = mapx[0].value
    ynew = mapy[1].value
    return [xnew, ynew]


def data_to_mappixel(smap, x, y):
    import astropy.units as u
    '''
    Convert data (world) coordinates in MkPlot to pixel coordinates in smap by using
    `~astropy.wcs.WCS.wcs_pix2world`.

    :param smap: sunpy map
    :param x: Data coordinates of the CTYPE1 axis. (Normally solar-x)
    :param y: Data coordinates of the CTYPE2 axis. (Normally solar-y)
    :return: pixel coordinates
    '''
    xynew = smap.data_to_pixel(x * u.arcsec, y * u.arcsec)
    xnew = xynew[0].value
    ynew = xynew[1].value
    return [xnew, ynew]


def polsfromfitsheader(header):
    '''
    get polarisation information from fits header
    :param header: fits header
    :return pols: polarisation stokes
    '''
    try:
        stokeslist = ['{}'.format(int(ll)) for ll in
                      (header["CRVAL4"] + np.arange(header["NAXIS4"]) * header["CDELT4"])]
        stokesdict = {'1': 'I', '2': 'Q', '3': 'U', '4': 'V', '-1': 'RR', '-2': 'LL', '-3': 'RL', '-4': 'LR',
                      '-5': 'XX', '-6': 'YY', '-7': 'XY',
                      '-8': 'YX'}
        pols = map(lambda x: stokesdict[x], stokeslist)
    except:
        print("error in fits header!")
    return pols


def headerfix(header):
    ## this code fix the header problem of fits out from CASA 5.4+ which leads to a streched solar image
    import copy
    hdr = copy.copy(header)
    for hd in header:
        if hd.upper().startswith('PC'):
            if not hd.upper().startswith('PC0'):
                hd_ = 'PC0' + hd.upper().replace('PC', '')
                hdr.pop(hd)
                hdr[hd_] = header[hd]
    return hdr


def freqsfromfitsheader(header):
    '''
    get frequency in GHz from fits header
    :param header: fits header
    :return pols: polarisation stokes
    '''
    try:
        freqs = ['{:.3f}'.format(ll / 1e9) for ll in
                 (header["CRVAL3"] + np.arange(header["NAXIS3"]) * header["CDELT3"])]
        return freqs
    except:
        print("error in fits header!")
        raise ValueError


def transfitdict2DF(datain, gaussfit=True, getcentroid=False):
    '''
    convert the results from pimfit or pmaxfit tasks to pandas DataFrame structure.
    :param datain: The component list from pimfit or pmaxfit tasks
    :param gaussfit: True if the results is from pimfit, otherwise False.
    :param getcentroid: If True returns the centroid
    :return: the pandas DataFrame structure.
    '''
    import pandas as pd

    ra2arcsec = 180. * 3600. / np.pi
    dspecDF0 = pd.DataFrame()
    for ll in datain['timestamps']:
        tidx = datain['timestamps'].index(ll)
        if datain['succeeded'][tidx]:
            pols = datain['outputs'][tidx].keys()
            dspecDFlist = []
            dspecDF = pd.DataFrame()
            for ppit in pols:
                dspecDFtmp = pd.DataFrame()
                shape_latitude = []
                shape_longitude = []
                shape_latitude_err = []
                shape_longitude_err = []
                shape_majoraxis = []
                shape_minoraxis = []
                shape_positionangle = []
                peak = []
                beam_major = []
                beam_minor = []
                beam_positionangle = []
                freqstrs = []
                fits_local = []
                for comp in datain['outputs'][tidx][ppit]['results'].keys():
                    if comp.startswith('component'):
                        if gaussfit:
                            majoraxis = datain['outputs'][tidx][ppit]['results'][comp]['shape']['majoraxis']['value']
                            minoraxis = datain['outputs'][tidx][ppit]['results'][comp]['shape']['minoraxis']['value']
                            positionangle = datain['outputs'][tidx][ppit]['results'][comp]['shape']['positionangle'][
                                'value']
                            bmajor = datain['outputs'][tidx][ppit]['results'][comp]['beam']['beamarcsec']['major'][
                                'value']
                            bminor = datain['outputs'][tidx][ppit]['results'][comp]['beam']['beamarcsec']['minor'][
                                'value']
                            bpositionangle = \
                                datain['outputs'][tidx][ppit]['results'][comp]['beam']['beamarcsec']['positionangle'][
                                    'value']
                            shape_majoraxis.append(majoraxis)
                            shape_minoraxis.append(minoraxis)
                            shape_positionangle.append(positionangle)
                            beam_major.append(bmajor)
                            beam_minor.append(bminor)
                            beam_positionangle.append(bpositionangle)
                            fluxpeak = datain['outputs'][tidx][ppit]['results'][comp]['peak']['value']
                        else:
                            fluxpeak = datain['outputs'][tidx][ppit]['results'][comp]['flux']['value'][0]
                        if getcentroid:
                            mkey = 'centroid'
                        else:
                            mkey = 'shape'
                        longitude = datain['outputs'][tidx][ppit]['results'][comp][mkey]['direction']['m0'][
                                        'value'] * ra2arcsec
                        latitude = datain['outputs'][tidx][ppit]['results'][comp][mkey]['direction']['m1'][
                                       'value'] * ra2arcsec
                        longitude_err = \
                            datain['outputs'][tidx][ppit]['results'][comp]['shape']['direction']['error']['longitude'][
                                'value']
                        latitude_err = \
                            datain['outputs'][tidx][ppit]['results'][comp]['shape']['direction']['error']['latitude'][
                                'value']
                        shape_longitude.append(longitude)
                        shape_latitude.append(latitude)
                        shape_longitude_err.append(longitude_err)
                        shape_latitude_err.append(latitude_err)
                        peak.append(fluxpeak)
                        freqstrs.append('{:.3f}'.format(
                            datain['outputs'][tidx][ppit]['results'][comp]['spectrum']['frequency']['m0']['value']))
                        fits_local.append(datain['imagenames'][tidx].split('/')[-1])
                if gaussfit:
                    dspecDFtmp['shape_latitude{}'.format(ppit)] = pd.Series(shape_latitude)
                    dspecDFtmp['shape_longitude{}'.format(ppit)] = pd.Series(shape_longitude)
                    dspecDFtmp['shape_latitude_err{}'.format(ppit)] = pd.Series(shape_latitude_err)
                    dspecDFtmp['shape_longitude_err{}'.format(ppit)] = pd.Series(shape_longitude_err)
                    dspecDFtmp['peak{}'.format(ppit)] = pd.Series(peak)
                    dspecDFtmp['shape_majoraxis{}'.format(ppit)] = pd.Series(shape_majoraxis)
                    dspecDFtmp['shape_minoraxis{}'.format(ppit)] = pd.Series(shape_minoraxis)
                    dspecDFtmp['shape_positionangle{}'.format(ppit)] = pd.Series(shape_positionangle)
                    dspecDFtmp['beam_major{}'.format(ppit)] = pd.Series(beam_major)
                    dspecDFtmp['beam_minor{}'.format(ppit)] = pd.Series(beam_minor)
                    dspecDFtmp['beam_positionangle{}'.format(ppit)] = pd.Series(beam_positionangle)
                else:
                    dspecDFtmp['shape_latitude{}'.format(ppit)] = pd.Series(shape_latitude)
                    dspecDFtmp['shape_longitude{}'.format(ppit)] = pd.Series(shape_longitude)
                    dspecDFtmp['shape_latitude_err{}'.format(ppit)] = pd.Series(shape_latitude_err)
                    dspecDFtmp['shape_longitude_err{}'.format(ppit)] = pd.Series(shape_longitude_err)
                    dspecDFtmp['peak{}'.format(ppit)] = pd.Series(peak)
                dspecDFtmp['freqstr'.format(ppit)] = pd.Series(freqstrs)
                dspecDFtmp['fits_local'.format(ppit)] = pd.Series(fits_local)
                dspecDFlist.append(dspecDFtmp)
            for DFidx, DFit in enumerate(dspecDFlist):
                if DFidx == 0:
                    dspecDF = dspecDFlist[0]
                else:
                    dspecDF = pd.merge(dspecDF.copy(), DFit, how='outer', on=['freqstr', 'fits_local'])
            dspecDF0 = dspecDF0.append(dspecDF, ignore_index=True)

    return dspecDF0


def getcolctinDF(dspecDF, col):
    '''
    return the count of how many times of the element starts with col occurs in columns of dspecDF
    :param dspecDF:
    :param col: the start string
    :return: the count and items
    '''
    itemset1 = col
    itemset2 = dspecDF.columns.tolist()
    items = []
    for ll in itemset2:
        if ll.startswith(itemset1):
            items.append(ll)
    itemct = [ll.startswith(itemset1) for ll in itemset2].count(True)
    columns = items
    return [itemct, columns]


def dspecDFfilter(dspecDF, pol):
    '''
    filter the unselect polarisation from dspecDF
    :param dspecDF: the original dspecDF
    :param pol: selected polarisation, dtype = string
    :return: the output dspecDF
    '''
    colnlistcom = ['shape_latitude', 'shape_longitude', 'peak', 'shape_latitude_err', 'shape_longitude_err']
    colnlistgaus = ['shape_majoraxis', 'shape_minoraxis', 'shape_positionangle', 'beam_major', 'beam_minor',
                    'beam_positionangle']
    ## above are the columns to filter
    colnlistess = dspecDF.columns.tolist()
    if getcolctinDF(dspecDF, 'peak')[0] > 0:
        for ll in colnlistcom + colnlistgaus:
            colinfo = getcolctinDF(dspecDF, ll)
            if colinfo[0] > 0:
                for cc in colinfo[1]:
                    if cc in colnlistess:
                        colnlistess.remove(cc)
        dspecDF1 = dspecDF.copy()[colnlistess]
        for ll in colnlistcom:
            dspecDF1[ll] = dspecDF.copy()[ll + pol]
        if getcolctinDF(dspecDF, 'shape_majoraxis')[0] > 0:
            for ll in colnlistgaus:
                dspecDF1[ll] = dspecDF.copy()[ll + pol]
        print('dspedDF is filtered')
        return dspecDF1
    else:
        print('dspedDF no need filter')
        return dspecDF


def dspecDF2text(DFfile, outfile=None):
    if DFfile:
        if os.path.exists(DFfile):
            if outfile:
                with open(DFfile, 'rb') as f:
                    dspecDF0 = pickle.load(f)
                dspecDF0.drop(['dspec', 'fits_global', 'fits_local'], axis=1, inplace=True)
                dspecDF0.to_csv(outfile, sep='\t')
            else:
                raise ValueError('provide output file name!')
        else:
            raise ValueError('input file "{}" does not existed!'.format(DFfile))
    else:
        raise ValueError('provide input file name!')


def smapmeshgrid2(smap, angle=None, rescale=1.0, origin=1):
    import astropy.units as u
    if angle is None:
        mrot = smap.rotation_matrix
    else:
        sin = np.sin(angle)
        cos = np.cos(angle)
        mrot = np.array([[cos, -sin], [sin, cos]])
    ref_pix = smap.reference_pixel
    ref_coord = smap.reference_coordinate
    scale = smap.scale
    XX, YY = np.meshgrid(np.arange(smap.data.shape[1] * rescale) / rescale,
                         np.arange(smap.data.shape[0] * rescale) / rescale)
    x, y = XX * u.pix, YY * u.pix
    try:
        x = (x - ref_pix[0] + origin * u.pix) * scale[0] + ref_coord.x
        y = (y - ref_pix[1] + origin * u.pix) * scale[1] + ref_coord.y
    except:
        x = (x - ref_pix[0] + origin * u.pix) * scale[0] + ref_coord.Tx
        y = (y - ref_pix[1] + origin * u.pix) * scale[1] + ref_coord.Ty
    xnew = mrot[0, 0] * x + mrot[0, 1] * y
    ynew = mrot[1, 0] * x + mrot[1, 1] * y
    return xnew, ynew


def map2wcsgrids(snpmap, cell=True, antialiased=True):
    '''

    :param snpmap:
    :param cell: if True, return the coordinates of the pixel centers. if False, return the coordinates of the pixel boundaries
    :return:
    '''
    # embed()
    import astropy.units as u
    ny, nx = snpmap.data.shape
    x0, x1 = snpmap.xrange.to(u.arcsec).value
    y0, y1 = snpmap.yrange.to(u.arcsec).value
    dx = snpmap.scale.axis1.to(u.arcsec / u.pix).value
    dy = snpmap.scale.axis2.to(u.arcsec / u.pix).value

    if cell:
        mapx, mapy = np.linspace(x0, x1, nx), np.linspace(y0, y1, ny)
        mapx = np.tile(mapx, ny).reshape(ny, nx)
        mapy = np.tile(mapy, nx).reshape(nx, ny).transpose()
    else:
        nx += 1
        ny += 1
        mapx, mapy = np.linspace(x0 - dx / 2.0, x1 + dx / 2.0, nx), np.linspace(y0 - dy / 2.0, y1 + dy / 2.0, ny)
        mapx = np.tile(mapx, ny).reshape(ny, nx)
        mapy = np.tile(mapy, nx).reshape(nx, ny).transpose()
    return mapx, mapy


def smapradialfilter(sunpymap, factor=5, grid=None):
    from sunpy import map as smap
    if grid:
        x, y = grid
    else:
        x, y = smapmeshgrid2(sunpymap)
    r = sunpymap.rsun_obs
    rr = np.sqrt(x * x + y * y)
    maskout = rr > r
    try:
        sunpymap.data[maskout] = sunpymap.data[maskout] * np.exp(factor * (rr[maskout] / r - 1))
    except:
        data = sunpymap.data.copy()
        data[maskout] = data[maskout] * np.exp(factor * (rr[maskout] / r - 1))
        sunpymap = smap.Map(data, sunpymap.meta)
    return sunpymap


def regridimage(values, x, y, grid=None, resize=[1.0, 1.0]):
    '''
    re-grid the data on a regular grid with uneven grid spacing to an uniform grid
    :param values: The image data on the regular grid
    :param x: the points defining the regular grid in x
    :param y: the points defining the regular grid in y
    :param grid: new uniform mesh grid [gridx,gridy]
    :param resize: list of re-size ratio factors of x and y. if resize is not [1.0,1.0], grid is neglected.
    :return: re-gridded image
    '''
    from scipy.interpolate import RegularGridInterpolator
    ny, nx = values.shape
    if grid and resize == [1.0, 1.0]:
        gridx, gridy = grid
    else:
        gridx, gridy = np.meshgrid(np.linspace(x[0], x[-1], nx * resize[0]), np.linspace(y[0], y[-1], ny * resize[1]))
    ny, nx = gridx.shape
    rgi = RegularGridInterpolator(points=(y, x), values=values, bounds_error=False)
    datanew = rgi(np.stack(np.stack((gridy.ravel(), gridx.ravel()), axis=-1))).reshape(ny, nx)
    if grid:
        return datanew
    else:
        return [datanew, gridx, gridy]


def regridspec(spec, x, y, nxmax=None, nymax=None, interp=False):
    '''
    :param spec: ndarray of float or complex, shape (npol,nbl,nf,nt) Data values.
    :param x: Data point x coordinates.
    :param y: Data point y coordinates.
    :param nxmax:
    :param nymax:
    :return:
    '''

    npol, nbl, nf, nt = spec.shape
    if interp:
        if nxmax:
            if nt > nxmax:
                nt = nxmax
        if nymax:
            if nf > nymax:
                nf = nymax
        specnew = np.zeros((npol, nbl, nf, nt))
        tt = np.linspace(x[0], x[-1], nt)
        ff = np.linspace(y[0], y[-1], nf)
        grid_x, grid_y = np.meshgrid(tt, ff)
        for p in range(npol):
            for b in range(nbl):
                specnew[p, b, :, :] = regridimage(spec[p, b, :, :], x, y, grid=[grid_x, grid_y])
    else:
        xstep, ystep = 1, 1
        if nxmax:
            if nt > nxmax:
                import math
                xstep = int(math.ceil(float(nt) / nxmax))
        if nymax:
            if nf > nymax:
                ystep = int(float(nf) / nymax)
        specnew = spec[:, :, ::ystep, ::xstep]
    return [specnew, xstep, ystep]


def get_contour_data(X, Y, Z, levels=[0.5, 0.7, 0.9]):
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm
    from bokeh.models import (ColumnDataSource)
    try:
        cs = plt.contour(X, Y, Z, levels=(np.array(levels) * np.nanmax(Z)).tolist(), cmap=cm.Greys_r)
        # dx = X[0,1]-X[0,0]
        # dy = Y[1,0]-Y[0,0]
        xs = []
        ys = []
        xt = []
        yt = []
        col = []
        text = []
        isolevelid = 0
        for isolevel in cs.collections:
            # isocol = isolevel.get_color()[0]
            # thecol = 3 * [None]
            theiso = '{:.0f}%'.format(cs.get_array()[isolevelid] / Z.max() * 100)
            isolevelid += 1
            # for i in range(3):
            # thecol[i] = int(255 * isocol[i])
            thecol = '#%02x%02x%02x' % (220, 220, 220)
            # thecol = '#03FFF9'

            for path in isolevel.get_paths():
                v = path.vertices
                # x = v[:, 0]+dx
                # y = v[:, 1]+dy
                x = v[:, 0]
                y = v[:, 1]
                xs.append(x.tolist())
                ys.append(y.tolist())
                xt.append(x[len(x) / 2])
                yt.append(y[len(y) / 2])
                text.append(theiso)
                col.append(thecol)
        for coll in cs.collections:
            coll.remove()

        source = ColumnDataSource(data={'xs': xs, 'ys': ys, 'line_color': col, 'xt': xt, 'yt': yt, 'text': text})
    except:
        source = ColumnDataSource(data={'xs': [], 'ys': [], 'line_color': [], 'xt': [], 'yt': [], 'text': []})

    return source


# def twoD_Gaussian((x, y), amplitude, xo, yo, sigma_x, sigma_y, theta, offset):
#     xo = float(xo)
#     yo = float(yo)
#     a = (np.cos(theta) ** 2) / (2 * sigma_x ** 2) + (np.sin(theta) ** 2) / (2 * sigma_y ** 2)
#     b = -(np.sin(2 * theta)) / (4 * sigma_x ** 2) + (np.sin(2 * theta)) / (4 * sigma_y ** 2)
#     c = (np.sin(theta) ** 2) / (2 * sigma_x ** 2) + (np.cos(theta) ** 2) / (2 * sigma_y ** 2)
#     g = offset + amplitude * np.exp(- (a * ((x - xo) ** 2) + 2 * b * (x - xo) * (y - yo) + c * ((y - yo) ** 2)))
#     return g.ravel()


def c_correlate(a, v, returnx=False):
    a = (a - np.mean(a)) / (np.std(a) * len(a))
    v = (v - np.mean(v)) / np.std(v)
    if returnx:
        return np.arange(len(a)) - np.floor(len(a) / 2.0), np.correlate(a, v, mode='same')
    else:
        return np.correlate(a, v, mode='same')


def c_correlateX(a, v, returnx=False, returnav=False, s=0):
    '''

    :param a:
    :param v: a and v can be a dict in following format {'x':[],'y':[]}. The length of a and v can be different.
    :param returnx:
    :return:
    '''
    from scipy.interpolate import splev, splrep
    import numpy.ma as ma

    if isinstance(a, dict):
        max_a = np.nanmax(a['x'])
        min_a = np.nanmin(a['x'])
        max_v = np.nanmax(v['x'])
        min_v = np.nanmin(v['x'])
        max_ = min(max_a, max_v)
        min_ = max(min_a, min_v)
        if not max_ > min_:
            print('the x ranges of a and v have no overlap.')
            return None
        a_x = ma.masked_outside(a['x'].copy(), min_, max_)
        if isinstance(a['y'], np.ma.core.MaskedArray):
            a_y = ma.masked_array(a['y'].copy(), a['y'].mask | a_x.mask)
            a_x = ma.masked_array(a_x, a_y.mask)
        else:
            a_y = ma.masked_array(a['y'].copy(), a_x.mask)
        v_x = ma.masked_outside(v['x'].copy(), min_, max_)
        if isinstance(v['y'], np.ma.core.MaskedArray):
            v_y = ma.masked_array(v['y'].copy(), v['y'].mask | v_x.mask)
            v_x = ma.masked_array(v_x, v_y.mask)
        else:
            v_y = ma.masked_array(v['y'].copy(), v_x.mask)

        dx_a = np.abs(np.nanmean(np.diff(a_x)))
        dx_v = np.abs(np.nanmean(np.diff(v_x)))
        if dx_a >= dx_v:
            v_ = v_y.compressed()
            x_ = v_x.compressed()
            tck = splrep(a_x.compressed(), a_y.compressed(), s=s)
            ys = splev(x_, tck)
            a_ = ys

        elif dx_a < dx_v:
            a_ = a_y.compressed()
            x_ = a_x.compressed()
            tck = splrep(v_x.compressed(), v_y.compressed(), s=s)
            ys = splev(x_, tck)
            v_ = ys

    else:
        a_ = a.copy()
        v_ = v.copy()
        x_ = None
    a_ = (a_ - np.nanmean(a_)) / (np.nanstd(a_) * len(a_))
    v_ = (v_ - np.nanmean(v_)) / np.nanstd(v_)
    if returnx:
        if x_ is None:
            return np.arange(len(a_)) - np.floor(len(a_) / 2.0), np.correlate(a_, v_, mode='same')
        else:
            return (np.arange(len(a_)) - np.floor(len(a_) / 2.0)) * np.nanmean(np.diff(x_)), np.correlate(a_, v_,
                                                                                                          mode='same'), x_, a_, v_
    else:
        return np.correlate(a_, v_, mode='same')


def XCorrMap(z, x, y, doxscale=True):
    '''
    get the cross correlation map along y axis
    :param z: data
    :param x: x axis
    :param y: y axis
    :return:
    '''
    from tqdm import tqdm
    from scipy.interpolate import splev, splrep
    if doxscale:
        xfit = np.linspace(x[0], x[-1], 10 * len(x) + 1)
        zfit = np.zeros((len(y), len(xfit)))
        for yidx1, yq in enumerate(y):
            xx = x
            yy = z[yidx1, :]
            s = len(yy)  # (len(yy) - np.sqrt(2 * len(yy)))*2
            tck = splrep(xx, yy, s=s)
            ys = splev(xfit, tck)
            zfit[yidx1, :] = ys
    else:
        xfit = x
        zfit = z
    ny, nxfit = zfit.shape
    ccpeak = np.empty((ny - 1, ny - 1))
    ccpeak[:] = np.nan
    ccmax = ccpeak.copy()
    ya = ccpeak.copy()
    yv = ccpeak.copy()
    yidxa = ccpeak.copy()
    yidxv = ccpeak.copy()
    for idx1 in tqdm(range(1, ny)):
        for idx2 in range(0, idx1):
            lightcurve1 = zfit[idx1, :]
            lightcurve2 = zfit[idx2, :]
            ccval = c_correlate(lightcurve1, lightcurve2)
            if sum(lightcurve1) != 0 and sum(lightcurve2) != 0:
                cmax = np.amax(ccval)
                cpeak = np.argmax(ccval) - (nxfit - 1) / 2
            else:
                cmax = 0
                cpeak = 0
            ccmax[idx2, idx1 - 1] = cmax
            ccpeak[idx2, idx1 - 1] = cpeak
            ya[idx2, idx1 - 1] = y[idx1 - 1]
            yv[idx2, idx1 - 1] = y[idx2]
            yidxa[idx2, idx1 - 1] = idx1 - 1
            yidxv[idx2, idx1 - 1] = idx2
            if idx1 - 1 != idx2:
                ccmax[idx1 - 1, idx2] = cmax
                ccpeak[idx1 - 1, idx2] = cpeak
                ya[idx1 - 1, idx2] = y[idx2]
                yv[idx1 - 1, idx2] = y[idx1 - 1]
                yidxa[idx1 - 1, idx2] = idx2
                yidxv[idx1 - 1, idx2] = idx1 - 1

    return {'zfit': zfit, 'ccmax': ccmax, 'ccpeak': ccpeak, 'x': x, 'nx': len(x), 'xfit': xfit, 'nxfit': nxfit, 'y': y,
            'ny': ny, 'yv': yv, 'ya': ya,
            'yidxv': yidxv, 'yidxa': yidxa}


class ButtonsPlayCTRL():
    '''
    Produce A play/stop button widget for bokeh plot

    '''
    __slots__ = ['buttons']

    def __init__(self, plot_width=None, *args, **kwargs):
        from bokeh.models import Button
        BUT_first = Button(label='<<', width=plot_width, button_type='primary')
        BUT_prev = Button(label='|<', width=plot_width, button_type='warning')
        BUT_play = Button(label='>', width=plot_width, button_type='success')
        BUT_next = Button(label='>|', width=plot_width, button_type='warning')
        BUT_last = Button(label='>>', width=plot_width, button_type='primary')
        self.buttons = [BUT_first, BUT_prev, BUT_play, BUT_next, BUT_last]

        # class FileDialog():  #     '''  #     produce a file dialog button widget for bokeh plot  #     '''  #     import Tkinter  #     import tkFileDialog  #  #     def __init__(self, plot_width=30, labels={'dir':'...','open':'open','save':'save'}, *args,  #                  **kwargs):  #         from bokeh.models import Button  #         buttons = {}  #         for k,v in labels.items():  #             buttons[k] = Button(label=v, width=plot_width)  #         self.buttons = buttons  #  #     def askdirectory(self):  #         tkRoot = Tkinter.Tk()  #         tkRoot.withdraw()  # Close the root window  #         in_path = tkFileDialog.askdirectory()  #         tkRoot.destroy()  #         if in_path:  #             return in_path  #  #     def askopenfilename(self):  #         tkRoot = Tkinter.Tk()  #         tkRoot.withdraw()  # Close the root window  #         in_path = tkFileDialog.askopenfilename()  #         tkRoot.destroy()  #         if in_path:  #             return in_path  #  #     def asksaveasfilename(self):  #         tkRoot = Tkinter.Tk()  #         tkRoot.withdraw()  # Close the root window  #         in_path = tkFileDialog.asksaveasfilename()  #         tkRoot.destroy()  #         if in_path:  #             return in_path
