'''
Originally written by Dr. Peijin Zhang @peijin94
'''

import subprocess
from typing import List, Optional, Union
import time
from datetime import datetime
from functools import wraps

def format_duration(seconds):
    """Format seconds into appropriate unit string."""
    if seconds < 60:
        return f"{seconds:.2f} seconds"
    minutes = seconds / 60
    if minutes < 60:
        return f"{minutes:.2f} minutes"
    hours = minutes / 60
    if hours < 24:
        return f"{hours:.2f} hours"
    days = hours / 24
    return f"{days:.2f} days"


def runtime_report(func):
    """Decorator to report runtime of a function and log completion time."""
    import time

    @wraps(func)
    def wrapper(*args, **kwargs):
        start = time.time()
        result = func(*args, **kwargs)
        end = time.time()
        duration = end - start
        now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        print(
            f"'{func.__name__}' completed at {now}; "
            f"runtime: {format_duration(duration)}"
        )
        return result
    return wrapper

class WSClean:
    def __init__(self, vis: str):
        """
        Initialize WSClean wrapper

        Parameters:
        -----------
        vis : str
            Input measurement set path
        """
        self.vis = vis
        self.params = {
            'size': [1024, 1024],  # default size
            'scale': "2.5asec",  # default pixel scale
            'weight_briggs': 0.0,  # default Briggs robust
            'niter': 0,  # default no cleaning
            'mgain': 1.0,  # default mgain
            'gain': 0.1,  # default gain
            'interval': [],  # default interval
            'data_column': "DATA",  # default data column
            'name': "wsclean",  # default output name
            'multiscale': False,  # default no multiscale
            'multiscale_gain': 0.1,  # default multiscale gain
            'multiscale_scales': [0, 5, 12.5],  # default multiscale scales
            'multiscale_scale_bias': 0.6,  # default multiscale scale bias. A lower bias will give more focus to larger scales. Default: 0.6
            'multiscale_shape': 'tapered-quadratic',  # default multiscale shape. Sets the shape function used during multi-scale clean. Either 'tapered-quadratic' (default) or 'gaussian'.
            'maxuvw_m': None,  # max uvw in meters
            'minuvw_m': None,  # min uvw in meters
            'maxuv_l': None,  # max uvw in lambda
            'minuv_l': None,  # min uvw in lambda
            'no_update_model': False,  # default update model
            'no_negative': False,  # default allow negative
            'beam_size': None,  # if not None, the beam size in arcsec
            'beam_fitting_size': 5,  # if not None, the beam fitting size in PSFs. Default: 5 PSFs
            'circular_beam': False,  # default no circular beam
            'theoretic_beam': False,  # Use the theoretical beam size instead of the fitted beam size
            'local_rms': False,  # Instead of using a single RMS for auto thresholding/masking, use a spatially varying RMS image
            'local_rms_strength': 1.0,  # default local rms strength
            'local_rms_window': 25,  # Size of window for creating the RMS background map, in number of PSFs. Default: 25 psfs
            'fits_mask': None,  # If not None, the FITS mask file to use
        }

    def setup(self, **kwargs):
        """
        Set WSClean parameters

        Parameters:
        -----------
        size : int or list, optional
            Image size in pixels [width, height] or single value for square
        scale : str, optional
            Pixel scale (e.g., "2.5asec")
        weight_briggs : float, optional
            Briggs robust parameter
        niter : int, optional
            Number of clean iterations
        mgain : float, optional
            Major cycle gain
        data_column : str, optional
            Data column to image
        name : str, optional
            Output name prefix
        multiscale : bool, optional
            Enable multiscale clean
        auto_mask : float, optional
            Auto-masking threshold
        auto_threshold : float, optional
            Auto threshold value
        local_rms : bool, optional
            Enable local RMS
        no_update_model : bool, optional
            Disable model data updates
        intervals_out : int, optional
            Number of time intervals
        spws : list or str, optional
            Spectral windows to image
        pol : str, optional
            Polarization to image
        no_negative : bool, optional
            Prevent negative components
        beam_size : float, optional
            Beam size in arcsec
        quiet : bool, optional
            Suppress output
        circular_beam : bool, optional
            Use circular beam
        """
        # Handle size parameter specially
        if 'size' in kwargs:
            size = kwargs['size']
            if isinstance(size, (int, float)):
                self.params['size'] = [int(size), int(size)]
            else:
                self.params['size'] = [int(size[0]), int(size[1])]

        # Update all other parameters
        for key, value in kwargs.items():
            if key != 'size':
                self.params[key] = value

        return self

    def build_command(self) -> str:
        """Build wsclean command"""
        cmd = ['wsclean']

        # Add basic parameters
        cmd.extend(['-size', str(self.params['size'][0]), str(self.params['size'][1])])
        cmd.extend(['-scale', self.params['scale']])
        cmd.extend(['-weight', 'briggs', str(self.params['weight_briggs'])])

        if self.params['niter'] > 0:
            cmd.extend(['-niter', str(self.params['niter'])])

        if len(self.params['interval'])==2:
            cmd.extend(['-interval', str(self.params['interval'][0]), str(self.params['interval'][1])])

        if self.params['multiscale']:
            cmd.append('-multiscale')
            cmd.extend(['-multiscale-gain', str(self.params['multiscale_gain'])])
            if len(self.params['multiscale_scales']) > 0:
                scales = ','.join(map(str, self.params['multiscale_scales']))
                cmd.extend(['-multiscale-scales', scales])
            cmd.extend(['-multiscale-scale-bias', str(self.params['multiscale_scale_bias'])])
            cmd.extend(['-multiscale-shape', self.params['multiscale_shape']])


        if self.params['gain'] != 0.1:
            cmd.extend(['-gain', str(self.params['gain'])])

        if self.params['mgain'] != 1.0:
            cmd.extend(['-mgain', str(self.params['mgain'])])

        if self.params['maxuvw_m'] is not None:
            cmd.extend(['-maxuvw-m', str(self.params['maxuvw_m'])])

        if self.params['minuvw_m'] is not None:
            cmd.extend(['-minuvw-m', str(self.params['minuvw_m'])])

        if self.params['maxuv_l'] is not None:
            cmd.extend(['-maxuv-l', str(self.params['maxuv_l'])])

        if self.params['minuv_l'] is not None:
            cmd.extend(['-minuv-l', str(self.params['minuv_l'])])

        if 'data_column' in self.params:
            if self.params['data_column'].lower().startswith('corrected'):
                self.params['data_column'] = 'CORRECTED_DATA'
            cmd.extend(['-data-column', self.params['data_column']])

        if 'pol' in self.params:
            cmd.extend(['-pol', self.params['pol']])

        if 'fits_mask' in self.params and self.params['fits_mask'] is not None:
            cmd.extend(['-fits-mask', self.params['fits_mask']])

        if 'auto_mask' in self.params:
            cmd.extend(['-auto-mask', str(self.params['auto_mask'])])

        if 'auto_threshold' in self.params:
            cmd.extend(['-auto-threshold', str(self.params['auto_threshold'])])

        if self.params['local_rms']:
            cmd.append('-local-rms')
            cmd.extend(['-local-rms-strength', str(self.params['local_rms_strength'])])
            cmd.extend(['-local-rms-window', str(self.params['local_rms_window'])])

        if self.params['no_update_model']:
            cmd.append('-no-update-model-required')

        if self.params['no_negative']:
            cmd.append('-no-negative')

        if self.params['beam_size'] is not None:
            cmd.extend(['-beam-size', str(self.params['beam_size'])])

        if self.params['beam_fitting_size'] is not None:
            cmd.extend(['-beam-fitting-size', str(self.params['beam_fitting_size'])])

        if self.params['theoretic_beam']:
            cmd.append('-theoretic-beam')

        if 'intervals_out' in self.params:
            cmd.extend(['-intervals-out', str(self.params['intervals_out'])])

        if self.params['quiet']:
            cmd.append('-quiet')

        if self.params['circular_beam']:
            cmd.append('-circular-beam')

        if 'spws' in self.params:
            spws = self.params['spws']
            if isinstance(spws, list):
                spws = ','.join(map(str, spws))
            if spws:
                cmd.extend(['-spws', spws])

        cmd.extend(['-name', self.params['name']])
        cmd.append(self.vis)

        return ' '.join(cmd)

    @runtime_report
    def run(self, dryrun: bool = False) -> int:
        """
        Run wsclean command

        Parameters:
        -----------
        dryrun : bool, optional
            If True, only print the command without executing

        Returns:
        --------
        int
            Return code from wsclean execution
        """
        cmd = self.build_command()

        if dryrun:
            print(f"Would run: {cmd}")
            return 0

        print(f"Running: {cmd}")
        process = subprocess.run(cmd, shell=True)
        return process.returncode


# Example usage:
if __name__ == "__main__":
    wsclean = WSClean("UDB20241215.ms")

    # Configure all parameters in one call
    wsclean.setup(
        size=1024,
        scale="2.5asec",
        weight_briggs=0.0,
        niter=1000,
        multiscale=True,
        multiscale_gain=0.1,
        multiscale_scales=[0, 5, 12.5],
        mgain=0.8,
        maxuvw_m=1500,
        minuvw_m=20,
        data_column="DATA",
        pol="xx",
        auto_mask=7,
        auto_threshold=2,
        local_rms=True,
        no_update_model=True,
        no_negative=True,
        spws=[4, 5, 6, 7, 8, 9],
        intervals_out=4,
        name="eovsa"
    )

    # Run WSClean (with dryrun=True to just print the command)
    wsclean.run(dryrun=True)
