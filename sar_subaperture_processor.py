import os
import sys
import numpy as np
import rasterio
from PIL import Image, ImageDraw, ImageFont
import imageio
from tqdm import tqdm
from rasterio.transform import Affine
from pathlib import Path


class SARSubapertureProcessor:
    def __init__(self, slc_path, structure):
        """
        Initialize with an SLC file.

        Parameters:
        slc_path : str
            Path to the SLC file
        structure : str 
            Either 'IQ' for I/Q bands or 'AP' for amplitude bands 
        """
        self.slc_path = slc_path
        self.structure = structure.upper()
        self.slc = None
        self.transform = None
        self.crs = None
        self.subapertures = None

    
    def read_slc(self):
        """
        Read SLC data from file and return slc, transform, crs
        """

        with rasterio.open(self.slc_path) as src:
            if self.structure == 'IQ': 
                I = src.read(1).astype(np.float32)
                Q = src.read(2).astype(np.float32)
                self.slc= I + 1j*Q
        
            elif self.structure == 'AP':
                amp = src.read(1).astype(np.float32)
                phase = src.read(2).astype(np.float32)
                I = amp * np.cos(phase)
                Q = amp * np.sin(phase)
                self.slc= I + 1j*Q

            self.transform= src.transform
            self.crs = src.crs

        return self.slc, self.transform, self.crs


def generate_subapertures(self, win_frac, step_frac, filter='raised_cosine'):
    """
    Generate azimuth sub apertures using rolling windows FFT method.

    Parameters:
        win_frac : float
            Window size as fraction of bandwidth
        step_frac : float
            Step size as fraction of bandwidth
        filter : str
            filter used: 'rectangular', 'raised_cosine', or 'hamming'
    """
    self.win_frac = win_frac  # needed for downsampling

    az_fft = np.fft.fftshift(np.fft.fft(self.slc, axis=0), axes=0)
    n_az = az_fft.shape[0]
    freqs = np.linspace(-0.5, 0.5, n_az, endpoint=False)
    win_half = win_frac / 2
    starts = np.arange(-0.5, 0.5, step_frac)
    self.subapertures = []

    for s in starts:
        center = (s + win_half) % 1.0 - 0.5
        if (center - win_half) < -0.5 or (center + win_half) > 0.5:
            continue
        df = np.abs(((freqs - center + 0.5) % 1.0) - 0.5)
        w = np.zeros_like(freqs, dtype=np.float32)
        inside = df <= win_half

        if filter == 'rectangular':
            w[inside] = 1.0
        elif filter == 'raised_cosine':
            w[inside] = 0.5 * (1 + np.cos(np.pi * df[inside] / win_half))
        elif filter == 'hamming':
            w[inside] = 0.54 + 0.46 * np.cos(np.pi * df[inside] / win_half)
        else:
            raise ValueError(f"Unknown filter: {filter}. Choose 'rectangular', 'raised_cosine' or 'hamming'.")

        azw = az_fft * w[:, None]
        sub = np.fft.ifft(np.fft.ifftshift(azw, axes=0), axis=0)
        self.subapertures.append(sub)

    return self

    def save_subs_intensity(self, out_dir, prefix="sub", bit_depth=16, georeference=True, scale="log"):
        """
        Save Subapertures as GeoTiff files

        Parameters:
        out_dir : str
            Ouput directory.
        prefix : str
            Prefix for output filenames.
        bit_depth : int
            Either 8 or 16 bit output.
        scale : str
            Either 'linear' or 'log/db".
        georeference : bool
            Include georeferencing (by default only for 16bit images).
        """
        Path(out_dir).mkdir(parents=True, exist_ok=True)

        if bit_depth == 8:
            dtype= 'uint8'
            max_val = 255
            georeference = False
        elif bit_depth == 16:
            dtype = 'float32'
            max_val = 65535

        for i,sub in enumerate(self.subapertures):
            intensity = np.abs(sub) ** 2

            if scale.lower() in ['log', 'db']:
                intensity = np.log1p(intensity)
            #normalize intensity     need to add if statement.
            #intensity = intensity - intensity.min()
            #if intensity.max() > 0:
            #    intensity = intensity / intensity.max()
            #intensity = (intensity * max_val).astype(dtype)

            fname = f"{out_dir}/{prefix}_{i:03d}.tif"

            use_transform = self.transform if georeference and bit_depth == 16 else Affine.identity()
            use_crs = self.crs if georeference and bit_depth == 16 else None

            with rasterio.open(
                fname,
                "w",
                driver="GTiff",
                height=intensity.shape[0],
                width=intensity.shape[1],
                count=1,
                dtype=dtype,
                transform=use_transform,
                crs=use_crs
            ) as dst:
                dst.write(intensity, 1)
        return self
    
    def save_gif(self,output_path, scale='log', fps=5, loop=0):
        """ export sub aperture as a GIF.

        Parameters:
        output_path : str
            Full path for output GIF file.
        scale : str
            Either 'linear' or 'log/db".
        fps : int
            Frames per second
        loop : int
            Number of loops (O = infinite)
        """
        Path(output_path).parent.mkdir(parents=True, exist_ok=True)

        duration= int(1000/fps)
        total_duration=len(self.subapertures) / fps

        gif_frames= []
        for sub in self.subapertures:
            intensity = np.abs(sub) **2

            if scale.lower() in ['log', 'db']:
                intensity = np.log1p(intensity)
            #Normalize intensity
            intensity = intensity - intensity.min()
            if intensity.max() > 0:
                intensity = intensity / intensity.max()
            intensity = (intensity * 255).astype('uint8')

            gif_frames.append(intensity)

        pil_frames = [Image.fromarray(frame, mode='L') for frame in gif_frames]
        
        #save as GIF
        pil_frames[0].save(
            output_path,
            save_all=True,
            append_images=pil_frames[1:],
            duration=duration,
            loop=loop
        )

    def save_subs_complex(self, out_dir,prefix='sub'):
        """
        Save Subapertures as complex GeoTiffs

        Parameters:
        output_path : str
            Full path for output GIF file.
        prefix : str
            Prefix for output filenames.
        """

        for i, sub in enumerate(self.subapertures):
            I= np.real(sub).astype('float32')
            Q= np.imag(sub).astype('float32')

            fname = f"{out_dir}/{prefix}_{i}.tif"

            use_transform = self.transform
            use_crs = self.crs

            with raserio.open(
                fname,
                "w",
                driver="GTiff",
                height=I.shape[0],
                width=I.shape[1],
                count=2,
                dtype="float32",
                transform= use_transform,
                crs=use_crs
            ) as dst:
                dst.write(I, 1)
                dst.write(Q, 2)
    
    def save_subs_intensity_w_downsample(self, out_dir, 
    prefix="sub", bit_depth=16, georeference=True, 
    scale="log",downsample_coef=None,downsample_range=0):

        """
        Save Subapertures as GeoTiff files with downsampling

        Parameters:
        out_dir : str
            Ouput directory.
        prefix : str
            Prefix for output filenames.
        bit_depth : int
            Either 8 or 16 bit output.
        scale : str
            Either 'linear' or 'log/db".
        georeference : bool
            Include georeferencing (by default only for 16bit images).
        downsample_coef : float
            Fraction of azimuth bandwidth used [0-1]

        downsample_range : bool
            If 0 => downsample range to match azimuth.
            Default 0 => non-square pixels (but physically correct).

        """
        Path(out_dir).mkdir(parents=True, exist_ok=True)
        az_ds = int(round(1 / downsample_coef))

        if downsample_coef is None:
            downsample_coef = self.win_frac

        if bit_depth == 8:
            dtype= 'uint8'
            max_val = 255
            georeference = False
        elif bit_depth == 16:
            dtype = 'float32'
            max_val = 65535

        for i,sub in enumerate(self.subapertures):
            intensity = np.abs(sub) ** 2

            if scale.lower() in ['log', 'db']:
                intensity = np.log1p(intensity)

            fname = f"{out_dir}/{prefix}_{i:03d}.tif"
            #Azimuth downsampling
            if az_ds > 1:
                intensity = intensity[::az_ds, :]
            #Range downsampling (optional)
            if downsample_range and az_ds > 1:
                intensity = intensity[:, ::az_ds]

#            if georeference and bit_depth == 16:
#                use_transform = self.transform * Affine.scale(1, az_ds)
#                use_crs = self.crs
#            else:
#                use_transform = Affine.identity()
#                use_crs = None

            if georeference and bit_depth == 16:
                sx = az_ds if downsample_range else 1
                sy = az_ds
                use_transform = self.transform * Affine.scale(sx, sy)
                use_crs = self.crs
            else:
                use_transform = Affine.identity()
                use_crs = None



            with rasterio.open(
                fname,
                "w",
                driver="GTiff",
                height=intensity.shape[0],
                width=intensity.shape[1],
                count=1,
                dtype=dtype,
                transform=use_transform,
                crs=use_crs
            ) as dst:
                dst.write(intensity, 1)
        return self



