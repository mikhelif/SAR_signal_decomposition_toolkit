from setuptools import setup

setup(
    name="sar-subaperture-processor",
    version="0.1.0",
    description="Lightweight Python tool for generating SAR azimuth sub-aperture images from SLC data",
    author="mourad_ikhelif",
    url="https://github.com/mikhelif/sar-subaperture-processor",
    python_requires=">=3.8",
    py_modules=["sar_subaperture_processor"],
    install_requires=[
        "numpy",
        "rasterio",
        "pillow",
        "imageio",
        "tqdm",
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
    ],
)
