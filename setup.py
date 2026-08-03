from setuptools import setup

setup(
    name="SAR_signal_decomposition_toolkit",
    version="0.1.0",
    description="Python toolkit for SAR signal decomposition of SLC data, including azimuth sub-apertures and range-frequency sub-bands.",
    author="mourad_ikhelif",
    url="https://github.com/mikhelif/sar-subaperture-processor",
    python_requires=">=3.8",
    py_modules=["SARDecomp"],
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
