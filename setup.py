from setuptools import setup

setup(
    name="pepars",
    version="0.4.1",
    packages=[
        "pepars",
        "pepars.analysis",
        "pepars.plotting",
        "pepars.utils",
        "pepars.alignment",
        "pepars.fileio",
        "pepars.simulation"
    ],
    install_requires=[
        "numpy",
        "plotly",
        "pandas",
        "matplotlib",
        "biopython",
        "seaborn",
        "h5py",
        "scipy",
        "statsmodels"
    ]
)
