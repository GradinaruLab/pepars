from distutils.core import setup

setup(
    name="pepars",
    version="0.3",
    packages=[
        "pepars",
        "pepars.analysis",
        "pepars.plotting",
        "pepars.utils",
        "pepars.alignment",
        "pepars.fileio"
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
