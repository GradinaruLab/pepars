from distutils.core import setup

setup(
    name="peseq",
    version="0.1",
    packages=[
        "peseq",
        "peseq.analysis",
        "peseq.plotting",
        "peseq.utils",
        "peseq.alignment",
        "peseq.fileio"
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
