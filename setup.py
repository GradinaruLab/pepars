from distutils.core import setup

setup(
    name="PE_deep_seq",
    version="0.1",
    packages=[
        "PE_deep_seq",
        "PE_deep_seq.analysis",
        "PE_deep_seq.plotting",
        "PE_deep_seq.utils"
    ],
    install_requires=["numpy", "plotly"]
)
