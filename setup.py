from setuptools import setup, find_packages


setup(
    name="DPmoire",
    version="1.0.0",
    description="DPmoire is an open-source code for automatically train Machine Learning Forcefields for twisted bilayer materials",
    author="Jiaxuan Liu",
    python_requires=">=3.10",
    packages=find_packages(include=["DPmoire", "DPmoire.*"]),
    entry_points={
        # make the scripts available as command line scripts
        "console_scripts": [
            "DPmoireTrain = DPmoire.script.DPmoireTrain:main",
            "DPmoireTest = DPmoire.script.DPmoireTest:main"
        ]
    },
    install_requires=[
        "numpy==1.26.3",
        "ase",
        "tqdm",
        "pyyaml",
        "spglib"
    ],
    zip_safe=True,
)
