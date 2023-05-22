from setuptools import setup, find_packages


setup(
    name="DPmoire",
    version="0.1.4",
    description="DPmoire is an open-source code for automatically train Machine Learning Forcefields for twisted bilayer materials",
    author="jiaxuan Liu",
    python_requires=">=3.7",
    packages=find_packages(include=["DPmoire", "DPmoire.*"]),
    entry_points={
        # make the scripts available as command line scripts
        "console_scripts": [
            "DPmoireTrain = DPmoire.script.DPmoireTrain:main",
            "DPmoireTest = DPmoire.script.DPmoireTest:main"
        ]
    },
    install_requires=[
        "numpy",
        "ase",
        "tqdm",
        "nequip",
        "pyyaml",
    ],
    zip_safe=True,
)