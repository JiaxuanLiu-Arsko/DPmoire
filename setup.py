from setuptools import setup, find_packages


setup(
    name="MLmoire",
    version="0.1.4",
    description="MLmoire is an open-source code for automatically train Machine Learning Forcefields for twisted bilayer materials",
    author="jiaxuan Liu",
    python_requires=">=3.7",
    packages=find_packages(include=["MLmoire", "MLmoire.*"]),
    entry_points={
        # make the scripts available as command line scripts
        "console_scripts": [
            "MLmoireTrain = MLmoire.script.MLmoire:main",
            "MLmoireTest = MLmoire.script.MLmoire:main"
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