from setuptools import setup, find_packages


setup(
    name="ml4t",
    version="0.1.1",
    description="ml4t is an open-source code for automatically train Machine Learning Forcefields for twisted bilayer materials",
    author="jiaxuan Liu",
    python_requires=">=3.7",
    packages=find_packages(include=["ml4t", "ml4t.*"]),
    entry_points={
        # make the scripts available as command line scripts
        "console_scripts": [
            "ml4tTrain = ml4t.script.ml4tTrain:main",
            "ml4tTest = ml4t.script.ml4tTest:main"
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