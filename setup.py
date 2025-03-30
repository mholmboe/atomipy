from setuptools import setup, find_packages

setup(
    name="atomipy",
    version="0.2.0",
    packages=find_packages(),
    install_requires=[
        "numpy>=1.18.0",
    ],
    author="Michael Holmboe",
    author_email="michael.holmboe@umu.se",
    description="A Python toolbox for molecular structure analysis with periodic boundary conditions",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/mholmboe/atomipy",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Scientific/Engineering :: Physics",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Visualization", 
    ],
    python_requires=">=3.6",
    keywords="molecular dynamics, periodic boundary conditions, minerals, chemistry, physics",
)
