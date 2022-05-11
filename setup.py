import setuptools

with open("README.md", "r") as f:
    long_description = f.read()

setuptools.setup(
    name="EEGToolkit", 
    version="2.0.0",
    author="Axel Giottonini, Noah Kleinschmidt, Kalvin Dobler",
    author_email="axel.giottonini@unifr.ch, noah.kleinschmidt@students.unibe.ch, kalvin.dobler@unifr.ch",
    description="A package for EEG data analysis of reaction time-delay experiments.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/AxelGiottonini/AP-EEG.git",
    packages=setuptools.find_packages(),
    entry_points ={
            "console_scripts": [
                "EEGToolkit = EEGToolkit.EEGData:main"
            ]
        },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Information Analysis",
        "Topic :: Scientific/Engineering :: Visualization"
    ],
    python_requires=">=3.6",
)