# anonymize_bam.py - de-identify sequencing data

import setuptools

# Set __version__ for the project.
exec(open("./BAMboozle/version.py").read())

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="BAMboozle",
    version=__version__,
    author="GM Franceschini, E Frigoli",
    author_email="franceschini.gianmarco@gmail.com",
    description="remove genetic variation from bisulfite sequencing data",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/gmfranceschini/dataprivacy_bisulfite",
    packages=setuptools.find_packages(),
    install_requires=["pysam>=0.14.0"],
    entry_points={
        "console_scripts": [
            "BAMboozle = BAMboozle.BAMboozle:main",
        ],
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
    ],
    python_requires=">=3.6",
    license="GPL-3.0-or-later",
)
