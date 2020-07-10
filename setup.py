import setuptools

setuptools.setup(
    name="qiime2utils",
    version="v0.0.3",
    author="Vini Salazar",
    author_email="viniws@gmail.com",
    description="qiime2utils - Utility scripts for Qiime 2",
    long_description="qiime2utils is a command-line tool for manipulating Qiime 2 data.",
    url="https://github.com/vinisalazar/qiime2utils",
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3.6",
    ],
    packages=setuptools.find_packages(),
    scripts=["qiime2utils/qiime2utils"],
    include_package_data=True,
    python_requires=">=3.6",
    install_requires=("pytest", "pandas", "biopython"),
)
