from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="CDSmsa",
    version="0.0.2",
    author="ToHanwei",
    author_email="whanwei@foxmail.com",
    description="Coding sequence MSA with optional HMM profile building",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/ToHanwei/CDSmsa",
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.6",
    install_requires=[
        "biopython",
        "mafft",
        "hmmer",
    ],
    entry_points={
        'console_scripts': [
            'single_msa=CDSmsa.single_msa:main',
            'batch_msa=CDSmsa.batch_msa:main',
        ],
    },
) 