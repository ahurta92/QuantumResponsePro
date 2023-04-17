from setuptools import setup, find_packages

# Read the contents of the README file
with open("README.md", "r", encoding="utf-8") as f:
    long_description = f.read()

# Define the package metadata
setup(
    name="quantumresponsepro",
    version="0.1.0",
    author="Adrian Hurtado",
    author_email="ahurta92@gmail.com",
    description="A molecular analysis tool for computing and analyzing dynamic response properties using MADNESS and "
                "dalton quantum chemistry packages",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/yourusername/quantumresponsepro",
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Natural Language :: English",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
    ],
    keywords="quantum chemistry, molecular analysis, TDDFT, TDHF, MADNESS, dalton",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    python_requires=">=3.6, <4",
    install_requires=[
        "numpy",
        "scipy",
        "pandas",
        "jupyter",
        "matplotlib",
        "seaborn",
        'pytest',
        'pyarrow'
        # Add any other required dependencies here
    ],
    extras_require={
        "dev": [
            "pytest",
            # Add any other development dependencies here
        ],
    },
    entry_points={
        "console_scripts": [
            # Define any command-line scripts here, e.g.
            # "quantumresponsepro=quantumresponsepro.cli:main",
        ],
    },
)
