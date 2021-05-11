import os
import setuptools

# Get the package dependencies
project_folder = os.path.dirname(os.path.realpath(__file__))
requirements_file = os.path.join(project_folder,'requirements.txt')
install_requires = [] # Examples: ["gunicorn", "docutils>=0.3", "lxml==0.5a7"]
if os.path.isfile(requirements_file):
    with open(requirements_file) as f:
        install_requires = f.read().splitlines()

# Get the package description from README.md file
with open("README.md", "r") as fh:
    long_description = fh.read()

# Configure setup
setuptools.setup(
    name="hes_off",
    version="1.0.0",
    author="Roberto Agromayor",
    author_email="rober.agro@gmail.com",
    description="A package for the simulation of hybrid energy systems for off-shore oil and gas installations.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/RoberAgro/HES-OFF",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
    ],
    install_requires=install_requires,
    python_requires='>=3',
)