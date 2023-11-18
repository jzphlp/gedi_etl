# Import required functions
from setuptools import setup, find_packages

# Call setup function
setup(
    author="<your-name>",
    description="A package for converting imperial lengths and weights.",
    name="gedi_etl",
    packages=find_packages(include=["gedi_etl", "gedi_etl.*"]),
    version="0.1.0",
    #  install_requires=['numpy>=1.10', 'pandas'],
)



# in folder and in env : install with pip install -e .
# pip freeze > requirement.txt @from raw enviroment to create this pkg 

# after add readme and license we can upload to 
# >source dist or wheel dist

p#ython setup.py sdist /bdist_wheel
#twine upload to pytest or pypi

# testing along the way