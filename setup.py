from setuptools import setup
from setuptools import find_packages

with open("requirements.txt") as f:
    content = f.readlines()

requirements = [x.strip() for x in content]

setup(
    name="DDapp",
    version="0.0.1",
    author="Re",
    description="AI models and Tools for KFMap",
    packages=find_packages(),
    install_requires=requirements
)
