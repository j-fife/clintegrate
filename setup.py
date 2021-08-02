from setuptools import setup, find_packages

setup(
    name='clintegrate',
    version='1.0.6',
    description='Integrative Risk Predictors',
    author='James Fife',
    author_email='jamesdavidfife@gmail.com',
    url='https://github.com/j-fife/clintegrate',
    package_dir={"": "src"},
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    py_modules = ["clintegrate", "clintegrate.predictors", "clintegrate.predictors.risk"],
    packages=find_packages(where="src")
)
