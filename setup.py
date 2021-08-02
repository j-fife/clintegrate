from setuptools import setup, find_packages


setup(
    name='clintegrate',
    version='1.0.1',
    description='Package for clintegrate framework',
    author='James Fife',
    author_email='jamesdavidfife@gmail.com',
    url='https://github.com/j-fife/clintegrate',
    packages=find_packages(exclude=['test']),
    py_modules = ["clintegrate"]
)
