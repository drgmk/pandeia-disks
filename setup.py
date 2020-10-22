import setuptools

setuptools.setup(
    name="pandeia-disks",
    version="0.0.1",
    author="Grant M. Kennedy",
    author_email="g.kennedy@warwick.ac.uk",
    description="disk-related functions for pandeia",
    packages=['pandisk'],
    classifiers=['Programming Language :: Python :: 3'],
    install_requires = ['matplotlib','numpy'],
    python_requires='>=3.6',
)
