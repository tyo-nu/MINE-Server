from setuptools import setup
import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(name='mine-server',
      version='2.0.0',
      description='MINE API',
      long_description=long_description,
      long_description_content_type="text/markdown",
      author='Jonathan Strutz',
      author_email='jonstrutz11@gmail.com',
      license='MIT',
      packages=setuptools.find_packages(),
      install_requires=['pymongo', 'sphinxcontrib-httpdomain'],
      extras_require={},
      classifiers=[
          'Development Status :: 5 - Production/Stable',
          'Intended Audience :: Science/Research',
          'License :: OSI Approved :: MIT License',
          'Topic :: Scientific/Engineering :: Bio-Informatics',
          'Topic :: Scientific/Engineering :: Chemistry',
          'Programming Language :: Python :: 3',
      ],
      )
