Overview
========
This is the NRPSDesigner software, developed by the iGEM Team Heidelberg during iGEM 2013. It is capable of designing a non-ribosomal peptide synthetase from scratch given a desired peptide which needs to be produced.

Directory structure
-------------------
* __Django/__: contains the web interface, written in Python using the Django framework
* __NRPS library/__: contains the C++ code for the NRPSDesigner executable
* __sbspks__/: contains the scraping scripts that were used to scrape the SBSPKS database. These require Scrapy 0.16 as well as all packages needed for Django
* __nrps_designer.sql__: Database dump. Note that all users' passwords and email addresses have been removed and the users have been set to inactive.
* __nrpsdesigner__: Executable compiled from the C++ code. This was compiled on Debian Wheezy with GCC 4.8.1. You may need to recompile to match your system's library versions.

Django stuff
=============
Because we all LOVE packages (AND LOTS OF THEM), here is the list of stuff you will need to run everything. The Python packages can be installed with pip by running
`pip install package`
If you are on Mac OS X or Linux, you need to run this as root:
`sudo pip install package`

Main packages
--------------
* [__Python__](http://www.python.org) (2.7)
* [__Django__](http://www.djangoproject.com) (>= 1.5)
* [__MySQL__](http://www.mysql.com) (>= 5.5.3)
* [__MySQL - Python connector__](http://sourceforge.net/projects/mysql-python/)

Django packages
---------------
* [__django-registration__](https://bitbucket.org/ubernostrum/django-registration/): Handles user login etc.
* [__django-annoying__](http://skorokithakis.github.io/django-annoying/): Gibthon uses it...
* [__south__](http://south.aeracode.org/): handles database schema changes
* [__django-celery__](http://celery.github.io/django-celery/): integration of Celery into Django

Other stuff that is loaded by some of our python functions:
-----------------------------------------------------------
* [__Requests__](http://docs.python-requests.org/en/latest/):Great library to do simple http requests!
* [__RabbitMQ__](http://www.rabbitmq.com): Message broker required by Celery
* [__Biopython__](http://www.biopython.org) (>= 1.62): Just like bioperl but not in a crappy language!
* [__xhtml2pdf__](http://www.xhtml2pdf.com): Also necessary for Gibthon..
* [__BeautifulSoup__](http://www.crummy.com/software/BeautifulSoup/): Beautiful soup!
* [__celery__](http://celeryproject.org/): asynchronous scheduling framework
* [__kombu__](https://github.com/celery/kombu): database-based message passing framework
* [__python-openbabel__](http://www.openbabel.org): (>= 2.3.2) Python bindings to openbabel for 2D structure generation
* [__UNAFold__ and __MFold utils__](http://dinamelt.rit.albany.edu/download.php): required by Gibthon
* [__HMMER__](http://hmmer.janelia.org/) (= 3.0): Hidden Markov Model based sequence analysis software used by antiSMASH
* [__MUSCLE__](http://www.drive5.com/muscle/): multiple sequence alignment used by antiSMASH
* [__Clustal Omega__](http://www.clustal.org): multiple sequence alignment used by NRPSDesigner
* [__libSBOLpy__](https://github.com/SynBioDex/libSBOLpy): Python bindings to the libSBOLc library

C++ stuff
=============
* [__CMake__](http://www.cmake.org): build system used by NRPSDesigner
* [__MySQL C++ connector__](http://www.mysql.com): C++ library for MySQL access
* [__libXML__](http://www.xmlsoft.org): C library for XML handling
* [__libcurl__](http://curl.haxx.se): C library for HTTP
* [__Boost.program_options__](http://www.boost.org): C++ library for command-line option parsing
* [__libSBOLc__](https://github.com/SynBioDex/libSBOLc) (fixed version from https://github.com/ilia-kats/libSBOLc required): C library for SBOL IO
* [__GCC__](http://gcc.gnu.org) (>= 4.8) or a comparable compiler with C++11 support

Compilation
-----------
###*ix
For an out-of-tree build, create a directory named 'build'. Chdir into the build directory and run
`cmake path_to_NRPSDesigner_source`
This will locate all required libraries and create a Makefile. Run `make` to compile, `make install` to install.

####MacOSX
On MacOSX, both GCC 4.2 and Clang are available. Although GCC 4.2 is ancient and does not support C++11, Clang by default still links against the GCC C++ standard library. To successfully compile NRPSDesigner, Clang needs to be told to link against its own C++ standard library. This can be achieved by adding `-DCMAKE_CXX_FLAGS=-stdlib=libc++` to the CMake command-line or, if using the CMake GUI, adding `-stdlib=libc++`to `CMAKE_CXX_FLAGS`. Note that all C++ libraries the NRPSDesigner depends on, i.e. Boost.program_options and the MySQL C++ connector, need to also link against libc++.

###Windows
Launch the CMake graphical user interface, select the source and build directories (create a build directory for an out-of-tree build) and click `Configure`. CMake will now locate all required libraries. If it is unable to find some libraries, you will need to input the paths manually. Run `Generate` to generate a makefile in the format of your choice (we recommend using MinGW64 from http://sourceforge.net/projects/mingwbuilds/ inside an MSYS environment and use the MSYS Makefile generator). Note that no attempts to compile NRPSDesigner on Windows were made, so you might run into trouble.
