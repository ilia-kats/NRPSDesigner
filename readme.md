Overview
========
This is the NRPSDesigner software, developed by the iGEM Team Heidelberg during iGEM 2013. It is capable of designing a non-ribosomal peptide synthetase from scratch given a desired peptide which needs to be produced.

Directory structure
-------------------
* __Django/__: contains the web interface, written in Python using the Django framework
* __NRPS library/__: contains the C++ code for the NRPSDesigner executable
* __sbspks__/: contains the scraping scripts that were used to scrape the SBSPKS database. These require Scrapy 0.16 as well as all packages needed for Django
* __nrps_designer.sql__: Database dump.

Django stuff
=============
Because we all LOVE packages (AND LOTS OF THEM), here is the list of stuff you will need to run everything. The Python packages can be installed with pip by running
`pip install package`
If you are on Mac OS X or Linux, you need to run this as root:
`sudo pip install package`

Main packages
--------------
* __Python__ (2.7)
* __Django__ (>= 1.5)
* __MySQL__
* __MySQL - Python connector__

Django packages
---------------
* __django-registration__: Handles user login etc.
* __django-annoying__: Gibthon uses it...
* __south__: handles database schema changes
* __django-celery__: integration of Celery into Django

Other stuff that is loaded by some of our python functions:
-----------------------------------------------------------
* __Requests__:Great library to do simple http requests!
* __Biopython__ (>= 1.61): Just like bioperl but not in a crappy language!
* __xhtml2pdf__: Also necessary for Gibthon..
* __BeautifulSoup__: Beautiful soup!
* __celery__: asynchronous scheduling framework
* __kombu__: database-based message passing framework
* __python-openbabel__: (>= 2.3.2) Python bindings to openbabel for 2D structure generation

C++ stuff
=============
* __CMake__: build system used by NRPSDesigner
* __MySQL C++ connector__: C++ library for MySQL access
* __libXML__: C library for XML handling
* __libcurl__: C library for HTTP
* __Boost.program_options__: C++ library for command-line option parsing
* __GCC__ (>= 4.8) or a comparable compiler with C++11 support

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
