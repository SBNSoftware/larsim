# project name
project(geotk)

# ----------------------------------------------------------------------
file(GLOB GEOFILES Geo* )

art_make_library(SOURCE ${GEOFILES}
                 LIBRARY_NAME GeoTK
	         LIBRARIES 
                   ${ROOT_BASIC_LIB_LIST}
		 )

# ----------------------------------------------------------------------
# Install this project's headers, sources, library, and executables:

install_fhicl()
install_headers()
install_source()




