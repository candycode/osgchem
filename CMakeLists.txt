cmake_minimum_required(VERSION 3.5)
project(osgchem)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
PROJECT( osgob )

if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    include_directories( /usr/local/openbabel/include/openbabel /usr/local/osg/include)
    link_directories( /usr/local/openbabel/lib /usr/local/osg/lib)
endif()

#IF( WIN32 )
#    IF( MINGW OR MSYS )
#        INCLUDE_DIRECTORIES( C:/mingw/devlib/release/OpenSceneGraph/include )
#        LINK_DIRECTORIES( C:/mingw/devlib/release/OpenSceneGraph/lib )
#    ELSE( MINGW OR MSYS )
#        INCLUDE_DIRECTORIES( C:/vclibs/OpenSceneGraph-svn-latest/include )
#        LINK_DIRECTORIES( C:/vclibs/OpenSceneGraph-svn-latest/lib )
#    ENDIF( MINGW OR MSYS )
#ENDIF( WIN32 )
#
#IF( WIN32 )
#    IF( MINGW OR MSYS )
#        SET( OPENBABEL_INCLUDE_DIR C:/mingw/openbabel-2.1.0/include/openbabel-2.0 ) #<<<<<!!!!
#        SET( OPENBABEL_LINK_DIR C:/mingw/openbabel-2.1.0/lib )
#        #SET( OPENBABEL_INCLUDE_DIR C:/ugo/mingw/include/openbabel-2.0 ) #<<<<<!!!!
#        SET( OBLIB openbabel )
#    ELSE( MINGW OR MSYS )
#        SET( OPENBABEL_INCLUDE_DIR C:/vclibs/openbabel-2.1.1/include ) #<<<<<!!!!
#        SET( OPENBABEL_LINK_DIR C:/vclibs/openbabel-2.1.1/lib )
#        SET( OBLIB openbabel )
#    ENDIF( MINGW OR MSYS )
#ELSE ( WIN32 )
#    SET( OPENBABEL_INCLUDE_DIR /usr/local/include/openbabel-2.0 )
#    SET( OBLIB openbabel )
#ENDIF( WIN32 )

#IF( UNIX )
#    SET( OPENBABEL_INCLUDE_DIR /usr/local/openbabel/include )
#    SET( OPENBABEL_LINK_DIR /usr/local/openbabel/lib )
#    SET( OBLIB openbabel )
#    INCLUDE_DIRECTORIES( /usr/local/osg/include )
#    LINK_DIRECTORIES( /usr/local/osg/lib64 )
#ENDIF( UNIX )

include_directories( ${OPENBABEL_INCLUDE_DIR} )
link_directories( ${OPENBABEL_LINK_DIR} )

#IF( MINGW OR MSYS )
#    SET( PLUGIN mingw_osgdb_molecule )
#ELSE( MINGW OR MSYS )
#    SET( PLUGIN osgdb_molecule )
#ENDIF( MINGW OR MSYS )

set( PLUGIN_SRCS src/ReaderWriterMolecule.cpp src/osgob.cpp src/ElementTable.cpp src/ElementTable.h src/default_atom_colors.cpp )
set( SRCS  src/osgob.cpp src/ElementTable.cpp src/default_atom_colors.cpp src/ElementTable.h )

add_library( osgchem SHARED ${PLUGIN_SRCS} )
#add_executable( osgob ${SRCS} )

target_link_libraries( osgchem
        openbabel
        OpenThreads osg osgDB osgUtil osgViewer osgText osgGA osgSim )