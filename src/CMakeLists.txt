PROJECT( osgob )

IF( APPLE )
  INCLUDE_DIRECTORIES( /usr/local/include )
  LINK_DIRECTORIES( /usr/local/lib )
ENDIF( APPLE )

IF( WIN32 )
  IF( MINGW OR MSYS )
    INCLUDE_DIRECTORIES( C:/mingw/devlib/release/OpenSceneGraph/include )
    LINK_DIRECTORIES( C:/mingw/devlib/release/OpenSceneGraph/lib )
  ELSE( MINGW OR MSYS )
    INCLUDE_DIRECTORIES( C:/vclibs/OpenSceneGraph-svn-latest/include )
    LINK_DIRECTORIES( C:/vclibs/OpenSceneGraph-svn-latest/lib )
  ENDIF( MINGW OR MSYS )
ENDIF( WIN32 )

IF( WIN32 )
  IF( MINGW OR MSYS )
    SET( OPENBABEL_INCLUDE_DIR C:/mingw/openbabel-2.1.0/include/openbabel-2.0 ) #<<<<<!!!! 
    SET( OPENBABEL_LINK_DIR C:/mingw/openbabel-2.1.0/lib )
    #SET( OPENBABEL_INCLUDE_DIR C:/ugo/mingw/include/openbabel-2.0 ) #<<<<<!!!! 
    SET( OBLIB openbabel ) 	
  ELSE( MINGW OR MSYS )
    SET( OPENBABEL_INCLUDE_DIR C:/vclibs/openbabel-2.1.1/include ) #<<<<<!!!! 
    SET( OPENBABEL_LINK_DIR C:/vclibs/openbabel-2.1.1/lib )
    SET( OBLIB openbabel ) 	
  ENDIF( MINGW OR MSYS )
ELSE ( WIN32 )
  SET( OPENBABEL_INCLUDE_DIR /usr/local/include/openbabel-2.0 )
  SET( OBLIB openbabel )
ENDIF( WIN32 )    

IF( UNIX )
  SET( OPENBABEL_INCLUDE_DIR /usr/local/openbabel/include )
  SET( OPENBABEL_LINK_DIR /usr/local/openbabel/lib )
  SET( OBLIB openbabel )
  INCLUDE_DIRECTORIES( /usr/local/osg/include )
  LINK_DIRECTORIES( /usr/local/osg/lib64 ) 
ENDIF( UNIX )

INCLUDE_DIRECTORIES( ${OPENBABEL_INCLUDE_DIR} )
LINK_DIRECTORIES( ${OPENBABEL_LINK_DIR} )

IF( MINGW OR MSYS )
  SET( PLUGIN mingw_osgdb_molecule )
ELSE( MINGW OR MSYS )
  SET( PLUGIN osgdb_molecule )
ENDIF( MINGW OR MSYS )

SET( PLUGIN_SRCS ReaderWriterMolecule.cpp osgob.cpp ElementTable.cpp ElementTable.h default_atom_colors.cpp )
SET( SRCS  osgob.cpp ElementTable.cpp default_atom_colors.cpp ElementTable.h )

#ADD_EXECUTABLE( osgob ${SRCS} )

ADD_LIBRARY( ${PLUGIN} SHARED ${PLUGIN_SRCS} )

TARGET_LINK_LIBRARIES( ${PLUGIN}
                       ${OBLIB}
                       OpenThreads osg osgDB osgUtil osgViewer osgText osgGA osgSim )
#TARGET_LINK_LIBRARIES( osgob
#                       ${OBLIB}
#                       OpenThreads osg osgDB osgUtil osgViewer osgText osgGA osgSim )                       
                       
                      
                       
