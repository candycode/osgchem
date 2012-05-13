#ifndef OSGOB_H_
#define OSGOB_H_

#include <string>
#include <osg/PrimitiveSet>

namespace osg
{
	class Node;
	class Group;
}

enum MoleculeDisplayStyle { MOLDISPSTYLE_WIREFRAME,
	                        MOLDISPSTYLE_LIQUORICE,
							MOLDISPSTYLE_BALLANDSTICK,
							MOLDISPSTYLE_SPACEFILL };

enum SSAOMode { SSAO_SIMPLE,
                SSAO_TRACE_FRAG,
                SSAO_TRACE_FRAG_MRT,
                SSAO_TRACE_VERT,
                SSAO_TRACE_VERT_MRT };

osg::Geode* CreateSphereNode( double radius,
							  int slices,
							  int stacks,
							  osg::PrimitiveSet::Mode primitiveMode = osg::PrimitiveSet::QUADS,
							  double minTheta = -180.0,
							  double maxTheta = 180.0,
							  double minPhi   = -90.0,
							  double maxPhi   = 90.0 );

osg::Node* CreateCylinderNode( double radius, double height, bool open, int slices, int stacks );
void InitAtomColors( const char* fname );
void InitAtomGeometries( int slices = 16, int stacks = 8, float scaling = 1.0f );
void InitAtomDiscGeometries( int slices = 12 );
void InitAtomTextureGeometries();

osg::Group* CreateMoleculeGeometry2( const std::string& fname, const std::string& format, bool billboard = false, bool texture = false );
osg::Group* CreateMoleculeBackBoneGeometry( const std::string& fname, const std::string& format );
osg::Group* CreatePointSpriteMoleculeGeometry( const std::string& fname, const std::string& format );
osg::Group* CreateShaderPointSpriteMoleculeGeometry( const std::string& fname, const std::string& format );

//------------------------------------------------------------------------------
osg::Group* CreatePointRayTracedMoleculeGeometry( const std::string& fname,
						 				          const std::string& format,
												  MoleculeDisplayStyle ds,
												  float radScale = 1.0f );
osg::Group* CreatePointRayTracedMoleculeGeometry2( const std::string& fname,
						 				          const std::string& format,
												  MoleculeDisplayStyle ds,
												  float radScale = 1.0f );

namespace OpenBabel
{
	class OBMol;
}
osg::Group* CreateAnimatedPointRayTracedMoleculeGeometry2( OpenBabel::OBMol* mol,
														   MoleculeDisplayStyle ds,
														   float radScale = 1.0f,
														   bool ssao = false,
                                                           SSAOMode ssaoMode = SSAO_SIMPLE ); 

osg::Group* CreateRayTracedMoleculeGeometry3( const std::string& fname,
						 				          const std::string& format,
												  MoleculeDisplayStyle ds,
												  float radScale = 1.0f );
osg::Group* CreatePointSpriteMoleculeGeometry( const std::string& fname,
											   const std::string& format,
											   MoleculeDisplayStyle ds );
osg::Group* CreateBillboardMoleculeGeometry( const std::string& fname,
											 const std::string& format,
											 MoleculeDisplayStyle ds );
osg::Group* CreateRayTracedMoleculeGeometry( const std::string& fname,
						 				     const std::string& format,
											 MoleculeDisplayStyle ds,
											 float radScale = 1.0f );
osg::Group* CreateRayTracedMoleculeGeometry2( const std::string& fname,
						 				      const std::string& format,
											  MoleculeDisplayStyle ds,
											  float radScale = 1.0f );
osg::Group* CreateMoleculeGeometry( const std::string& fname,
								    const std::string& format,
									MoleculeDisplayStyle ds,
								    bool billboard = false,
									bool texture = false );




#endif
