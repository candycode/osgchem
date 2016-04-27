#include <osg/Geode>
#include <osg/ShapeDrawable>
#include <osg/Material>
#include <osg/Texture2D>
#include <osg/Geometry>
#include <osg/CullFace>
#include <osg/MatrixTransform>
#include <osg/PositionAttitudeTransform>
#include <osg/BlendFunc>
#include <osg/ClearNode>
#include <osgUtil/Optimizer>
#include <osg/PointSprite>
#include <osg/Point>
#include <osg/Program>
#include <osg/Shader>
#include <osg/Uniform>
#include <osg/BoundingBox> 

#include <osgDB/ReadFile>

#include <osgSim/Impostor>
#include <osgSim/InsertImpostorsVisitor>
#include <osgGA/EventVisitor>
#include <osgUtil/CullVisitor>

#include <iostream>
#include <cassert>

#include <osgViewer/ViewerEventHandlers>
#include <osgGA/StateSetManipulator>

#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/obmolecformat.h>

#include <sstream>
#include <iterator>
#include <cmath>
#include <set>

using namespace OpenBabel;

#include "ElementTable.h"
#include "BSpline.h"
#include "default_atom_colors.h"
#include "osgob.h"

//----------------------------------------------------------------------------------
/// Creates a sphere as an indexed triangle or quad vertex array.
/// @param radius sphere radius
/// @param slices number of slices equal to the number of rotation steps around
///        the Y axis
/// @param stacks number of stacks along the Y axis
/// @param primitiveMode osg::PrimitiveSet::QUADS or osg::PrimitiveSet::TRIANGLES
/// @param minTheta start azimuth angle in degrees (rotation around Y axis)
/// @param maxTheta end azimuth angle in degrees (rotation around Y axis)
/// @param minPhi start elevation angle in degrees (rotation around X axis)
/// @param maxPhi end elevation angle in degrees (rotation around X axis)
/// @return osg::Geode containing and osg::Geometry representing the sphere
osg::Geode* CreateSphereNode( double radius, int slices, int stacks,
							  osg::PrimitiveSet::Mode primitiveMode,
							  double minTheta, double maxTheta,
							  double minPhi, double maxPhi )
{

#ifndef M_PI //not in the C/C++ standards but defined in MS VC up to version 2003
  #define M_PI 3.1415926535897932384626433832795
#endif
	assert( radius > 0.0 );
	assert( slices > 1 );
	assert( stacks > 1 );
	assert( minTheta <= maxTheta );
	assert( minPhi < maxPhi );
	assert( minPhi >= -90.0 && minPhi < 90.0 );
	assert( maxPhi > -90.0 && maxPhi <= 90.0 );
	minPhi = std::max( minPhi, -90.0 );
	const bool closed = !( int( maxTheta - minTheta ) % 360 ) ||  minTheta == maxTheta;
	const bool topCap = maxPhi == 90.0;
	const bool bottomCap = minPhi == -90.0;

	if( minTheta == maxTheta )
	{
		minTheta = 0.0;
		maxTheta = 2 * M_PI;
	}

	// convert to radians
	minTheta = M_PI * minTheta / 180.0;
	maxTheta = M_PI * maxTheta / 180.0;
	minPhi = M_PI * minPhi / 180.0;
	maxPhi = M_PI * maxPhi / 180.0;

	// geode
	osg::Geode* geode = new osg::Geode;
	osg::ref_ptr< osg::Geometry > geom = new osg::Geometry;

	// drawable (geometry)
	geode->addDrawable( geom.get() );

	// coordinates
	osg::Vec3Array* coords = new osg::Vec3Array;
	// normals
	osg::Vec3Array* normals = new osg::Vec3Array;

	geom->setVertexArray( coords );
	geom->setNormalArray( normals );

	geom->setNormalBinding( osg::Geometry::BIND_PER_VERTEX );

	// indices for corrdinates and normals
	osg::DrawElementsUInt* indices = new osg::DrawElementsUInt( primitiveMode, 0);
	geom->addPrimitiveSet( indices );

	// create points
	const int T = slices;
	const int P = stacks;
	const double r = radius;

	const double dT = ( maxTheta - minTheta ) / ( closed ? T : T - 1 );
	const double dP = ( maxPhi - minPhi ) / ( topCap && bottomCap ? P : P - 1 );

	for( int i = bottomCap ? 1 : 0; i != P; ++i )
	{
		const double phi = minPhi + dP * i;

		for( int j = 0; j != T; ++j )
		{
			const double theta = minTheta + dT * j;
			osg::Vec3d p( cos( phi ) * cos( theta ), sin( phi ), cos( phi ) * sin( theta ) );
			coords->push_back( p * r );
			p.normalize(); // ??
			normals->push_back( p );
		}
	}

	// add north pole point
	int offset = coords->size();;

	if( bottomCap )
	{
		coords->push_back( osg::Vec3d( 0, -r, 0 ) );
		normals->push_back( osg::Vec3d( 0, -1, 0 ) );
	}

	// add south pole point
	if( topCap )
	{
		coords->push_back( osg::Vec3d( 0, r, 0 ) );
		normals->push_back( osg::Vec3d( 0, 1, 0 ) );
	}

	// generate indices
	int maxPhiIndex = P - 2;
	if( topCap && !bottomCap ) ++maxPhiIndex;
	for( int i = 0; i != maxPhiIndex; ++i )
	{
		for( int j = 0; j != T; ++j )
		{
			const int id = j < T - 1 ? j + 1 : 0;
			if( j == T - 1 && !closed ) continue;
			indices->push_back( i * T + id );

			if( primitiveMode == osg::PrimitiveSet::QUADS )
			{
				indices->push_back( i *  T + j );
				indices->push_back( ( i + 1 ) *  T + j );
				indices->push_back( ( i + 1 ) *  T + id );
			}
			else if( primitiveMode == osg::PrimitiveSet::TRIANGLES )
			{
				if( i % 2 )
				{
					if( j % 2 )
					{
						indices->push_back( i * T + id );
						indices->push_back( i *  T + j );
						indices->push_back( ( i + 1 ) *  T + j );

						indices->push_back( ( i + 1 ) *  T + j );
						indices->push_back( ( i + 1 ) *  T + id );
						indices->push_back( i * T + id );
					}
					else
					{
						indices->push_back( i * T + id );
						indices->push_back( i *  T + j );
						indices->push_back( ( i + 1 ) *  T + id );

						indices->push_back( ( i + 1 ) *  T + j );
						indices->push_back( ( i + 1 ) *  T + id );
						indices->push_back( i * T + j);
					}
				}
				else
				{
					if( !( j % 2 ) )
					{
						indices->push_back( i * T + id );
						indices->push_back( i *  T + j );
						indices->push_back( ( i + 1 ) *  T + j );

						indices->push_back( ( i + 1 ) *  T + j );
						indices->push_back( ( i + 1 ) *  T + id );
						indices->push_back( i * T + id );
					}
					else
					{
						indices->push_back( i * T + id );
						indices->push_back( i *  T + j );
						indices->push_back( ( i + 1 ) *  T + id );

						indices->push_back( ( i + 1 ) *  T + j );
						indices->push_back( ( i + 1 ) *  T + id );
						indices->push_back( i * T + j);
					}
				}
			}
		}
	}

	// generate indices for bottom cap
	if( bottomCap )
	{
		osg::DrawElementsUInt* indicesBottom = new osg::DrawElementsUInt( osg::PrimitiveSet::TRIANGLE_FAN, 0);
		geom->addPrimitiveSet( indicesBottom );

		indicesBottom->push_back( offset );
		for( int i = 0; i != T; ++i ) indicesBottom->push_back( i );
		if( closed ) indicesBottom->push_back( 0 );
	}

	// generate indices for top cap
	if( topCap )
	{
		osg::DrawElementsUInt* indicesTop = new osg::DrawElementsUInt( osg::PrimitiveSet::TRIANGLE_FAN, 0);
		geom->addPrimitiveSet( indicesTop );

		indicesTop->push_back( bottomCap ? offset + 1 : offset );
		for( int i = 0; i != T; ++i ) indicesTop->push_back( offset - i - ( bottomCap ? 1 : 0 ) );
		if( closed ) indicesTop->push_back( offset - 1 );
	}

	return geode;
}


//----------------------------------------------------------------------------------
osg::Geometry* CreateDiscGeometry( double radius, int slices )
{

#ifndef M_PI //not in the C/C++ standards but defined in MS VC up to version 2003
  #define M_PI 3.1415926535897932384626433832795
#endif

	osg::ref_ptr< osg::Geometry > geom = new osg::Geometry;

	osg::Vec3Array* coords = new osg::Vec3Array;
	osg::Vec3Array* normals = new osg::Vec3Array;

	geom->setVertexArray( coords );
	geom->setNormalArray( normals );

	geom->setNormalBinding( osg::Geometry::BIND_PER_VERTEX );

	osg::DrawElementsUInt* indices = new osg::DrawElementsUInt( osg::PrimitiveSet::TRIANGLE_FAN, 0);
	geom->addPrimitiveSet( indices );

	/////
	coords->push_back( osg::Vec3d( 0, 0, 0 ) );
	normals->push_back( osg::Vec3d( 0, -1, 0 ) );
	indices->push_back( 0 );
	const int T = slices;
	const double r = radius;
	const double dT = 2 * M_PI / T;
	for( int j = 0; j != T; ++j )
	{
		const double theta = dT * j;
		osg::Vec3d p( cos( theta ), 0.0, sin( theta ) );
		coords->push_back( p * r );
		normals->push_back( p );
		indices->push_back( j + 1 );
	}
	indices->push_back( 1 );
	/////
	return geom.release();
}



//----------------------------------------------------------------------------------
/// Creates a cylinder as a quad vertex array centered at 'center'.
/// @param radius cylinder radius
/// @param open open/closed switch
/// @param slices number of slices equal to the number of rotation steps around
///        the Y axis
/// @param stacks number of stacks along the Y axis
/// @center cylinder center; the geometry span-range along the Y axis is
/// [center - height/2, center + height/2]
/// @return Geode containing a cylinder aligned along the Y axis
osg::Geode* CreateCylinderNode( double radius, double height, bool open, int slices, int stacks,
				const osg::Vec3d& center = osg::Vec3d( 0, 0, 0 ) )
{

#ifndef M_PI //not in the C/C++ standards but defined in MS VC up to version 2003
  #define M_PI 3.1415926535897932384626433832795
#endif

	osg::Geode* geode = new osg::Geode;
	osg::ref_ptr< osg::Geometry > geom = new osg::Geometry;

	geode->addDrawable( geom.get() );

	osg::Vec3Array* coords = new osg::Vec3Array;
	osg::Vec3Array* normals = new osg::Vec3Array;

	geom->setVertexArray( coords );
	geom->setNormalArray( normals );

	geom->setNormalBinding( osg::Geometry::BIND_PER_VERTEX );

	osg::DrawElementsUInt* indices = new osg::DrawElementsUInt( osg::PrimitiveSet::QUADS, 0);
	geom->addPrimitiveSet( indices );


	/////
	const int T = slices;
	const int P = stacks;
	const double r = radius;
	const double h = height;

	const double dT = 2 * M_PI / T;
	const double dh = h / ( P - 1 );

	for( int i = 0; i != P; ++i )
	{
		const double y = -0.5 * h + dh * i;

		for( int j = 0; j != T; ++j )
		{
			const double theta = dT * j;
			osg::Vec3d p( r * cos( theta ), y, r * sin( theta ) );
			coords->push_back( p + center );
			osg::Vec3d n( cos( theta ), 0, sin( theta ) );
			normals->push_back( n );
		}
	}


	for( int i = 0; i != P - 1; ++i )
	{
		for( int j = 0; j != T; ++j )
		{
			const int id = j < T - 1 ? j + 1 : 0;
			indices->push_back( i * T + id );
			indices->push_back( i *  T + j );
			indices->push_back( ( i + 1 ) *  T + j );
			indices->push_back( ( i + 1 ) *  T + id );
		}
	}

	if( !open )
	{
		osg::DrawElementsUInt* indicesTop = new osg::DrawElementsUInt( osg::PrimitiveSet::TRIANGLE_FAN, 0);
		geom->addPrimitiveSet( indicesTop );

		osg::DrawElementsUInt* indicesBottom = new osg::DrawElementsUInt( osg::PrimitiveSet::TRIANGLE_FAN, 0);
		geom->addPrimitiveSet( indicesBottom );

		int offset = coords->size();

		coords->push_back( osg::Vec3d( 0, -0.5 * h, 0 ) + center );
		normals->push_back( osg::Vec3d( 0, -1, 0 ) );

		coords->push_back( osg::Vec3d( 0, 0.5 * h, 0 )  + center );
		normals->push_back( osg::Vec3d( 0, 1, 0 ) );

		indicesBottom->push_back( offset );
		for( int i = 0; i != T; ++i ) indicesBottom->push_back( i );
		indicesBottom->push_back( 0 );

		indicesTop->push_back( offset + 1 );
		for( int i = 0; i != T; ++i ) indicesTop->push_back( offset - i - 1 );
		indicesTop->push_back( offset - 1 );
	}

	/////
	return geode;
}


//----------------------------------------------------------------------------------
/// Generates a transformation matrix that rotates, translates and scales an object
/// centered at the origin and enclosed inside the unit box aligned along a specific
/// axis. The scaling-translation-rotation and rotation-only matrices are returned
/// @param UP up vector (direction of alignment; i.e. (0, 1, 0) for an y-aligned cylinder
/// @param P0 first point
/// @param P1 second point
/// @param R returned rotation matrix that aligns the UP vector to (P1-P0)
/// @return transformation matrix that will position an object at (P0+P1)/2, rotated
/// such that the passed UP axis is directed along the P1-P0 direction and scaled
/// by length( P1 - P0 ) along the UP direction.
inline osg::Matrixd ComputeSTRMatrix( const osg::Vec3d& UP,
									  const osg::Vec3d& P0,
									  const osg::Vec3d& P1,
									  osg::Matrixd& R)
{
    const double d = ( P1 - P0 ).length();
    assert( d > 0.0 );
    osg::Matrixd S;
    S.makeScale( 1,  d, 1 );
    R.makeRotate( UP, ( P1 - P0 ) * 1.0 / d );
    osg::Matrixd M = S * R;
    M.setTrans( ( P1 + P0 ) * 0.5 );
	return M;
}



//------------------------------------------------------------------------------
osg::Geometry* CreateExtrudedSurfaceGeometry( osg::Vec3Array* section, osg::Vec3Array* path )
{

	// NULL ?
	assert( section );
	assert( path );
	assert( path->size() > 2 );

	// compute tangents
	osg::ref_ptr< osg::Vec3Array > tangents = new osg::Vec3Array;
	osg::Vec3d t = ( *path )[ 1 ] - ( *path )[ 0 ];
	t.normalize();
	tangents->push_back( t );
	for( int i = 1; i != path->size() - 1; ++i )
	{
		osg::Vec3d t = ( *path )[ i + 1 ] - ( *path )[ i - 1 ];
		t.normalize();
		tangents->push_back( t );
	}
	t = ( *path )[ path->size() - 1 ] - ( *path )[ path->size() - 2 ];
	t.normalize();
	tangents->push_back( t );

	// compute directions
	osg::ref_ptr< osg::Vec3Array > directions = new osg::Vec3Array;
	for( int i = 0; i != path->size() - 1; ++i )
	{
		t = ( *path )[ i + 1 ] - ( *path )[ i ];
		t.normalize();
		directions->push_back( t );
	}
	directions->push_back( tangents->back() );
	std::cout << '!' << std::endl;

	// compute rotated section
	osg::Vec3d normal = ( *section )[ 1 ] ^ ( *section )[ 0 ];
	if( normal * directions->front() < 0.0 ) normal = -normal;
	normal.normalize();
	osg::Matrixd r;
	r.makeRotate( normal, osg::Vec3d( directions->front() ) );

	osg::Matrix m( r );
	m.setTrans( path->front() );
	osg::ref_ptr< osg::Vec3Array > rsection = new osg::Vec3Array;
	osg::ref_ptr< osg::Vec3Array > cross = new osg::Vec3Array;
	for( osg::Vec3Array::iterator i = section->begin(); i != section->end(); ++i )
	{
		const osg::Vec3d t = *i * m;
		const osg::Vec3d c = *i * r;
		rsection->push_back( t );
		cross->push_back( c );
	}
	normal = normal * r;
	// for each path point:
	// 1) compute plane through point oriented in tangent direction
	// 2) compute new section by intersecting line from previous
	//    section points with plane
	// 3) add points to geometry
	std::vector< osg::ref_ptr< osg::Vec3Array > > sections;
	sections.push_back( rsection );

	for( int i = 1; i != path->size(); ++i )
	{
		osg::ref_ptr< osg::Vec3Array > newSection = new osg::Vec3Array;
		osg::Matrixd M;
		M.makeRotate( ( *tangents )[ i - 1 ], ( *tangents )[ i ] );
		osg::Matrixd Ti;
		osg::Matrixd T;
		Ti.setTrans( -( *path )[ i - 1 ] );
		T.setTrans( ( *path )[ i ] );
		osg::ref_ptr< osg::Vec3Array > prevSection = sections.back();
		for( osg::Vec3Array::iterator si = prevSection->begin(); si != prevSection->end(); ++si )
		{
			const osg::Vec3d startPoint = *si;
			newSection->push_back( ( *si ) * Ti * M * T );
		}
		sections.push_back( newSection );
	}
	std::cout << '!' << std::endl;
	osg::ref_ptr< osg::Vec3Array > vertices = new osg::Vec3Array;
	std::back_insert_iterator< osg::Vec3Array > ii( *vertices );
	std::cout << '!' << std::endl;
	std::vector< osg::ref_ptr< osg::Vec3Array > >::iterator si = sections.begin();
	std::cout << '!' << std::endl;
	for( ; si != sections.end(); ++si ) std::copy( ( *si )->begin(), ( *si )->end(), ii );
	osg::ref_ptr< osg::Geometry > geom = new osg::Geometry;
	geom->setVertexArray( vertices.get() );
	std::cout << '!' << std::endl;
	// create indices
	const int stripes = section->size() - 1;
	for( int s = 0; s != stripes ; ++s )
	{
		osg::DrawElementsUInt* indices = new osg::DrawElementsUInt( osg::PrimitiveSet::QUADS, 0);
		geom->addPrimitiveSet( indices );
		for( int p = 0; p != path->size() - 1; ++p )
		{
			indices->push_back( p * section->size() + s + 1 );
			indices->push_back( p * section->size() + s );
			indices->push_back( ( p + 1 ) * section->size() + s );
			indices->push_back( ( p + 1 ) * section->size() + s + 1 );
		}
	}
	std::cout << '!' << std::endl;
	return geom.release();
}

std::vector< osg::ref_ptr< osg::Geode > > gAtoms;
std::vector< osg::Vec4d > gAtomColors;

//----------------------------------------------------------------------------------
/// Create an osg::Group containing two geodes representing two halves of a chemical bond and
/// optionally two half spherical caps for liquorice type representations.
/// @param P0 center of first atom
/// @param atomicNum0 atomic number of first atom; used to select the proper state set
/// @param P1 center of second atom
/// @param atomicNum1 atomic number of second atom; used to select the proper state set
/// @param atom0 close bottom end of cylinder with spherical cap
/// @param atom1 close top end of cylinder with spherical cap
/// @return group with two cylinders representing two halves of a bond and optionally
/// two half spheres that act as caps for the generated cylinder.
osg::Group* CreateBondGroup( const osg::Vec3d& P0,
							 int atomicNum0,
							 const osg::Vec3d& P1,
							 int atomicNum1,
							 bool atom0,
							 bool atom1 )
{
    const osg::Vec3d half = ( P0 + P1   ) * 0.5;
    const osg::Vec3d c0   = ( P0 + half ) * 0.5;
    const osg::Vec3d c1   = ( half + P1 ) * 0.5;

    static const double radius = 0.2;
    static osg::ref_ptr< osg::Geode > cylinder = CreateCylinderNode( radius, 1, true, 16, 2, osg::Vec3d( 0, 0, 0 ) );
    static osg::ref_ptr< osg::Geode > bottomSphere = CreateSphereNode( radius, 16, 16, osg::PrimitiveSet::QUADS, -180, 180, -90, 0 );
	static osg::ref_ptr< osg::Geode > topSphere = CreateSphereNode( radius, 16, 16, osg::PrimitiveSet::QUADS, -180, 180, 0, 90 );

    osg::ref_ptr< osg::Group > group = new osg::Group;

	osg::Matrixd R;
    // bond 0
	osg::Matrixd t0 = ComputeSTRMatrix( osg::Vec3d( 0, 1, 0 ), P0, half, R );
    osg::ref_ptr< osg::Geode > b0 = new osg::Geode;
    b0->setStateSet( gAtoms[ atomicNum0 ]->getStateSet() );
    osg::Geometry* geom0 = cylinder->getDrawable( 0 )->asGeometry();
    osg::Geometry* g0 = new osg::Geometry;
    b0->addDrawable( geom0 );
    osg::MatrixTransform* bt0 = new osg::MatrixTransform;
    bt0->setMatrix( t0 );
    bt0->addChild( b0.get() );
    group->addChild( bt0 );

    // bond 1
    osg::Matrixd t1 = ComputeSTRMatrix( osg::Vec3d( 0, 1, 0 ), half, P1, R );
    osg::ref_ptr< osg::Geode > b1 = new osg::Geode;
    b1->setStateSet( gAtoms[ atomicNum1 ]->getStateSet() );
    osg::Geometry* geom1 = cylinder->getDrawable( 0 )->asGeometry();
    osg::Geometry* g1 = new osg::Geometry;
    b1->addDrawable( geom1 );
    osg::MatrixTransform* bt1 = new osg::MatrixTransform;
    bt1->setMatrix( t1 );
    bt1->addChild( b1.get() );
    group->addChild( bt1 );

    if( atom0 )
    {
	    // sphere 0
	    osg::Geode* a0 = new osg::Geode;
		a0->setStateSet( gAtoms[ atomicNum0 ]->getStateSet() );
	    osg::Geometry* geom0 = bottomSphere->getDrawable( 0 )->asGeometry();
	    osg::Geometry* a0geom = new osg::Geometry;
	    a0->addDrawable( geom0 );
		osg::MatrixTransform* at0 = new osg::MatrixTransform;
		osg::Matrixd AM0 = R;
		AM0.setTrans( P0 );
		at0->setMatrix( AM0 );
	    at0->addChild( a0 );
	    group->addChild( at0 );
	}

	if( atom1 )
	{
		// sphere 1
	    osg::Geode* a1 = new osg::Geode;
	    a1->setStateSet( gAtoms[ atomicNum1 ]->getStateSet() );
		osg::Geometry* geom1 = topSphere->getDrawable( 0 )->asGeometry();
	    osg::Geometry* a1geom = new osg::Geometry;
	    a1->addDrawable( geom1 );
	   	osg::MatrixTransform* at1 = new osg::MatrixTransform;
		osg::Matrixd AM1 = R;
		AM1.setTrans( P1 );
		at1->setMatrix( AM1 );
	    at1->addChild( a1 );
	    group->addChild( at1 );
	}

    return group.release();

}

//----------------------------------------------------------------------------------
osg::Group* CreateBondGroup2( const osg::Vec3d& P0,
							  int atomicNum0,
							  const osg::Vec3d& P1,
							  int atomicNum1,
							  bool atom0,
							  bool atom1 )
{
    const osg::Vec3d half = ( P0 + P1   ) * 0.5;
    const osg::Vec3d c0   = ( P0 + half ) * 0.5;
    const osg::Vec3d c1   = ( half + P1 ) * 0.5;

	osg::Matrixd R;
    // bond 0
	osg::Matrixd t0 = ComputeSTRMatrix( osg::Vec3d( 0, 1, 0 ), P0, half, R );

    static const double radius = 0.2;
	const double l = ( P1 - P0 ).length();
	osg::ref_ptr< osg::Vec3Array > vertices  = new osg::Vec3Array;
	osg::ref_ptr< osg::Vec2Array > texcoords = new osg::Vec2Array;
	vertices->push_back( osg::Vec3d( -.2, -.5 * l, 0 ) * R );
	vertices->push_back( osg::Vec3d(  .2, -.5 * l, 0 ) * R );
	vertices->push_back( osg::Vec3d(  .2,  .5 * l, 0 ) * R );
	vertices->push_back( osg::Vec3d( -.2,  .5 * l, 0 ) * R );
	osg::Vec3d normal = osg::Vec3d( 0, 0, 1 ) * R;

	texcoords->push_back( osg::Vec2d( 0, 0 ) );
	texcoords->push_back( osg::Vec2d( 1, 0 ) );
	texcoords->push_back( osg::Vec2d( 1, 1 ) );
	texcoords->push_back( osg::Vec2d( 0, 1 ) );

	osg::ref_ptr< osg::Geometry > geom = new osg::Geometry;
	geom->setVertexArray( vertices.get() );
	const int textureUnit = 0;
	geom->setTexCoordArray( textureUnit, texcoords.get() );
	geom->addPrimitiveSet( new osg::DrawArrays( osg::PrimitiveSet::QUADS, 0, 4 ) );

	osg::ref_ptr< osg::Billboard > b = new osg::Billboard;
	b->setMode(osg::Billboard::AXIAL_ROT );
	osg::Vec3d V = P0 - P1;
	V.normalize();
    b->setAxis( V );
    b->setNormal( normal );
    b->addDrawable( geom.get(), ( P0 + P1 ) * 0.5 );
    osg::ref_ptr< osg::Group > group = new osg::Group;

	osg::StateSet* set = b->getOrCreateStateSet();
	static osg::Texture2D* tex = 0;
	if( !tex )
	{
		tex = new osg::Texture2D();
		tex->setImage(osgDB::readImageFile( "./cylinder.rgb" ) );
	}
	//osg::TexEnv* blendTexEnv = new osg::TexEnv;
	//blendTexEnv->setMode(osg::TexEnv::BLEND);
    //set->setTextureAttribute( 0, blendTexEnv, osg::StateAttribute::ON );
	set->setTextureAttributeAndModes( 0, tex, osg::StateAttribute::ON );
	//set->setAttributeAndModes( new osg::BlendFunc, osg::StateAttribute::ON );
	set->setMode(GL_LIGHTING, osg::StateAttribute::OFF);

	group->addChild( b.get() );

	if( atom0 )
	{
		osg::Billboard* a1 = new osg::Billboard;
		a1->setMode( osg::Billboard::POINT_ROT_EYE );
		osg::Geometry* geom = new osg::Geometry;
		a1->addDrawable( geom, P0 );
		geom->setStateSet( gAtoms[ atomicNum0 ]->getStateSet() );
		osg::Geometry* g = gAtoms[ atomicNum0 ]->getDrawable( 0 )->asGeometry();
		geom->setColorArray( g->getColorArray() );
		geom->setColorBinding( osg::Geometry::BIND_OVERALL );
		osg::Vec3Array* coord = new osg::Vec3Array;
		coord->push_back( osg::Vec3d( -radius, 0, -radius ) );
		coord->push_back( osg::Vec3d( radius, 0, -radius ) );
		coord->push_back( osg::Vec3d( radius, 0, radius ) );
		coord->push_back( osg::Vec3d( -radius, 0, radius ) );
		geom->setVertexArray( coord );
		geom->addPrimitiveSet( new osg::DrawArrays( osg::PrimitiveSet::QUADS, 0, 4 ) );
		osg::Vec2Array* tex = new osg::Vec2Array;
		tex->push_back( osg::Vec2d( 0, 0 ) );
		tex->push_back( osg::Vec2d( 1, 0 ) );
		tex->push_back( osg::Vec2d( 1, 1 ) );
		tex->push_back( osg::Vec2d( 0, 1 ) );
		geom->setTexCoordArray( 0, tex );
		group->addChild( a1 );
	}

	if( atom1 )
	{
		osg::Billboard* a1 = new osg::Billboard;
		a1->setMode( osg::Billboard::POINT_ROT_EYE );
		osg::Geometry* geom = new osg::Geometry;
		a1->addDrawable( geom, P1 );
		geom->setStateSet( gAtoms[ atomicNum1 ]->getStateSet() );
		osg::Geometry* g = gAtoms[ atomicNum1 ]->getDrawable( 0 )->asGeometry();
		geom->setColorArray( g->getColorArray() );
		geom->setColorBinding( osg::Geometry::BIND_OVERALL );
		osg::Vec3Array* coord = new osg::Vec3Array;
		coord->push_back( osg::Vec3d( -radius, 0, -radius ) );
		coord->push_back( osg::Vec3d( radius, 0, -radius ) );
		coord->push_back( osg::Vec3d( radius, 0, radius ) );
		coord->push_back( osg::Vec3d( -radius, 0, radius ) );
		geom->setVertexArray( coord );
		geom->addPrimitiveSet( new osg::DrawArrays( osg::PrimitiveSet::QUADS, 0, 4 ) );
		osg::Vec2Array* tex = new osg::Vec2Array;
		tex->push_back( osg::Vec2d( 0, 0 ) );
		tex->push_back( osg::Vec2d( 1, 0 ) );
		tex->push_back( osg::Vec2d( 1, 1 ) );
		tex->push_back( osg::Vec2d( 0, 1 ) );
		geom->setTexCoordArray( 0, tex );
		group->addChild( a1 );
	}

	return group.release();
}
//----------------------------------------------------------------------------------

osg::Group* CreateBondGroup3( const osg::Vec3d& P0,
							  int atomicNum0,
							  const osg::Vec3d& P1,
							  int atomicNum1,
							  bool atom0,
							  bool atom1 )
{
    const osg::Vec3d half = ( P0 + P1   ) * 0.5;
    const osg::Vec3d c0   = ( P0 + half ) * 0.5;
    const osg::Vec3d c1   = ( half + P1 ) * 0.5;

	osg::Matrixd R;
    // bond 0
	osg::Matrixd t0 = ComputeSTRMatrix( osg::Vec3d( 0, 1, 0 ), P0, half, R );

    static const double radius = 0.2;
	const double l = ( P1 - P0 ).length();
	osg::ref_ptr< osg::Vec3Array > vertices  = new osg::Vec3Array;
	osg::ref_ptr< osg::Vec2Array > texcoords = new osg::Vec2Array;
	vertices->push_back( osg::Vec3d( -radius, -.5 * l, 0 ) * R );
	vertices->push_back( osg::Vec3d(  radius, -.5 * l, 0 ) * R );
	vertices->push_back( osg::Vec3d(  radius,  .5 * l, 0 ) * R );
	vertices->push_back( osg::Vec3d( -radius,  .5 * l, 0 ) * R );
	osg::Vec3d normal = osg::Vec3d( 0, 0, 1 ) * R;

	texcoords->push_back( osg::Vec2d( 0, 0 ) );
	texcoords->push_back( osg::Vec2d( 1, 0 ) );
	texcoords->push_back( osg::Vec2d( 1, 1 ) );
	texcoords->push_back( osg::Vec2d( 0, 1 ) );

	osg::ref_ptr< osg::Geometry > geom = new osg::Geometry;
	geom->setVertexArray( vertices.get() );
	const int textureUnit = 0;
	geom->setTexCoordArray( textureUnit, texcoords.get() );
	geom->addPrimitiveSet( new osg::DrawArrays( osg::PrimitiveSet::QUADS, 0, 4 ) );

	osg::ref_ptr< osg::Billboard > b = new osg::Billboard;
	b->setMode(osg::Billboard::AXIAL_ROT );
	osg::Vec3d V = P0 - P1;
	V.normalize();
    b->setAxis( V );
    b->setNormal( normal );
    b->addDrawable( geom.get(), ( P0 + P1 ) * 0.5 );
    osg::ref_ptr< osg::Group > group = new osg::Group;

	osg::StateSet* set = b->getOrCreateStateSet();
	osg::Program* program = new osg::Program;
    program->setName( "cylinder" );
	program->addShader( osg::Shader::readShaderFile( osg::Shader::FRAGMENT, "cylinder2.frag" ) );
	program->addShader( osg::Shader::readShaderFile( osg::Shader::VERTEX, "billboard2.vert" ) );
    set->setAttributeAndModes( program, osg::StateAttribute::ON );
	set->addUniform( new osg::Uniform( "radius", float( radius ) ) );
	set->addUniform( new osg::Uniform( "length", float( l ) ) );
	set->addUniform( new osg::Uniform( "invObjectMatrix", osg::Matrixd::inverse( R ) ) );
	set->addUniform( new osg::Uniform( "objectMatrix", R ) );

	//osg::Texture2D *tex = new osg::Texture2D();
 //   tex->setImage(osgDB::readImageFile( "./atom.png" ) );
 //   set->setTextureAttributeAndModes( 0, tex, osg::StateAttribute::ON );
	//set->addUniform( new osg::Uniform( "texture", tex ) );

	//set->setMode(GL_BLEND, osg::StateAttribute::ON);
	set->setMode(GL_LIGHTING, osg::StateAttribute::OFF);

	group->addChild( b.get() );
	return group.release();
}


//----------------------------------------------------------------------------------
osg::Group* CreateBondGroup4( const osg::Vec3d& P0,
							  int atomicNum0,
							  const osg::Vec3d& P1,
							  int atomicNum1,
							  bool atom0,
							  bool atom1 )
{
    const osg::Vec3d half = ( P0 + P1   ) * 0.5;
    const osg::Vec3d c0   = ( P0 + half ) * 0.5;
    const osg::Vec3d c1   = ( half + P1 ) * 0.5;

	osg::Matrixd R;
    // bond 0
	osg::Matrixd t0 = ComputeSTRMatrix( osg::Vec3d( 0, 1, 0 ), P0, half, R );

    static const double radius = 0.2;
	const double l = ( P1 - P0 ).length();

	// compute billboard
	const double r = radius;
	osg::BoundingBox aabox;
	aabox.expandBy( osg::Vec3d( -r, -.5 * l, -r ) * R );
	aabox.expandBy( osg::Vec3d(  r, -.5 * l, -r ) * R );
	aabox.expandBy( osg::Vec3d(  r,  .5 * l, -r ) * R );
	aabox.expandBy( osg::Vec3d( -r,  .5 * l, -r ) * R );
	aabox.expandBy( osg::Vec3d( -r, -.5 * l,  r ) * R );
	aabox.expandBy( osg::Vec3d(  r, -.5 * l,  r ) * R );
	aabox.expandBy( osg::Vec3d(  r,  .5 * l,  r ) * R );
	aabox.expandBy( osg::Vec3d( -r,  .5 * l,  r ) * R );
	const double hex = .5 *( aabox.xMax() - aabox.xMin() );
	const double hey = .5 *( aabox.yMax() - aabox.yMin() );
	const double h = std::max( hex, hey );

	osg::ref_ptr< osg::Vec3Array > vertices = new osg::Vec3Array;
	vertices->push_back( osg::Vec3d( -h, 0, -h ) );
	vertices->push_back( osg::Vec3d(  h, 0, -h ) );
	vertices->push_back( osg::Vec3d(  h, 0,  h ) );
	vertices->push_back( osg::Vec3d( -h, 0,  h ) );

	osg::ref_ptr< osg::Vec2Array > texcoords = new osg::Vec2Array;
	texcoords->push_back( osg::Vec2d( 0, 0 ) );
	texcoords->push_back( osg::Vec2d( 1, 0 ) );
	texcoords->push_back( osg::Vec2d( 1, 1 ) );
	texcoords->push_back( osg::Vec2d( 0, 1 ) );

	osg::ref_ptr< osg::Geometry > geom = new osg::Geometry;
	geom->setVertexArray( vertices.get() );
	const int textureUnit = 0;
	geom->setTexCoordArray( textureUnit, texcoords.get() );
	geom->addPrimitiveSet( new osg::DrawArrays( osg::PrimitiveSet::QUADS, 0, 4 ) );

	osg::ref_ptr< osg::Billboard > b = new osg::Billboard;
	b->setMode( osg::Billboard::POINT_ROT_EYE );
	b->addDrawable( geom.get(), half );
	osg::StateSet* set = b->getOrCreateStateSet();
	set->setMode(GL_LIGHTING, osg::StateAttribute::OFF);
	osg::Texture2D *tex = new osg::Texture2D();
    tex->setImage(osgDB::readImageFile( "./cylinder.rgb" ) );
    set->setTextureAttributeAndModes( 0, tex, osg::StateAttribute::ON );
	osg::ref_ptr< osg::Vec4Array > colors = new osg::Vec4Array;
	colors->push_back( osg::Vec4d( 1, 1, 1, 1 ) );
	geom->setColorArray( colors.get() );
	geom->setColorBinding( osg::Geometry::BIND_OVERALL );

	osg::ref_ptr< osg::Group > group = new osg::Group;
	group->addChild( b.get() );

	// shaders
	osg::ref_ptr< osg::Program > program = new osg::Program;
    program->setName( "cylinder" );
	program->addShader( osg::Shader::readShaderFile( osg::Shader::FRAGMENT, "sphere.frag" ) );
	program->addShader( osg::Shader::readShaderFile( osg::Shader::VERTEX, "billboard.vert" ) );
    set->setAttributeAndModes( program.get(), osg::StateAttribute::ON );
	set->addUniform( new osg::Uniform( "radius", float( radius ) ) );
	set->addUniform( new osg::Uniform( "length", float( l ) ) );
	set->addUniform( new osg::Uniform( "objectMatrix", /*osg::Matrixd::inverse*/( R ) ) );
	set->addUniform( new osg::Uniform( "invObjectMatrix", osg::Matrixd::inverse( R ) ) );
	set->addUniform( new osg::Uniform( "center", half ) );

	return group.release();



	if( atom0 )
	{
		osg::Billboard* a1 = new osg::Billboard;
		a1->setMode( osg::Billboard::POINT_ROT_EYE );
		osg::Geometry* geom = new osg::Geometry;
		a1->addDrawable( geom, P0 );
		geom->setStateSet( gAtoms[ atomicNum0 ]->getStateSet() );
		osg::Geometry* g = gAtoms[ atomicNum0 ]->getDrawable( 0 )->asGeometry();
		geom->setColorArray( g->getColorArray() );
		geom->setColorBinding( osg::Geometry::BIND_OVERALL );
		osg::Vec3Array* coord = new osg::Vec3Array;
		coord->push_back( osg::Vec3d( -radius, 0, -radius ) );
		coord->push_back( osg::Vec3d( radius, 0, -radius ) );
		coord->push_back( osg::Vec3d( radius, 0, radius ) );
		coord->push_back( osg::Vec3d( -radius, 0, radius ) );
		geom->setVertexArray( coord );
		geom->addPrimitiveSet( new osg::DrawArrays( osg::PrimitiveSet::QUADS, 0, 4 ) );
		osg::Vec2Array* tex = new osg::Vec2Array;
		tex->push_back( osg::Vec2d( 0, 0 ) );
		tex->push_back( osg::Vec2d( 1, 0 ) );
		tex->push_back( osg::Vec2d( 1, 1 ) );
		tex->push_back( osg::Vec2d( 0, 1 ) );
		geom->setTexCoordArray( 0, tex );
		group->addChild( a1 );
	}

	if( atom1 )
	{
		osg::Billboard* a1 = new osg::Billboard;
		a1->setMode( osg::Billboard::POINT_ROT_EYE );
		osg::Geometry* geom = new osg::Geometry;
		a1->addDrawable( geom, P1 );
		geom->setStateSet( gAtoms[ atomicNum1 ]->getStateSet() );
		osg::Geometry* g = gAtoms[ atomicNum1 ]->getDrawable( 0 )->asGeometry();
		geom->setColorArray( g->getColorArray() );
		geom->setColorBinding( osg::Geometry::BIND_OVERALL );
		osg::Vec3Array* coord = new osg::Vec3Array;
		coord->push_back( osg::Vec3d( -radius, 0, -radius ) );
		coord->push_back( osg::Vec3d( radius, 0, -radius ) );
		coord->push_back( osg::Vec3d( radius, 0, radius ) );
		coord->push_back( osg::Vec3d( -radius, 0, radius ) );
		geom->setVertexArray( coord );
		geom->addPrimitiveSet( new osg::DrawArrays( osg::PrimitiveSet::QUADS, 0, 4 ) );
		osg::Vec2Array* tex = new osg::Vec2Array;
		tex->push_back( osg::Vec2d( 0, 0 ) );
		tex->push_back( osg::Vec2d( 1, 0 ) );
		tex->push_back( osg::Vec2d( 1, 1 ) );
		tex->push_back( osg::Vec2d( 0, 1 ) );
		geom->setTexCoordArray( 0, tex );
		group->addChild( a1 );
	}

	return group.release();
}

//----------------------------------------------------------------------------------
osg::Group* CreateBondGroup5( const osg::Vec3d& P0,
							  int atomicNum0,
							  const osg::Vec3d& P1,
							  int atomicNum1,
							  bool atom0,
							  bool atom1 )
{
    const osg::Vec3d half = ( P0 + P1   ) * 0.5;
    const osg::Vec3d c0   = ( P0 + half ) * 0.5;
    const osg::Vec3d c1   = ( half + P1 ) * 0.5;

	osg::Matrixd R;
    // bond 0
	osg::Matrixd t0 = ComputeSTRMatrix( osg::Vec3d( 0, 1, 0 ), P0, half, R );

    static const double radius = 0.2;
	const double l = ( P1 - P0 ).length();

	// compute billboard
	const double r = radius;
	osg::BoundingBox aabox;
	aabox.expandBy( osg::Vec3d( -r, -.5 * l, -r ) * R );
	aabox.expandBy( osg::Vec3d(  r, -.5 * l, -r ) * R );
	aabox.expandBy( osg::Vec3d(  r,  .5 * l, -r ) * R );
	aabox.expandBy( osg::Vec3d( -r,  .5 * l, -r ) * R );
	aabox.expandBy( osg::Vec3d( -r, -.5 * l,  r ) * R );
	aabox.expandBy( osg::Vec3d(  r, -.5 * l,  r ) * R );
	aabox.expandBy( osg::Vec3d(  r,  .5 * l,  r ) * R );
	aabox.expandBy( osg::Vec3d( -r,  .5 * l,  r ) * R );
	const double hex = .5 *( aabox.xMax() - aabox.xMin() );
	const double hey = .5 *( aabox.yMax() - aabox.yMin() );
	const double hez = .5 *( aabox.zMax() - aabox.zMin() );
	const double h = std::max( hex, std::max( hey, hez ) );

	osg::ref_ptr< osg::Vec3Array > vertices = new osg::Vec3Array;
	vertices->push_back( osg::Vec3d( -h, 0, -h ) );
	vertices->push_back( osg::Vec3d(  h, 0, -h ) );
	vertices->push_back( osg::Vec3d(  h, 0,  h ) );
	vertices->push_back( osg::Vec3d( -h, 0,  h ) );

	osg::ref_ptr< osg::Vec2Array > texcoords = new osg::Vec2Array;
	texcoords->push_back( osg::Vec2d( 0, 0 ) );
	texcoords->push_back( osg::Vec2d( 1, 0 ) );
	texcoords->push_back( osg::Vec2d( 1, 1 ) );
	texcoords->push_back( osg::Vec2d( 0, 1 ) );

	osg::ref_ptr< osg::Geometry > geom = new osg::Geometry;
	geom->setVertexArray( vertices.get() );
	const int textureUnit = 0;
	geom->setTexCoordArray( textureUnit, texcoords.get() );
	geom->addPrimitiveSet( new osg::DrawArrays( osg::PrimitiveSet::QUADS, 0, 4 ) );

	osg::ref_ptr< osg::Billboard > b = new osg::Billboard;
	b->setMode( osg::Billboard::POINT_ROT_EYE );
	b->addDrawable( geom.get(), half );
	osg::StateSet* set = b->getOrCreateStateSet();
	set->setMode(GL_LIGHTING, osg::StateAttribute::OFF);
	osg::Texture2D *tex = new osg::Texture2D();
    tex->setImage(osgDB::readImageFile( "./cylinder.rgb" ) );
    set->setTextureAttributeAndModes( 0, tex, osg::StateAttribute::ON );
	osg::ref_ptr< osg::Vec4Array > colors = new osg::Vec4Array;
	colors->push_back( osg::Vec4d( 1, 1, 1, 1 ) );
	geom->setColorArray( colors.get() );
	geom->setColorBinding( osg::Geometry::BIND_OVERALL );

	osg::ref_ptr< osg::Group > group = new osg::Group;
	group->addChild( b.get() );

	// shaders
	osg::ref_ptr< osg::Program > program = new osg::Program;
    program->setName( "cylinder" );
	program->addShader( osg::Shader::readShaderFile( osg::Shader::FRAGMENT, "cylinder2.frag" ) );
	program->addShader( osg::Shader::readShaderFile( osg::Shader::VERTEX, "billboard2.vert" ) );
    set->setAttributeAndModes( program.get(), osg::StateAttribute::ON );
	set->addUniform( new osg::Uniform( "radius", float( radius ) ) );
	set->addUniform( new osg::Uniform( "length", float( l ) ) );
	set->addUniform( new osg::Uniform( "objectMatrix", /*osg::Matrixd::inverse*/( R ) ) );
	set->addUniform( new osg::Uniform( "invObjectMatrix", osg::Matrixd::inverse( R ) ) );
	set->addUniform( new osg::Uniform( "center", half ) );

	return group.release();



	if( atom0 )
	{
		osg::Billboard* a1 = new osg::Billboard;
		a1->setMode( osg::Billboard::POINT_ROT_EYE );
		osg::Geometry* geom = new osg::Geometry;
		a1->addDrawable( geom, P0 );
		geom->setStateSet( gAtoms[ atomicNum0 ]->getStateSet() );
		osg::Geometry* g = gAtoms[ atomicNum0 ]->getDrawable( 0 )->asGeometry();
		geom->setColorArray( g->getColorArray() );
		geom->setColorBinding( osg::Geometry::BIND_OVERALL );
		osg::Vec3Array* coord = new osg::Vec3Array;
		coord->push_back( osg::Vec3d( -radius, 0, -radius ) );
		coord->push_back( osg::Vec3d( radius, 0, -radius ) );
		coord->push_back( osg::Vec3d( radius, 0, radius ) );
		coord->push_back( osg::Vec3d( -radius, 0, radius ) );
		geom->setVertexArray( coord );
		geom->addPrimitiveSet( new osg::DrawArrays( osg::PrimitiveSet::QUADS, 0, 4 ) );
		osg::Vec2Array* tex = new osg::Vec2Array;
		tex->push_back( osg::Vec2d( 0, 0 ) );
		tex->push_back( osg::Vec2d( 1, 0 ) );
		tex->push_back( osg::Vec2d( 1, 1 ) );
		tex->push_back( osg::Vec2d( 0, 1 ) );
		geom->setTexCoordArray( 0, tex );
		group->addChild( a1 );
	}

	if( atom1 )
	{
		osg::Billboard* a1 = new osg::Billboard;
		a1->setMode( osg::Billboard::POINT_ROT_EYE );
		osg::Geometry* geom = new osg::Geometry;
		a1->addDrawable( geom, P1 );
		geom->setStateSet( gAtoms[ atomicNum1 ]->getStateSet() );
		osg::Geometry* g = gAtoms[ atomicNum1 ]->getDrawable( 0 )->asGeometry();
		geom->setColorArray( g->getColorArray() );
		geom->setColorBinding( osg::Geometry::BIND_OVERALL );
		osg::Vec3Array* coord = new osg::Vec3Array;
		coord->push_back( osg::Vec3d( -radius, 0, -radius ) );
		coord->push_back( osg::Vec3d( radius, 0, -radius ) );
		coord->push_back( osg::Vec3d( radius, 0, radius ) );
		coord->push_back( osg::Vec3d( -radius, 0, radius ) );
		geom->setVertexArray( coord );
		geom->addPrimitiveSet( new osg::DrawArrays( osg::PrimitiveSet::QUADS, 0, 4 ) );
		osg::Vec2Array* tex = new osg::Vec2Array;
		tex->push_back( osg::Vec2d( 0, 0 ) );
		tex->push_back( osg::Vec2d( 1, 0 ) );
		tex->push_back( osg::Vec2d( 1, 1 ) );
		tex->push_back( osg::Vec2d( 0, 1 ) );
		geom->setTexCoordArray( 0, tex );
		group->addChild( a1 );
	}

	return group.release();
}

//----------------------------------------------------------------------------------
#include "cylinder3.vert.h"
#include "cylinder3.frag.h"
#include "sphere3.frag.h"
#include "sphere3.vert.h"
// Add quads used as ray targets. Best results are obtained with two quads per primitive
// faster results with one quad put in front of the primitive.
osg::Group* CreateBondGroup6( const osg::Vec3d& P0,
							  int atomicNum0,
							  const osg::Vec3d& P1,
							  int atomicNum1,
							  bool atom0,
							  bool atom1 )
{
    const osg::Vec3d half = ( P0 + P1   ) * 0.5;
    const osg::Vec3d c0   = ( P0 + half ) * 0.5;
    const osg::Vec3d c1   = ( half + P1 ) * 0.5;

	osg::Matrixd R;
    // bond 0
	osg::Matrixd t0 = ComputeSTRMatrix( osg::Vec3d( 0, 1, 0 ), P0, half, R );

    static const double radius = 0.2;
	const double l = ( P1 - P0 ).length();

	osg::Vec3d axis = osg::Vec3d( 0, 1, 0 ) * R;
	axis.normalize();
	osg::Matrixd Q;
	Q.makeIdentity();
	// | a d e g |
	// | d b f h |
	// | e f c i |
	// | g h i j |
	Q( 1, 1 ) = 0;
	Q( 3, 3 ) = -radius * radius;
	/*Q( 1, 1 ) = 0;
	Q( 0, 0 ) = 4;
	Q( 2, 2 ) = 4;
	Q( 1, 3 ) = -.5;
    Q( 3, 1 ) = Q( 1, 3 );
	Q( 3, 3 ) = 0;*/

	osg::Matrixd QT = osg::Matrixd::inverse( R ) * Q * R ;

	// compute billboard
	const double r = radius;
	osg::BoundingBox aabox;
	aabox.expandBy( osg::Vec3d( -r, -.5 * l, -r ) * R );
	aabox.expandBy( osg::Vec3d(  r, -.5 * l, -r ) * R );
	aabox.expandBy( osg::Vec3d(  r,  .5 * l, -r ) * R );
	aabox.expandBy( osg::Vec3d( -r,  .5 * l, -r ) * R );
	aabox.expandBy( osg::Vec3d( -r, -.5 * l,  r ) * R );
	aabox.expandBy( osg::Vec3d(  r, -.5 * l,  r ) * R );
	aabox.expandBy( osg::Vec3d(  r,  .5 * l,  r ) * R );
	aabox.expandBy( osg::Vec3d( -r,  .5 * l,  r ) * R );
	const double hex = .5 *( aabox.xMax() - aabox.xMin() );
	const double hey = .5 *( aabox.yMax() - aabox.yMin() );
	const double hez = .5 *( aabox.zMax() - aabox.zMin() );
	const double h = std::max( hex, std::max( hey, hez ) );
	//const double h = aabox.radius();

	osg::ref_ptr< osg::Vec3Array > vertices = new osg::Vec3Array;

	// option 1 two back/front slabs
	// best accuracy
	//vertices->push_back( osg::Vec3d( -h, -h, -h ) );
	//vertices->push_back( osg::Vec3d(  h, -h, -h ) );
	//vertices->push_back( osg::Vec3d(  h, -h,  h ) );
	//vertices->push_back( osg::Vec3d( -h, -h,  h ) );
	//vertices->push_back( osg::Vec3d( -h, h, -h ) );
	//vertices->push_back( osg::Vec3d(  h, h, -h ) );
	//vertices->push_back( osg::Vec3d(  h, h,  h ) );
	//vertices->push_back( osg::Vec3d( -h, h,  h ) );


    // option 2 one quad
	// faster
	// put quad in front: billboards are automatically rotated such that -y is parallel to z in eye
	// coordinates
	vertices->push_back( osg::Vec3d( -h, -h, -h ) );
	vertices->push_back( osg::Vec3d(  h, -h, -h ) );
	vertices->push_back( osg::Vec3d(  h, -h,  h ) );
	vertices->push_back( osg::Vec3d( -h, -h,  h ) );

	//osg::ref_ptr< osg::Vec2Array > texcoords = new osg::Vec2Array;
	//texcoords->push_back( osg::Vec2d( 0, 0 ) );
	//texcoords->push_back( osg::Vec2d( 1, 0 ) );
	//texcoords->push_back( osg::Vec2d( 1, 1 ) );
	//texcoords->push_back( osg::Vec2d( 0, 1 ) );
	//texcoords->push_back( osg::Vec2d( 0, 0 ) );
	//texcoords->push_back( osg::Vec2d( 1, 0 ) );
	//texcoords->push_back( osg::Vec2d( 1, 1 ) );
	//texcoords->push_back( osg::Vec2d( 0, 1 ) );

	osg::ref_ptr< osg::Geometry > geom = new osg::Geometry;
	geom->setVertexArray( vertices.get() );
	const int textureUnit = 0;
	//geom->setTexCoordArray( textureUnit, texcoords.get() );
	geom->addPrimitiveSet( new osg::DrawArrays( osg::PrimitiveSet::QUADS, 0, 4 ) );

	osg::ref_ptr< osg::Billboard > b = new osg::Billboard;
	b->setMode( osg::Billboard::POINT_ROT_EYE );
	b->addDrawable( geom.get(), half );
	osg::StateSet* set = b->getOrCreateStateSet();
	set->setMode(GL_LIGHTING, osg::StateAttribute::OFF);
	//osg::Texture2D *tex = new osg::Texture2D();
    //tex->setImage(osgDB::readImageFile( "./cylinder.png" ) );
    //set->setTextureAttributeAndModes( 0, tex, osg::StateAttribute::ON );
	osg::ref_ptr< osg::Vec4Array > colors = new osg::Vec4Array;
	colors->push_back( osg::Vec4d( 1, 1, 1, 1 ) );
	geom->setColorArray( colors.get() );
	geom->setColorBinding( osg::Geometry::BIND_OVERALL );

	osg::ref_ptr< osg::Group > group = new osg::Group;
	group->addChild( b.get() );

	// shaders
	static osg::ref_ptr< osg::Program > program;
	if( program == 0 )
	{
		program = new osg::Program;
		program->setName( "cylinder" );
		//program->addShader( osg::Shader::readShaderFile( osg::Shader::FRAGMENT, "cylinder3.frag" ) );
		//program->addShader( osg::Shader::readShaderFile( osg::Shader::VERTEX, "billboard3.vert" ) );
		program->addShader( new osg::Shader( osg::Shader::FRAGMENT, BILLBOARD_CYLINDER3FRAG ) );
		program->addShader( new osg::Shader( osg::Shader::VERTEX,   BILLBOARD_CYLINDER3VERT ) );
	}
	set->setAttributeAndModes( program.get(), osg::StateAttribute::ON );
	set->addUniform( new osg::Uniform( "hlength", float( 0.5 * l ) ) );
	set->addUniform( new osg::Uniform( "Q", QT ) );
	set->addUniform( new osg::Uniform( "axis", axis ) );
	set->addUniform( new osg::Uniform( "bcolor", gAtomColors[ atomicNum0 ] ) );
	set->addUniform( new osg::Uniform( "tcolor", gAtomColors[ atomicNum1 ] ) );
	//return group.release();


	static osg::ref_ptr< osg::Program > aprogram;
	if( aprogram == 0 )
	{
		aprogram = new osg::Program;
		aprogram->setName( "sphere" );
		aprogram->addShader( new osg::Shader( osg::Shader::FRAGMENT, SPHEREFRAG ) );
		aprogram->addShader( new osg::Shader( osg::Shader::VERTEX, SPHEREVERT ) );
	}

	if( atom0 )
	{
		osg::Billboard* a1 = new osg::Billboard;
		osg::StateSet* set = a1->getOrCreateStateSet();
		set->setAttributeAndModes( aprogram.get(), osg::StateAttribute::ON );
	    set->addUniform( new osg::Uniform( "radius", float( radius ) ) );
		set->setMode(GL_LIGHTING, osg::StateAttribute::OFF);

		a1->setMode( osg::Billboard::POINT_ROT_EYE );
		osg::Geometry* geom = new osg::Geometry;
		a1->addDrawable( geom, P0 );

		osg::Geometry* g = gAtoms[ atomicNum0 ]->getDrawable( 0 )->asGeometry();
		geom->setColorArray( g->getColorArray() );
		geom->setColorBinding( osg::Geometry::BIND_OVERALL );
		osg::Vec3Array* coord = new osg::Vec3Array;
		// put quad in front: billboards are automatically rotated such that -y is parallel to z in eye
		// coordinates
		coord->push_back( osg::Vec3d( -radius, -radius, -radius ) );
		coord->push_back( osg::Vec3d( radius, -radius, -radius ) );
		coord->push_back( osg::Vec3d( radius, -radius, radius ) );
		coord->push_back( osg::Vec3d( -radius, -radius, radius ) );
		geom->setVertexArray( coord );
		geom->addPrimitiveSet( new osg::DrawArrays( osg::PrimitiveSet::QUADS, 0, 4 ) );
		//osg::Vec2Array* tex = new osg::Vec2Array;
		//tex->push_back( osg::Vec2d( 0, 0 ) );
		//tex->push_back( osg::Vec2d( 1, 0 ) );
		//tex->push_back( osg::Vec2d( 1, 1 ) );
		//tex->push_back( osg::Vec2d( 0, 1 ) );
		//geom->setTexCoordArray( 0, tex );
		group->addChild( a1 );
	}

	if( atom1 )
	{

		osg::Billboard* a1 = new osg::Billboard;
		osg::StateSet* set = a1->getOrCreateStateSet();
		set->setAttributeAndModes( aprogram.get(), osg::StateAttribute::ON );
	    set->addUniform( new osg::Uniform( "radius", float( radius ) ) );
		set->setMode(GL_LIGHTING, osg::StateAttribute::OFF);

		a1->setMode( osg::Billboard::POINT_ROT_EYE );
		osg::Geometry* geom = new osg::Geometry;
		a1->addDrawable( geom, P1 );

		osg::Geometry* g = gAtoms[ atomicNum1 ]->getDrawable( 0 )->asGeometry();
		geom->setColorArray( g->getColorArray() );
		geom->setColorBinding( osg::Geometry::BIND_OVERALL );
		osg::Vec3Array* coord = new osg::Vec3Array;
		// put quad in front: billboards are automatically rotated such that -y is parallel to z in eye
		// coordinates
		coord->push_back( osg::Vec3d( -radius, -radius, -radius ) );
		coord->push_back( osg::Vec3d(  radius, -radius, -radius ) );
		coord->push_back( osg::Vec3d(  radius, -radius,  radius ) );
		coord->push_back( osg::Vec3d( -radius, -radius,  radius ) );
		geom->setVertexArray( coord );
		geom->addPrimitiveSet( new osg::DrawArrays( osg::PrimitiveSet::QUADS, 0, 4 ) );
		//osg::Vec2Array* tex = new osg::Vec2Array;
		//tex->push_back( osg::Vec2d( 0, 0 ) );
		//tex->push_back( osg::Vec2d( 1, 0 ) );
		//tex->push_back( osg::Vec2d( 1, 1 ) );
		//tex->push_back( osg::Vec2d( 0, 1 ) );
		//geom->setTexCoordArray( 0, tex );
		group->addChild( a1 );
	}

	return group.release();
}



//------------------------------------------------------------------------------
//#include "cylinder4.vert.h"
//#include "cylinder4.frag.h"
// BOX
osg::Group* CreateBondGroup7( const osg::Vec3d& P0,
							  int atomicNum0,
							  const osg::Vec3d& P1,
							  int atomicNum1,
							  bool atom0,
							  bool atom1 )
{
    const osg::Vec3d half = ( P0 + P1   ) * 0.5;

	osg::Matrixd R;
    // bond 0
	osg::Matrixd t0 = ComputeSTRMatrix( osg::Vec3d( 0, 1, 0 ), P0, half, R );
    static const double radius = 0.2;
	const double r = radius;
	const double l = ( P1 - P0 ).length();

	osg::Vec3d axis = osg::Vec3d( 0, 1, 0 ) * R;
	axis.normalize();
	osg::Matrixd Q;
	Q.makeIdentity();
	// | a d e g |
	// | d b f h |
	// | e f c i |
	// | g h i j |
	Q( 1, 1 ) = 0;
	Q( 3, 3 ) = -radius * radius;
	osg::Matrixd QT = osg::Matrixd::inverse( R ) * Q * R ;

	osg::ref_ptr< osg::Group > group = new osg::Group;
	// OPTION 1 oriented box: pass center to shader
	//osg::Vec3 c0( -r, -.5*l, -r ), c1( r, -.5*l, -r ), c2( r, .5*l, -r ), c3( -r, .5*l, -r );
	//osg::Vec3 c4( -r, -.5*l, r ), c5( r, -.5*l, r ), c6( r, .5*l, r ), c7( -r, .5*l, r );
	//osg::ref_ptr< osg::Geometry > geom = new osg::Geometry;
	//osg::ref_ptr< osg::Vec3Array > coord = new osg::Vec3Array;
	////back
	//coord->push_back( c0 * R + half ); coord->push_back( c3 * R + half ); coord->push_back( c2 * R + half ); coord->push_back( c1 * R + half );
	////left
	//coord->push_back( c0 * R + half ); coord->push_back( c4 * R + half ); coord->push_back( c7 * R + half ); coord->push_back( c3 * R + half );
	////right
	//coord->push_back( c5 * R + half ); coord->push_back( c1 * R + half ); coord->push_back( c2 * R + half ); coord->push_back( c6 * R + half );
	//// front
	//coord->push_back( c4 * R + half ); coord->push_back( c5 * R + half ); coord->push_back( c6 * R + half ); coord->push_back( c7 * R + half );
	////top
	//coord->push_back( c7 * R + half ); coord->push_back( c6 * R + half ); coord->push_back( c2 * R + half ); coord->push_back( c3 * R + half );
	////bottom
	//coord->push_back( c4 * R + half ); coord->push_back( c0 * R + half ); coord->push_back( c1 * R + half ); coord->push_back( c5 * R + half );
	//geom->setVertexArray( coord.get() );
	//geom->addPrimitiveSet( new osg::DrawArrays( osg::PrimitiveSet::QUADS, 0, 24 ) );
	//osg::ref_ptr< osg::Geode > b = new osg::Geode;
	//b->addDrawable( geom.get() );
	//group->addChild( b.get() );

	// OPTION 2 shape drawable; centered in zero and moved/rotated with transform
	osg::ref_ptr< osg::Box > box = new osg::Box;
	box->setHalfLengths( osg::Vec3( r, .5 * l, r ) );
	box->setRotation( R.getRotate() );
	box->setCenter( half );
    osg::ref_ptr< osg::Geode > b = new osg::Geode;
	b->addDrawable( new osg::ShapeDrawable( box.get() ) );
	group->addChild( b.get() );

	osg::StateSet* set = b->getOrCreateStateSet();
	set->setMode(GL_LIGHTING, osg::StateAttribute::OFF);
	set->setAttributeAndModes( new osg::CullFace );

	// shaders
	static osg::ref_ptr< osg::Program > program;
	if( program == 0 )
	{
		program = new osg::Program;
		program->setName( "cylinder" );
		program->addShader( new osg::Shader( osg::Shader::FRAGMENT, BILLBOARD_CYLINDER3FRAG ) );
		program->addShader( new osg::Shader( osg::Shader::VERTEX,   BILLBOARD_CYLINDER3VERT ) );
	}
	set->setAttributeAndModes( program.get(), osg::StateAttribute::ON );
	set->addUniform( new osg::Uniform( "hlength", float( 0.5 * l ) ) );
	set->addUniform( new osg::Uniform( "Q", QT ) );
	set->addUniform( new osg::Uniform( "axis", axis ) );

	//set->addUniform( new osg::Uniform( "center", half ) );

	set->addUniform( new osg::Uniform( "bcolor", gAtomColors[ atomicNum0 ] ) );
	set->addUniform( new osg::Uniform( "tcolor", gAtomColors[ atomicNum1 ] ) );

	return group.release();
}


//----------------------------------------------------------------------------------
#include "pointcylinder.vert.h"
#include "pointcylinder.frag.h"
#include "pointsphere.frag.h"
#include "pointsphere.vert.h"
// Add quads used as ray targets. Best results are obtained with two quads per primitive
// faster results with one quad put in front of the primitive.
osg::Group* CreateBondGroup8( const osg::Vec3d& P0,
							  int atomicNum0,
							  const osg::Vec3d& P1,
							  int atomicNum1,
							  bool atom0,
							  bool atom1 )
{
	const osg::Vec3d half = ( P0 + P1 ) * 0.5;
    const osg::Vec3d thalf = ( half + P1 ) * 0.5;
	const osg::Vec3d bhalf = ( P0 + half ) * 0.5;
	osg::Matrixd R;
    // bond 0
	osg::Matrixd t0 = ComputeSTRMatrix( osg::Vec3d( 0, 1, 0 ), P0, half, R );

    static const double radius = 0.2;
	const double l = ( P1 - P0 ).length();

	osg::Vec3d axis = osg::Vec3d( 0, 1, 0 ) * R;
	axis.normalize();
	osg::Matrixd Q;
	Q.makeIdentity();
	// | a d e g |
	// | d b f h |
	// | e f c i |
	// | g h i j |
	// ax^2 + by^2 + cz^2 + 2dxy +2exz + 2fyz + 2gx + 2hy + 2iz + j = 0;
	Q( 1, 1 ) = 0;
	Q( 3, 3 ) = -radius * radius;
	osg::Matrixd QT = osg::Matrixd::inverse( R ) * Q * R ; // v * T^-1 * Q * T * vt where: T, Q are 4x4 matrices and v is a ROW vector.


    //PROBLEM: need to split the bond into two cylinders and render
	// the two cylinders with two separate points if not the point doesn't get
	// rendered.

	osg::ref_ptr< osg::Vec3Array > vertices = new osg::Vec3Array;
   	vertices->push_back( bhalf );
	vertices->push_back( thalf );
	osg::ref_ptr< osg::Geometry > geom = new osg::Geometry;
	osg::ref_ptr< osg::Vec4Array > colors = new osg::Vec4Array;
	colors->push_back( gAtomColors[ atomicNum0 ] );
	colors->push_back( gAtomColors[ atomicNum1 ] );
	geom->setColorArray( colors.get() );
	geom->setColorBinding( osg::Geometry::BIND_PER_VERTEX );
	geom->setVertexArray( vertices.get() );
	geom->addPrimitiveSet( new osg::DrawArrays( osg::PrimitiveSet::POINTS, 0, 2 ) );
	osg::ref_ptr< osg::Geode > b = new osg::Geode;
	b->addDrawable( geom.get() );
	osg::StateSet* set = b->getOrCreateStateSet();
	set->setMode(GL_LIGHTING, osg::StateAttribute::OFF);
	osg::ref_ptr< osg::Group > group = new osg::Group;
	group->addChild( b.get() );
	//osg::Point *point = new osg::Point();
	//point->setSize( 32.0f );
	//set->setAttributeAndModes(point);
	set->setMode( GL_VERTEX_PROGRAM_POINT_SIZE, osg::StateAttribute::ON );
	//set->setMode( GL_POINT_SMOOTH, osg::StateAttribute::OFF );
	// shaders
	static osg::ref_ptr< osg::Program > program;
	if( program == 0 )
	{
		program = new osg::Program;
		program->setName( "pointcylinder" );
		program->addShader( new osg::Shader( osg::Shader::FRAGMENT, POINTCYLINDERFRAG ) );
		program->addShader( new osg::Shader( osg::Shader::VERTEX,   POINTCYLINDERVERT ) );
	}
	set->setAttributeAndModes( program.get(), osg::StateAttribute::ON );
	set->addUniform( new osg::Uniform( "hlength", float( 0.25 * l ) ) );
	set->addUniform( new osg::Uniform( "Q", QT ) );
	set->addUniform( new osg::Uniform( "axis", axis ) );


	static osg::ref_ptr< osg::Program > aprogram;
	if( aprogram == 0 )
	{
		aprogram = new osg::Program;
		aprogram->setName( "pointsphere" );
		aprogram->addShader( new osg::Shader( osg::Shader::FRAGMENT, POINTSPHEREFRAG ) );
		aprogram->addShader( new osg::Shader( osg::Shader::VERTEX, POINTSPHEREVERT ) );
	}

	// WARNING: CANNOT ADD SINGLE POINT IN GEOMETRY: WILL BE CULLED
	if( atom0 || atom1 )
	{
		osg::ref_ptr< osg::Geode > geode = new osg::Geode;
		osg::StateSet* set = geode->getOrCreateStateSet();
		set->setAttributeAndModes( aprogram.get(), osg::StateAttribute::ON );
		set->addUniform( new osg::Uniform( "radius", float( radius ) ) );
		set->setMode( GL_LIGHTING, osg::StateAttribute::OFF );
		set->setMode( GL_VERTEX_PROGRAM_POINT_SIZE, osg::StateAttribute::ON );
		geom = new osg::Geometry;
		geode->addDrawable( geom.get() );
		osg::ref_ptr< osg::Vec3Array > coord = new osg::Vec3Array;
		coord->push_back( P0 );
		coord->push_back( P1 );
		geom->setVertexArray( coord.get() );
		osg::ref_ptr< osg::Vec4Array > colors = new osg::Vec4Array;
		colors->push_back( gAtomColors[ atomicNum0 ] );
		colors->push_back( gAtomColors[ atomicNum1 ] );
		geom->setColorArray( colors.get() );
		geom->addPrimitiveSet( new osg::DrawArrays( osg::PrimitiveSet::POINTS, 0, 2 ) );
		geom->setColorBinding( osg::Geometry::BIND_PER_VERTEX );
		group->addChild( geode.get() );
	}

	/*if( atom1 )
	{
		osg::ref_ptr< osg::Geode > geode = new osg::Geode;
		osg::StateSet* set = geode->getOrCreateStateSet();
		set->setAttributeAndModes( aprogram.get(), osg::StateAttribute::ON );
		set->addUniform( new osg::Uniform( "radius", float( radius ) ) );
		set->setMode( GL_LIGHTING, osg::StateAttribute::OFF );
		set->setMode( GL_VERTEX_PROGRAM_POINT_SIZE, osg::StateAttribute::ON );
		geom = new osg::Geometry;
		geode->addDrawable( geom.get() );
		osg::ref_ptr< osg::Vec3Array > coord = new osg::Vec3Array;
		coord->push_back( P1 );
		coord->push_back( P0 );
		geom->setVertexArray( coord.get() );
		osg::ref_ptr< osg::Vec4Array > colors = new osg::Vec4Array;
		colors->push_back( gAtomColors[ atomicNum1 ] );
		geom->setColorArray( colors.get() );
		geom->addPrimitiveSet( new osg::DrawArrays( osg::PrimitiveSet::POINTS, 0, 2 ) );
		geom->setColorBinding( osg::Geometry::BIND_OVERALL );
		group->addChild( geode.get() );
	}*/

	return group.release();
}


//------------------------------------------------------------------------------
//#include "cylinder4.vert.h"
//#include "cylinder4.frag.h"
#include <osg/LineWidth>
osg::Group* CreateBondGroupLOD( const osg::Vec3d& P0,
							  int atomicNum0,
							  const osg::Vec3d& P1,
							  int atomicNum1,
							  bool atom0,
							  bool atom1 )
{
    const osg::Vec3d half = ( P0 + P1   ) * 0.5;

	osg::Matrixd R;
    // bond 0
	osg::Matrixd t0 = ComputeSTRMatrix( osg::Vec3d( 0, 1, 0 ), P0, half, R );
    static const double radius = 0.2;
	const double r = radius;
	const double l = ( P1 - P0 ).length();

	osg::Vec3d axis = osg::Vec3d( 0, 1, 0 ) * R;
	axis.normalize();
	osg::Matrixd Q;
	Q.makeIdentity();
	// | a d e g |
	// | d b f h |
	// | e f c i |
	// | g h i j |
	Q( 1, 1 ) = 0;
	Q( 3, 3 ) = -radius * radius;
	osg::Matrixd QT = osg::Matrixd::inverse( R ) * Q * R ;

	osg::ref_ptr< osg::Group > group = new osg::Group;
	// OPTION 1 oriented box: pass center to shader
	//osg::Vec3 c0( -r, -.5*l, -r ), c1( r, -.5*l, -r ), c2( r, .5*l, -r ), c3( -r, .5*l, -r );
	//osg::Vec3 c4( -r, -.5*l, r ), c5( r, -.5*l, r ), c6( r, .5*l, r ), c7( -r, .5*l, r );
	//osg::ref_ptr< osg::Geometry > geom = new osg::Geometry;
	//osg::ref_ptr< osg::Vec3Array > coord = new osg::Vec3Array;
	////back
	//coord->push_back( c0 * R + half ); coord->push_back( c3 * R + half ); coord->push_back( c2 * R + half ); coord->push_back( c1 * R + half );
	////left
	//coord->push_back( c0 * R + half ); coord->push_back( c4 * R + half ); coord->push_back( c7 * R + half ); coord->push_back( c3 * R + half );
	////right
	//coord->push_back( c5 * R + half ); coord->push_back( c1 * R + half ); coord->push_back( c2 * R + half ); coord->push_back( c6 * R + half );
	//// front
	//coord->push_back( c4 * R + half ); coord->push_back( c5 * R + half ); coord->push_back( c6 * R + half ); coord->push_back( c7 * R + half );
	////top
	//coord->push_back( c7 * R + half ); coord->push_back( c6 * R + half ); coord->push_back( c2 * R + half ); coord->push_back( c3 * R + half );
	////bottom
	//coord->push_back( c4 * R + half ); coord->push_back( c0 * R + half ); coord->push_back( c1 * R + half ); coord->push_back( c5 * R + half );
	//geom->setVertexArray( coord.get() );
	//geom->addPrimitiveSet( new osg::DrawArrays( osg::PrimitiveSet::QUADS, 0, 24 ) );
	//osg::ref_ptr< osg::Geode > b = new osg::Geode;
	//b->addDrawable( geom.get() );
	//group->addChild( b.get() );

	// OPTION 2 shape drawable; centered in zero and moved/rotated with transform
	osg::ref_ptr< osg::Box > box = new osg::Box;
	box->setHalfLengths( osg::Vec3( r, .5 * l, r ) );
	box->setRotation( R.getRotate() );
	box->setCenter( half );
    osg::ref_ptr< osg::Geode > b = new osg::Geode;
	b->addDrawable( new osg::ShapeDrawable( box.get() ) );
	group->addChild( b.get() );

	osg::StateSet* set = b->getOrCreateStateSet();
	set->setMode( GL_LIGHTING, osg::StateAttribute::OFF );
	set->setAttributeAndModes( new osg::CullFace );

	// shaders
	static osg::ref_ptr< osg::Program > program;
	if( program == 0 )
	{
		program = new osg::Program;
		program->setName( "cylinder" );
		program->addShader( new osg::Shader( osg::Shader::FRAGMENT, BILLBOARD_CYLINDER3FRAG ) );
		program->addShader( new osg::Shader( osg::Shader::VERTEX,   BILLBOARD_CYLINDER3VERT ) );
	}
	set->setAttributeAndModes( program.get(), osg::StateAttribute::ON );
	set->addUniform( new osg::Uniform( "hlength", float( 0.5 * l ) ) );
	set->addUniform( new osg::Uniform( "Q", QT ) );
	set->addUniform( new osg::Uniform( "axis", axis ) );

	//set->addUniform( new osg::Uniform( "center", half ) );

	set->addUniform( new osg::Uniform( "bcolor", gAtomColors[ atomicNum0 ] ) );
	set->addUniform( new osg::Uniform( "tcolor", gAtomColors[ atomicNum1 ] ) );

	//
	osg::ref_ptr< osg::Geometry > geom = new osg::Geometry;
	osg::ref_ptr< osg::Vec3Array > coord = new osg::Vec3Array;
	coord->push_back( P0 ); coord->push_back( half );
	coord->push_back( half ); coord->push_back( P1 );
	geom->setVertexArray( coord.get() );
	geom->addPrimitiveSet( new osg::DrawArrays( osg::PrimitiveSet::LINES, 0, 4 ) );
	osg::ref_ptr< osg::Vec4Array > colors = new osg::Vec4Array;
	colors->push_back( gAtomColors[ atomicNum0 ] ); colors->push_back( gAtomColors[ atomicNum1 ] );
	geom->setColorArray( colors.get() );
	//geom->setColorBinding( osg::Geometry::BIND_PER_PRIMITIVE );
	osg::ref_ptr< osg::Geode > geode = new osg::Geode;
	geode->addDrawable( geom.get() );
	group->addChild( geode.get() );
	osg::StateSet* ss = geode->getOrCreateStateSet();
	ss->setMode( GL_LIGHTING, osg::StateAttribute::OFF );
	//ss->setAttributeAndModes( new osg::LineWidth( 4 ) );

	//group->setCenter( half );
	//group->setRadius( 0.5 * l );
	//group->setRange( 0, 0.0f, 12.0f * l );
	//group->setRange( 1, 12.0f * l, 1e7f );
	return group.release();
}
//------------------------------------------------------------------------------
void InitAtomColors( const char* fname )
{
	gAtomColors.clear();
	std::ifstream is( fname );
	if( !is )
	{
		const double* b = GetAtomColorsBegin();
		const double* e = GetAtomColorsEnd();
		b += 3;// skip atom 0
		for( ; b < e; b += 3 )
		{
			gAtomColors.push_back( osg::Vec4d( b[ 0 ], b[ 1 ], b[ 2 ], 2.0 ) * 0.5 +
				                   osg::Vec4d( .5, .5, .5, 0. ) );
		}
		return;
	}
	std::string lineBuffer;
	std::getline( is, lineBuffer );
	while( std::getline( is, lineBuffer ) )
    {
        double r, g, b;
		std::istringstream in( lineBuffer );
        in >> r >> g >> b;
		gAtomColors.push_back( osg::Vec4d( r, g, b, 1.0 ) );
        lineBuffer = "";
    }
}

//------------------------------------------------------------------------------
void InitAtomGeometries( int slices, int stacks, float scaling )
{
	gAtoms.clear();
	const int sz = GetElementTableSize();
	gAtoms.reserve( sz );
	const MolekelElement* et = GetElementTable();
	for( int i = 0; i != sz; ++i )
	{
		osg::Geode* g = CreateSphereNode( scaling * et[ i ].vdwRadius, slices, stacks );
		if( gAtomColors.size() )
		{
			osg::StateSet* ss = g->getOrCreateStateSet();
			osg::Material* m = new osg::Material;
			m->setDiffuse( osg::Material::FRONT_AND_BACK, gAtomColors[ std::min( i, int( gAtomColors.size() - 1 ) ) ] );
			//m->setSpecular( osg::Material::FRONT_AND_BACK, osg::Vec4d( 1, 1, 1, 1 ) );
			//m->setShininess( osg::Material::FRONT_AND_BACK, 20 );
			ss->setAttribute( m, osg::StateAttribute::ON );
			osg::CullFace* cull_face = new osg::CullFace;
			cull_face->setMode( osg::CullFace::BACK );
			ss->setAttribute( cull_face, osg::StateAttribute::ON );
		}
		gAtoms.push_back( g );
	}
}

//------------------------------------------------------------------------------
void InitAtomDiscGeometries( int slices )
{
	gAtoms.clear();
	const int sz = GetElementTableSize();
	gAtoms.reserve( sz );
	const MolekelElement* et = GetElementTable();
	for( int i = 0; i != sz; ++i )
	{
		osg::Billboard* g = new osg::Billboard;
		osg::Geometry* geom = CreateDiscGeometry( et[ i ].vdwRadius, slices );
		g->addDrawable( geom );
		if( gAtomColors.size() )
		{
			osg::StateSet* ss = g->getOrCreateStateSet();
			osg::Material* m = new osg::Material;
			m->setDiffuse( osg::Material::FRONT, gAtomColors[ i ] );
			ss->setAttribute( m, osg::StateAttribute::ON );
		}
		gAtoms.push_back( g );
	}
}

//------------------------------------------------------------------------------
void InitAtomTextureGeometries()
{
	osg::Image* img = osgDB::readImageFile( "./atom.rgb" );
	if( !img ) osg::notify( osg::WARN ) << "Error: cannot read from texture file" << std::endl;
	osg::ref_ptr< osg::Texture2D > texture = new osg::Texture2D;
	texture->setImage( img );
	/*texture->setFilter( osg::Texture::MIN_FILTER, osg::Texture::NEAREST_MIPMAP_NEAREST );
	texture->setFilter( osg::Texture::MAG_FILTER, osg::Texture::NEAREST_MIPMAP_NEAREST );
	texture->setWrap( osg::Texture::WRAP_S, osg::Texture::CLAMP );
	texture->setWrap( osg::Texture::WRAP_T, osg::Texture::CLAMP );*/
	//texture->setUseHardwareMipMapGeneration( true );

	//texture->setDataVariance( osg::Object::DYNAMIC );
	gAtoms.clear();
	const int sz = GetElementTableSize();
	gAtoms.reserve( sz );
	const MolekelElement* et = GetElementTable();
	for( int i = 0; i != sz; ++i )
	{
		osg::Billboard* g = new osg::Billboard;
		osg::Geometry* geom = new osg::Geometry;
		osg::Vec3Array* coord = new osg::Vec3Array;
		const double r = et[ i ].vdwRadius;
		coord->push_back( osg::Vec3d( -r, 0, -r ) );
		coord->push_back( osg::Vec3d( r, 0, -r ) );
		coord->push_back( osg::Vec3d( r, 0, r ) );
		coord->push_back( osg::Vec3d( -r, 0, r ) );
		geom->setVertexArray( coord );
		osg::Vec2Array* tex = new osg::Vec2Array;
		tex->push_back( osg::Vec2d( 0.005, 0.005 ) );
		tex->push_back( osg::Vec2d( 0.995, 0.005 ) );
		tex->push_back( osg::Vec2d( 0.995, 0.995 ) );
		tex->push_back( osg::Vec2d( 0.005, 0.995 ) );

		geom->setTexCoordArray( 0, tex );
		geom->addPrimitiveSet( new osg::DrawArrays( osg::PrimitiveSet::QUADS, 0, 4 ) );

		g->addDrawable( geom );
		if( gAtomColors.size() )
		{
			osg::Vec4Array* colorArray = new osg::Vec4Array;
			colorArray->push_back( gAtomColors[ std::min( i, int( gAtomColors.size() - 1 ) ) ] );
			// An index array for assigning vertices to colors (based on index in the array)
			//osg::TemplateIndexArray< GLuint, osg::Array::UIntArrayType, 4, 1 >* colorIndexArray;
			//colorIndexArray =
			//	new osg::TemplateIndexArray< GLuint, osg::Array::UIntArrayType, 4, 1 >;
			//colorIndexArray->push_back( 0 );
			geom->setColorArray( colorArray );
			//geom->setColorIndices( colorIndexArray );
			geom->setColorBinding( osg::Geometry::BIND_OVERALL );

			osg::StateSet* ss = g->getOrCreateStateSet();
			ss->setTextureAttributeAndModes( 0, texture.get(), osg::StateAttribute::ON );

			osg::AlphaFunc* alphaFunc = new osg::AlphaFunc;
			alphaFunc->setFunction( osg::AlphaFunc::GEQUAL, 0.05f );
			ss->setAttributeAndModes( alphaFunc, osg::StateAttribute::ON );
			ss->setMode( GL_LIGHTING, osg::StateAttribute::OFF );
			//ss->setAttributeAndModes( new osg::BlendFunc, osg::StateAttribute::ON );

		}
		gAtoms.push_back( g );
	}
}

//------------------------------------------------------------------------------
osg::Group* CreateMoleculeGeometry( const std::string& fname,
								    const std::string& format,
									MoleculeDisplayStyle ds,
									bool billboard,
									bool texture )
{
	const bool bonds = ds != MOLDISPSTYLE_SPACEFILL;
	const bool atoms = ds != MOLDISPSTYLE_WIREFRAME && ds != MOLDISPSTYLE_LIQUORICE;
	const bool liquorice = ds == MOLDISPSTYLE_LIQUORICE;
	const float liqRadius = 0.2f;
	OBConversion obConversion;
	if( !bonds ) obConversion.AddOption( "b", OBConversion::INOPTIONS );
    obConversion.SetInFormat( format.c_str() );
	OBMol molecule;
	OBMol* mol = &molecule;
	std::ifstream in( fname.c_str() );
	if( !in ) return 0;
	bool ok = obConversion.Read( &molecule, &in );
	if( !ok ) return 0;

	osg::ref_ptr< osg::Group > group = new osg::Group;
	const int numAtoms = mol->NumAtoms();
	//iterate over atoms and add one transform and one atom
	if( atoms )
	{
		for( int a = 0; a != numAtoms; ++a )
	    {
			osg::PositionAttitudeTransform* t = new osg::PositionAttitudeTransform;
			OBAtom* atom = mol->GetAtom( a + 1 );
			osg::Vec3d p( atom->GetX(), atom->GetY(), atom->GetZ() );
			t->setPosition( p );

			osg::ref_ptr< osg::Geode > g = gAtoms[ atom->GetAtomicNum() ];
			osg::Geometry* gg = new osg::Geometry;
			osg::Geometry* geom = gAtoms[ atom->GetAtomicNum() ]->getDrawable( 0 )->asGeometry();
			if( billboard )
			{
				osg::ref_ptr< osg::Billboard > b = new osg::Billboard;
				b->setMode( osg::Billboard::POINT_ROT_EYE );
				b->addDrawable( gg, p );
				g = b.get();
			}
			else
			{
				g = new osg::Geode;
				g->addDrawable( gg );
			}

			gg->setVertexArray( geom->getVertexArray() );
			if( billboard && texture )
			{
				gg->setTexCoordArray( 0, geom->getTexCoordArray( 0 ) );
				gg->setColorArray( geom->getColorArray() );
				gg->setColorBinding( osg::Geometry::BIND_OVERALL );

			}
			else
			{
				gg->setNormalArray( geom->getNormalArray() );
				gg->setNormalBinding( osg::Geometry::BIND_PER_VERTEX );
			}

			g->setStateSet( gAtoms[ atom->GetAtomicNum() ]->getStateSet() );
			for ( GLuint ipr=0; ipr < geom->getNumPrimitiveSets(); ipr++ )
			{
				osg::PrimitiveSet* prs = geom->getPrimitiveSet( ipr );
				gg->addPrimitiveSet( prs );
			}

			if( !billboard )
			{
				t->addChild( g.get() );
				group->addChild( t );
			}
			else group->addChild( g.get() );
		}
	}

	if( bonds )
	{
		std::set< OBAtom* > bonded;
		const int numBonds = mol->NumBonds();
		for( int b = 0; b != numBonds; ++b )
		{
			OBBond* bond = mol->GetBond( b );
			if( !bond ) continue;
			OBAtom* begin = bond->GetBeginAtom();
			OBAtom* end = bond->GetEndAtom();
			const bool a0 = liquorice && bonded.find( begin ) == bonded.end();
			const bool a1 = liquorice && bonded.find( end ) == bonded.end();
			osg::Vec3d p0( begin->GetX(), begin->GetY(), begin->GetZ() );
			osg::Vec3d p1( end->GetX(), end->GetY(), end->GetZ() );
			group->addChild( CreateBondGroup( p0, begin->GetAtomicNum(), p1, end->GetAtomicNum(), a0, a1 ) );
			if( liquorice )
			{
				bonded.insert( begin );
				bonded.insert( end   );
			}
		}
	}

	return group.release();
}


//------------------------------------------------------------------------------
osg::Group* CreateBillboardMoleculeGeometry( const std::string& fname,
											 const std::string& format,
											 MoleculeDisplayStyle ds )
{
	const bool bonds = ds != MOLDISPSTYLE_SPACEFILL;
	const bool atoms = ( ds != MOLDISPSTYLE_WIREFRAME ) && ( ds != MOLDISPSTYLE_LIQUORICE );
	const bool liquorice = ds == MOLDISPSTYLE_LIQUORICE;
	const float liqRadius = 0.2f;
	OBConversion obConversion;
	if( !bonds ) obConversion.AddOption( "b", OBConversion::INOPTIONS );
    obConversion.SetInFormat( format.c_str() );
	OBMol molecule;
	OBMol* mol = &molecule;
	std::ifstream in( fname.c_str() );
	if( !in ) return 0;
	bool ok = obConversion.Read( &molecule, &in );
	if( !ok ) return 0;

	osg::ref_ptr< osg::Group > group = new osg::Group;
	const int numAtoms = mol->NumAtoms();
	//iterate over atoms and add one transform and one atom
	if( atoms )
	{
		for( int a = 0; a != numAtoms; ++a )
	    {
			OBAtom* atom = mol->GetAtom( a + 1 );
			osg::Vec3d p( atom->GetX(), atom->GetY(), atom->GetZ() );
			osg::ref_ptr< osg::Billboard > b = new osg::Billboard;
			b->setMode( osg::Billboard::POINT_ROT_EYE );
			osg::Geometry* gg = new osg::Geometry;
			osg::Geometry* geom = gAtoms[ atom->GetAtomicNum() ]->getDrawable( 0 )->asGeometry();
			b->addDrawable( gg, p );
			gg->setVertexArray( geom->getVertexArray() );
			gg->setColorArray( geom->getColorArray() );
			gg->setTexCoordArray( 0, geom->getTexCoordArray( 0 ) );
			gg->setColorBinding( osg::Geometry::BIND_PER_VERTEX );
			//gg->setNormalArray( geom->getNormalArray() );
			//gg->setNormalBinding( osg::Geometry::BIND_PER_VERTEX );
			b->setStateSet( gAtoms[ atom->GetAtomicNum() ]->getStateSet() );
			for ( GLuint ipr=0; ipr < geom->getNumPrimitiveSets(); ipr++ )
			{
				osg::PrimitiveSet* prs = geom->getPrimitiveSet( ipr );
				gg->addPrimitiveSet( prs );
			}
			group->addChild( b.get() );
		}
	}

	if( bonds )
	{
		std::set< OBAtom* > bonded;
		const int numBonds = mol->NumBonds();
		for( int b = 0; b != numBonds; ++b )
		{
			OBBond* bond = mol->GetBond( b );
			if( !bond ) continue;
			OBAtom* begin = bond->GetBeginAtom();
			OBAtom* end = bond->GetEndAtom();
			const bool a0 = liquorice && bonded.find( begin ) == bonded.end();
			const bool a1 = liquorice && bonded.find( end ) == bonded.end();
			osg::Vec3d p0( begin->GetX(), begin->GetY(), begin->GetZ() );
			osg::Vec3d p1( end->GetX(), end->GetY(), end->GetZ() );
			group->addChild( CreateBondGroup2( p0, begin->GetAtomicNum(), p1, end->GetAtomicNum(), a0, a1 ) );
			if( liquorice )
			{
				bonded.insert( begin );
				bonded.insert( end   );
			}
		}
	}

	return group.release();
}
//------------------------------------------------------------------------------
osg::Group* CreateMoleculeGeometry2( const std::string& fname,
								     const std::string& format,
									 bool billboard,
								 	 bool texture )
{

	OBConversion obConversion;
	obConversion.AddOption( "b", OBConversion::INOPTIONS );
    obConversion.SetInFormat( format.c_str() );
	OBMol* mol = new OBMol;
	std::ifstream in( fname.c_str() );
	if( !in ) return 0;
	bool ok = obConversion.Read( mol, &in );
	if( !ok ) return 0;


	osg::ref_ptr< osg::Billboard > geode = new osg::Billboard;
	geode->setMode( osg::Billboard::POINT_ROT_EYE );
	geode->setStateSet( gAtoms[ 6 ]->getStateSet() );
	const int numAtoms = mol->NumAtoms();
	for( int a = 0; a != numAtoms; ++a )
    {
		OBAtom* atom = mol->GetAtom( a + 1 );
		osg::Vec3d p( atom->GetX(), atom->GetY(), atom->GetZ() );
		osg::Geometry* gg = new osg::Geometry;
		//gg->setUseVertexBufferObjects( true );
		osg::Geometry* geom = gAtoms[ atom->GetAtomicNum() ]->getDrawable( 0 )->asGeometry();
		geode->addDrawable( gg, p );
		gg->setVertexArray( geom->getVertexArray() );
		gg->setTexCoordArray( 0, geom->getTexCoordArray( 0 ) );
		gg->setColorArray( geom->getColorArray() );
		gg->setColorBinding( osg::Geometry::BIND_OVERALL );
		for ( GLuint ipr=0; ipr < geom->getNumPrimitiveSets(); ipr++ )
		{
			osg::PrimitiveSet* prs = geom->getPrimitiveSet( ipr );
			gg->addPrimitiveSet( prs );
		}
	}
	osg::ref_ptr< osg::Group > group = new osg::Group;
	group->addChild( geode.get() );
	return group.release();
}


//------------------------------------------------------------------------------
osg::StateSet* makeStateSet(float size)
{
    osg::StateSet *set = new osg::StateSet();

    /// Setup cool blending
    set->setMode(GL_BLEND, osg::StateAttribute::ON);
    //osg::BlendFunc *fn = new osg::BlendFunc();
    //fn->setFunction(osg::BlendFunc::SRC_ALPHA, osg::BlendFunc::DST_ALPHA);
    //set->setAttributeAndModes(fn, osg::StateAttribute::ON);
	osg::AlphaFunc* alphaFunc = new osg::AlphaFunc;
	alphaFunc->setFunction( osg::AlphaFunc::GEQUAL, 0.05f );
	set->setAttributeAndModes( alphaFunc, osg::StateAttribute::ON );

    /// Setup the point sprites
    osg::PointSprite *sprite = new osg::PointSprite();
    set->setTextureAttributeAndModes(0, sprite, osg::StateAttribute::ON);

    /// Give some size to the points to be able to see the sprite
    osg::Point *point = new osg::Point();
    point->setSize(size);
    point->setDistanceAttenuation( osg::Vec3( 0, 1 ,0 ) ); // 1/d
    set->setAttribute(point);

    /// Disable depth test to avoid sort problems and Lighting
    //set->setMode(GL_DEPTH_TEST, osg::StateAttribute::OFF);
    set->setMode(GL_LIGHTING, osg::StateAttribute::OFF);

    /// The texture for the sprites
    osg::Texture2D *tex = new osg::Texture2D();
    tex->setImage(osgDB::readImageFile("./atom.rgb"));
    set->setTextureAttributeAndModes(0, tex, osg::StateAttribute::ON);

    return set;
}

//------------------------------------------------------------------------------
osg::Group* CreatePointSpriteMoleculeGeometry( const std::string& fname,
											   const std::string& format,
											   MoleculeDisplayStyle ds )
{

	OBConversion obConversion;
	obConversion.AddOption( "b", OBConversion::INOPTIONS );
    obConversion.SetInFormat( format.c_str() );
	OBMol* mol = new OBMol;
	std::ifstream in( fname.c_str() );
	if( !in ) return 0;
	bool ok = obConversion.Read( mol, &in );
	if( !ok ) return 0;

	const int numAtoms = mol->NumAtoms();
	const float pointsize = 16;
	typedef std::map< int, osg::ref_ptr< osg::Geometry > > GEOMETRIES;
	GEOMETRIES geometries;
	const MolekelElement* et = GetElementTable();
	for( int a = 0; a != numAtoms; ++a )
    {
		OBAtom* atom = mol->GetAtom( a + 1 );
		osg::Vec3d p( atom->GetX(), atom->GetY(), atom->GetZ() );
		osg::Geometry* geom = gAtoms[ atom->GetAtomicNum() ]->getDrawable( 0 )->asGeometry();
		if( geometries.find( atom->GetAtomicNum() ) == geometries.end() )
		{
			geometries[ atom->GetAtomicNum() ] = new osg::Geometry;
			geometries[ atom->GetAtomicNum() ]->setVertexArray( new osg::Vec3Array );
			osg::ref_ptr< osg::Vec4Array > colors = new osg::Vec4Array;
			colors->push_back( gAtomColors[ atom->GetAtomicNum() ] );
			geometries[ atom->GetAtomicNum() ]->setColorArray( colors.get() );
			geometries[ atom->GetAtomicNum() ]->setColorBinding(osg::Geometry::BIND_OVERALL );
		}
		osg::Vec3Array* coords = static_cast< osg::Vec3Array* >( geometries[ atom->GetAtomicNum() ]->getVertexArray() );
		coords->push_back( p );
	}

	osg::ref_ptr< osg::Group > group = new osg::Group;
	for( GEOMETRIES::iterator i = geometries.begin(); i != geometries.end(); ++i )
	{
		osg::ref_ptr< osg::Geode > g = new osg::Geode;
		g->setStateSet( makeStateSet( et[ i->first ].vdwRadius * pointsize ) );
		i->second->addPrimitiveSet(
			new osg::DrawArrays( osg::PrimitiveSet::POINTS, 0, static_cast< osg::Vec3Array* >( i->second->getVertexArray() )->size() ) );
		g->addDrawable( i->second.get() );
		group->addChild( g.get() );
	}

	return group.release();
}



//------------------------------------------------------------------------------

static const char gVertexSource[] =
    "varying vec4 color;\n"
    "uniform float radius;\n"
    "void main(void)\n"
    "{\n"
        "    gl_Position = gl_ModelViewProjectionMatrix * gl_Vertex;\n"
	"    gl_PointSize = 100.0;\n"
    "}\n";

osg::StateSet* makeShaderStateSet( float size )
{
    osg::StateSet *set = new osg::StateSet();

    /// Setup cool blending
    set->setMode( GL_BLEND, osg::StateAttribute::ON );

	/// Add shaders
	set->setMode( GL_VERTEX_PROGRAM_POINT_SIZE, osg::StateAttribute::ON );
	osg::Program* program = new osg::Program;
    program->setName( "pointshader" );
    program->addShader( new osg::Shader( osg::Shader::VERTEX, gVertexSource ) );
    //program->addShader( new osg::Shader( osg::Shader::FRAGMENT, gFragSource ) );
	set->setAttributeAndModes( program, osg::StateAttribute::ON );
	set->addUniform( new osg::Uniform( "radius", size ) );

    //osg::BlendFunc *fn = new osg::BlendFunc();
    //fn->setFunction(osg::BlendFunc::SRC_ALPHA, osg::BlendFunc::DST_ALPHA);
    //set->setAttributeAndModes(fn, osg::StateAttribute::ON);
	osg::AlphaFunc* alphaFunc = new osg::AlphaFunc;
	alphaFunc->setFunction( osg::AlphaFunc::GEQUAL, 0.05f );
	set->setAttributeAndModes( alphaFunc, osg::StateAttribute::ON );

    /// Setup the point sprites
    osg::PointSprite *sprite = new osg::PointSprite();
    set->setTextureAttributeAndModes(0, sprite, osg::StateAttribute::ON);

    /// Give some size to the points to be able to see the sprite
    osg::Point *point = new osg::Point();
    //point->setSize(size);
    set->setAttribute(point);

    /// Disable depth test to avoid sort problems and Lighting
    //set->setMode(GL_DEPTH_TEST, osg::StateAttribute::OFF);
    set->setMode(GL_LIGHTING, osg::StateAttribute::OFF);

    /// The texture for the sprites
    osg::Texture2D *tex = new osg::Texture2D();
    tex->setImage(osgDB::readImageFile("./atom.png"));
    set->setTextureAttributeAndModes(0, tex, osg::StateAttribute::ON);

    return set;
}


osg::Group* CreateShaderPointSpriteMoleculeGeometry( const std::string& fname,
											         const std::string& format )
{

	OBConversion obConversion;
	obConversion.AddOption( "b", OBConversion::INOPTIONS );
    obConversion.SetInFormat( format.c_str() );
	OBMol* mol = new OBMol;
	std::ifstream in( fname.c_str() );
	if( !in ) return 0;
	bool ok = obConversion.Read( mol, &in );
	if( !ok ) return 0;

	const int numAtoms = mol->NumAtoms();
	osg::Geometry* gg = new osg::Geometry;
	osg::ref_ptr< osg::Geode > g = new osg::Geode;
	g->addDrawable( gg );
	osg::ref_ptr< osg::Vec3Array > centers = new osg::Vec3Array;
	osg::ref_ptr< osg::Vec4Array > colors = new osg::Vec4Array;
	gg->setVertexArray( centers.get() );
	gg->setColorArray( colors.get() );
	gg->setColorBinding(osg::Geometry::BIND_PER_VERTEX);
    gg->addPrimitiveSet(new osg::DrawArrays( osg::PrimitiveSet::POINTS, 0, numAtoms ));
	for( int a = 0; a != numAtoms; ++a )
    {
		OBAtom* atom = mol->GetAtom( a + 1 );
		osg::Vec3d p( atom->GetX(), atom->GetY(), atom->GetZ() );
		osg::Geometry* geom = gAtoms[ atom->GetAtomicNum() ]->getDrawable( 0 )->asGeometry();
		centers->push_back( p );
		colors->push_back( gAtomColors[ atom->GetAtomicNum() ] );
	}
	g->setStateSet( makeShaderStateSet( 100.0f ) );
	osg::ref_ptr< osg::Group > group = new osg::Group;
	group->addChild( g.get() );
	return group.release();
}




//------------------------------------------------------------------------------
#include "quadcylinder.vert.h"
#include "quadcylinder.frag.h"
osg::Group* CreateRayTracedMoleculeGeometry( const std::string& fname,
						 				     const std::string& format,
											 MoleculeDisplayStyle ds,
											 float radScale )
{

	const bool bonds = ds != MOLDISPSTYLE_SPACEFILL;
	const bool atoms = ds != MOLDISPSTYLE_WIREFRAME && ds != MOLDISPSTYLE_LIQUORICE;
	const bool liquorice = ds == MOLDISPSTYLE_LIQUORICE;
	const float liqRadius = 0.2f;
	OBConversion obConversion;
	if( !bonds ) obConversion.AddOption( "b", OBConversion::INOPTIONS );
    obConversion.SetInFormat( format.c_str() );
	OBMol* mol = new OBMol;
	std::ifstream in( fname.c_str() );
	if( !in ) return 0;
	bool ok = obConversion.Read( mol, &in );
	if( !ok ) return 0;

	const int numAtoms = mol->NumAtoms();
	typedef std::map< int, osg::ref_ptr< osg::Billboard > > BILLBOARDS;
	BILLBOARDS billboards;
	const MolekelElement* et = GetElementTable();

	static osg::ref_ptr< osg::Program > aprogram;
	if( aprogram == 0 )
	{
		aprogram = new osg::Program;
		aprogram->setName( "sphere" );
		//aprogram->addShader( osg::Shader::readShaderFile( osg::Shader::FRAGMENT, "sphere3.frag" ) );
		//aprogram->addShader( osg::Shader::readShaderFile( osg::Shader::VERTEX, "sphere3.vert" ) );
		aprogram->addShader( new osg::Shader( osg::Shader::FRAGMENT, SPHEREFRAG ) );
		aprogram->addShader( new osg::Shader( osg::Shader::VERTEX, SPHEREVERT ) );
	}

	osg::ref_ptr< osg::Group > group = new osg::Group;

	for( int a = 0; a != numAtoms && atoms; ++a )
    {
		OBAtom* atom = mol->GetAtom( a + 1 );
		osg::Vec3d p( atom->GetX(), atom->GetY(), atom->GetZ() );
		osg::ref_ptr< osg::Geometry > geom;
		if( billboards.find( atom->GetAtomicNum() ) == billboards.end() )
		{
			billboards[ atom->GetAtomicNum() ] = new osg::Billboard;
		    billboards[ atom->GetAtomicNum() ]->setMode( osg::Billboard::POINT_ROT_EYE );
			osg::ref_ptr< osg::Vec4Array > colors = new osg::Vec4Array;
			colors->push_back( gAtomColors[ atom->GetAtomicNum() ] );

			osg::StateSet* set = billboards[ atom->GetAtomicNum() ]->getOrCreateStateSet();
			set->setAttributeAndModes( aprogram.get(), osg::StateAttribute::ON );
			float radius = radScale * et[ atom->GetAtomicNum() ].vdwRadius;
			set->addUniform( new osg::Uniform( "radius", radius ) );
			set->setMode(GL_LIGHTING, osg::StateAttribute::OFF);

			geom = new osg::Geometry;
			if( gAtomColors.size() )
			{
				osg::ref_ptr< osg::Vec4Array > colors = new osg::Vec4Array;
				colors->push_back( gAtomColors[ std::min( size_t( atom->GetAtomicNum() ), gAtomColors.size() - 1 ) ] );
				geom->setColorArray( colors.get() );
				geom->setColorBinding( osg::Geometry::BIND_OVERALL );
			}
			osg::ref_ptr< osg::Vec3Array > coord = new osg::Vec3Array;
			// put quad in front: billboards are automatically rotated such that -y is parallel to z in eye
			// coordinates
			coord->push_back( osg::Vec3d( -radius, -radius, -radius ) );
			coord->push_back( osg::Vec3d(  radius, -radius, -radius ) );
			coord->push_back( osg::Vec3d(  radius, -radius,  radius ) );
			coord->push_back( osg::Vec3d( -radius, -radius,  radius ) );
			geom->setVertexArray( coord.get() );
			geom->addPrimitiveSet( new osg::DrawArrays( osg::PrimitiveSet::QUADS, 0, 4 ) );
			group->addChild( billboards[ atom->GetAtomicNum() ].get() );
		}
		else geom = billboards[ atom->GetAtomicNum() ]->getDrawable( 0 )->asGeometry();
		billboards[ atom->GetAtomicNum() ]->addDrawable( geom.get(), p );
	}

	// BONDS
	//if( bonds )
	//{
	//	std::set< OBAtom* > bonded;
	//	const int numBonds = mol->NumBonds();
	//	for( int b = 0; b != numBonds; ++b )
	//	{
	//		OBBond* bond = mol->GetBond( b );
	//		if( !bond ) continue;
	//		OBAtom* begin = bond->GetBeginAtom();
	//		OBAtom* end = bond->GetEndAtom();
	//		const bool a0 = liquorice && bonded.find( begin ) == bonded.end();
	//		const bool a1 = liquorice && bonded.find( end )   == bonded.end();
	//		osg::Vec3d p0( begin->GetX(), begin->GetY(), begin->GetZ() );
	//		osg::Vec3d p1( end->GetX(), end->GetY(), end->GetZ() );
	//		group->addChild( CreateBondGroup6( p0, begin->GetAtomicNum(), p1, end->GetAtomicNum(), a0, a1 ) );
	//		if( liquorice )
	//		{
	//			bonded.insert( begin );
	//			bonded.insert( end   );
	//		}
	//	}
	//}
	// BONDS
	if( bonds )
	{
		std::set< OBAtom* > bonded;
		const int numBonds = mol->NumBonds();
		osg::ref_ptr< osg::Billboard > geode = new osg::Billboard;
		geode->setMode( osg::Billboard::POINT_ROT_EYE );

		static osg::ref_ptr< osg::Program > bprogram;
		if( bprogram == 0 )
		{
			bprogram = new osg::Program;
			bprogram->setName( "quadcylinder" );
			bprogram->addShader( new osg::Shader( osg::Shader::FRAGMENT, QUADCYLINDERFRAG ) );
			bprogram->addShader( new osg::Shader( osg::Shader::VERTEX, QUADCYLINDERVERT ) );
		}
		osg::StateSet* set = geode->getOrCreateStateSet();
		set->setMode( GL_CULL_FACE, osg::StateAttribute::OFF );
		set->setAttributeAndModes( bprogram.get(), osg::StateAttribute::ON );
		set->addUniform( new osg::Uniform( "radius", liqRadius ) );
		set->setMode(GL_LIGHTING, osg::StateAttribute::OFF);
		group->addChild( geode.get() );

		for( int b = 0; b != numBonds; ++b )
		{
			OBBond* bond = mol->GetBond( b );
			if( !bond ) continue;
			OBAtom* begin = bond->GetBeginAtom();
			OBAtom* end = bond->GetEndAtom();
			const bool a0 = liquorice && bonded.find( begin ) == bonded.end();
			const bool a1 = liquorice && bonded.find( end ) == bonded.end();
			osg::Vec3d p0( begin->GetX(), begin->GetY(), begin->GetZ() );
			osg::Vec3d p1( end->GetX(), end->GetY(), end->GetZ() );
			const osg::Vec3d half = ( p0 + p1 ) * 0.5;
			const osg::Vec3d bh = ( p0 + half ) * 0.5;
			const osg::Vec3d th = ( half + p1 ) * 0.5;
			const osg::Vec3d bdir = half - p0;
			const osg::Vec3d tdir = p1 - half;

			const double l = bdir.length();
			const double r = 0.2;
			osg::Matrixd R;
			osg::Vec3d d = bdir;
			d.normalize();
			R.makeRotate( osg::Vec3d( 0., 1., 0. ), d );
			osg::BoundingBox bbox;
			bbox.expandBy( osg::Vec3d( -r, -.5*l, -r ) * R );
			bbox.expandBy( osg::Vec3d(  r, -.5*l, -r ) * R );
			bbox.expandBy( osg::Vec3d(  r,  .5*l, -r ) * R );
			bbox.expandBy( osg::Vec3d( -r,  .5*l, -r ) * R );
			bbox.expandBy( osg::Vec3d( -r, -.5*l,  r ) * R );
			bbox.expandBy( osg::Vec3d(  r, -.5*l,  r ) * R );
			bbox.expandBy( osg::Vec3d(  r,  .5*l, r ) * R );
			bbox.expandBy( osg::Vec3d( -r,  .5*l, r ) * R );
			const double br = bbox.radius();


			osg::ref_ptr< osg::Geometry >  geom1   = new osg::Geometry;
			osg::ref_ptr< osg::Vec4Array > colors1 = new osg::Vec4Array;
			osg::ref_ptr< osg::Vec3Array > quads1  = new osg::Vec3Array;
			osg::ref_ptr< osg::Vec3Array > dir1	   = new osg::Vec3Array;
			geom1->setVertexArray( quads1.get() );
			geom1->setColorArray( colors1.get() );
			geom1->setNormalArray( dir1.get() );
			geom1->setNormalBinding( osg::Geometry::BIND_PER_VERTEX );
			geom1->setColorBinding( osg::Geometry::BIND_PER_VERTEX );
			geom1->addPrimitiveSet( new osg::DrawArrays( osg::PrimitiveSet::QUADS, 0, 4 ) );
			quads1->push_back( osg::Vec3d( -br , -br, -br ) );
			quads1->push_back( osg::Vec3d(  br , -br, -br ) );
			quads1->push_back( osg::Vec3d(  br , -br,  br ) );
			quads1->push_back( osg::Vec3d( -br , -br,  br ) );
			dir1->push_back( bdir );
			dir1->push_back( bdir );
			dir1->push_back( bdir );
			dir1->push_back( bdir );
			colors1->push_back( gAtomColors[ begin->GetAtomicNum() ] );
			colors1->push_back( gAtomColors[ begin->GetAtomicNum() ] );
			colors1->push_back( gAtomColors[ begin->GetAtomicNum() ] );
			colors1->push_back( gAtomColors[ begin->GetAtomicNum() ] );

			osg::ref_ptr< osg::Geometry >  geom2   = new osg::Geometry;
			osg::ref_ptr< osg::Vec4Array > colors2 = new osg::Vec4Array;
			osg::ref_ptr< osg::Vec3Array > quads2  = new osg::Vec3Array;
			osg::ref_ptr< osg::Vec3Array > dir2	   = new osg::Vec3Array;
			geom2->setVertexArray( quads2.get() );
			geom2->setColorArray( colors2.get() );
			geom2->setNormalArray( dir2.get() );
			geom2->setNormalBinding( osg::Geometry::BIND_PER_VERTEX );
			geom2->setColorBinding( osg::Geometry::BIND_PER_VERTEX );
			geom2->addPrimitiveSet( new osg::DrawArrays( osg::PrimitiveSet::QUADS, 0, 4 ) );
			quads2->push_back( osg::Vec3d( -br , -br, -br ) );
			quads2->push_back( osg::Vec3d(  br , -br, -br ) );
			quads2->push_back( osg::Vec3d(  br , -br,  br ) );
			quads2->push_back( osg::Vec3d( -br , -br,  br ) );
			dir2->push_back( tdir );
			dir2->push_back( tdir );
			dir2->push_back( tdir );
			dir2->push_back( tdir );
			colors2->push_back( gAtomColors[ end->GetAtomicNum() ] );
			colors2->push_back( gAtomColors[ end->GetAtomicNum() ] );
			colors2->push_back( gAtomColors[ end->GetAtomicNum() ] );
			colors2->push_back( gAtomColors[ end->GetAtomicNum() ] );

			geode->addDrawable( geom1.get(), bh );
			geode->addDrawable( geom2.get(), th );

			if( liquorice )
			{
				bonded.insert( begin );
				bonded.insert( end   );
			}
		}
	}
	return group.release();
}

//------------------------------------------------------------------------------
// RAY TRACED QUADS
#include "quadcylinder2.vert.h"
#include "quadsphere.vert.h"
#include "quadcylinder2.frag.h"
#include "coolsphere3.frag.h"
osg::Group* CreateRayTracedMoleculeGeometry2( const std::string& fname,
						 				     const std::string& format,
											 MoleculeDisplayStyle ds,
											 float radScale )
{

	const bool bonds = ds != MOLDISPSTYLE_SPACEFILL;
	const bool atoms = ds != MOLDISPSTYLE_WIREFRAME;
	const bool liquorice = ds == MOLDISPSTYLE_LIQUORICE;
	const float liqRadius = radScale;
	OBConversion obConversion;
	if( !bonds ) obConversion.AddOption( "b", OBConversion::INOPTIONS );
    obConversion.SetInFormat( format.c_str() );
	OBMol* mol = new OBMol;
	std::ifstream in( fname.c_str() );
	if( !in ) return 0;
	bool ok = obConversion.Read( mol, &in );
	if( !ok ) return 0;

	const int numAtoms = mol->NumAtoms();
	typedef std::map< int, osg::ref_ptr< osg::Geode > > GEODES;
	GEODES geodes;
	const MolekelElement* et = GetElementTable();
	static osg::ref_ptr< osg::Program > aprogram;
	if( aprogram == 0 )
	{
		aprogram = new osg::Program;
		aprogram->setName( "sphere2" );
		aprogram->addShader( new osg::Shader( osg::Shader::FRAGMENT, COOLSPHEREFRAG ) );
		aprogram->addShader( new osg::Shader( osg::Shader::VERTEX, SPHEREVERT2 ) );
	}

	osg::ref_ptr< osg::Group > group = new osg::Group;

	for( int a = 0; a != numAtoms && atoms; ++a )
    {
		OBAtom* atom = mol->GetAtom( a + 1 );
		osg::Vec3d p( atom->GetX(), atom->GetY(), atom->GetZ() );
		float radius = liquorice ? liqRadius : radScale * et[ atom->GetAtomicNum() ].vdwRadius;
		if( geodes.find( atom->GetAtomicNum() ) == geodes.end() )
		{
			geodes[ atom->GetAtomicNum() ] = new osg::Geode;
			osg::StateSet* set = geodes[ atom->GetAtomicNum() ]->getOrCreateStateSet();
			set->setAttributeAndModes( aprogram.get(), osg::StateAttribute::ON );
			set->addUniform( new osg::Uniform( "radius", radius ) );
			set->setMode(GL_LIGHTING, osg::StateAttribute::OFF);

			osg::ref_ptr< osg::Geometry >  geom   = new osg::Geometry;
			osg::ref_ptr< osg::Vec4Array > colors = new osg::Vec4Array;
			osg::ref_ptr< osg::Vec3Array > quads  = new osg::Vec3Array;
			osg::ref_ptr< osg::Vec3Array > pos	  = new osg::Vec3Array;

			geom->setVertexArray( quads.get() );
			geom->setColorArray( colors.get() );
			const int texUnit = 1;
			geom->setTexCoordArray( texUnit, pos.get() );
			geom->setNormalBinding( osg::Geometry::BIND_PER_VERTEX );
			geom->setColorBinding( osg::Geometry::BIND_PER_VERTEX );
			geodes[ atom->GetAtomicNum() ]->addDrawable( geom.get() );

			quads->push_back( p + osg::Vec3d( -radius, -radius, 0. ) );
			quads->push_back( p + osg::Vec3d(  radius, -radius, 0. ) );
			quads->push_back( p + osg::Vec3d(  radius, radius, 0. ) );
			quads->push_back( p + osg::Vec3d( -radius, radius, 0. ) );
			pos->push_back( p );
			pos->push_back( p );
			pos->push_back( p );
			pos->push_back( p );
			colors->push_back( gAtomColors[ atom->GetAtomicNum() ] );
			colors->push_back( gAtomColors[ atom->GetAtomicNum() ] );
			colors->push_back( gAtomColors[ atom->GetAtomicNum() ] );
			colors->push_back( gAtomColors[ atom->GetAtomicNum() ] );

			group->addChild( geodes[ atom->GetAtomicNum() ].get() );
		}
		else
		{
			osg::Geometry* geom = geodes[ atom->GetAtomicNum() ]->getDrawable( 0 )->asGeometry();
			osg::ref_ptr< osg::Vec4Array > colors = static_cast< osg::Vec4Array* >( geom->getColorArray()  );
			osg::ref_ptr< osg::Vec3Array > quads  = static_cast< osg::Vec3Array* >( geom->getVertexArray() );
			osg::ref_ptr< osg::Vec3Array > pos	  = static_cast< osg::Vec3Array* >( geom->getTexCoordArray( 1 ) );
			quads->push_back( p + osg::Vec3d( -radius, -radius, 0. ) );
			quads->push_back( p + osg::Vec3d(  radius, -radius, 0. ) );
			quads->push_back( p + osg::Vec3d(  radius, radius, 0. ) );
			quads->push_back( p + osg::Vec3d( -radius, radius, 0. ) );
			pos->push_back( p );
			pos->push_back( p );
			pos->push_back( p );
			pos->push_back( p );
			colors->push_back( gAtomColors[ atom->GetAtomicNum() ] );
			colors->push_back( gAtomColors[ atom->GetAtomicNum() ] );
			colors->push_back( gAtomColors[ atom->GetAtomicNum() ] );
			colors->push_back( gAtomColors[ atom->GetAtomicNum() ] );
		}
	}

	if( atoms )
	{
		for( GEODES::iterator i = geodes.begin(); i != geodes.end(); ++i )
		{
			osg::Geode* g = i->second.get();
			osg::Geometry* geom = g->getDrawable( 0 )->asGeometry();
			geom->addPrimitiveSet( new osg::DrawArrays( osg::PrimitiveSet::QUADS, 0,
														geom->getVertexArray()->getNumElements() ) );
		}
	}

	// BONDS
	if( bonds )
	{
		std::set< OBAtom* > bonded;
		const int numBonds = mol->NumBonds();
		osg::ref_ptr< osg::Geode > geode = new osg::Geode;

		static osg::ref_ptr< osg::Program > bprogram;
		if( bprogram == 0 )
		{
			bprogram = new osg::Program;
			bprogram->setName( "quadcylinder2" );
			bprogram->addShader( new osg::Shader( osg::Shader::FRAGMENT, QUADCYLINDERFRAG2 ) );
			bprogram->addShader( new osg::Shader( osg::Shader::VERTEX, QUADCYLINDERVERT2 ) );
		}
		osg::StateSet* set = geode->getOrCreateStateSet();
		set->setMode( GL_CULL_FACE, osg::StateAttribute::OFF );
		set->setAttributeAndModes( bprogram.get(), osg::StateAttribute::ON );
		set->addUniform( new osg::Uniform( "radius", liqRadius ) );
		set->setMode(GL_LIGHTING, osg::StateAttribute::OFF);
		group->addChild( geode.get() );

		osg::ref_ptr< osg::Geometry >  geom   = new osg::Geometry;
		osg::ref_ptr< osg::Vec4Array > colors = new osg::Vec4Array;
		osg::ref_ptr< osg::Vec3Array > quads  = new osg::Vec3Array;
		osg::ref_ptr< osg::Vec3Array > dir	  = new osg::Vec3Array;
		osg::ref_ptr< osg::Vec3Array > pos	  = new osg::Vec3Array;
		geom->setVertexArray( quads.get() );
		geom->setColorArray( colors.get() );
		geom->setNormalArray( dir.get() );
		const int texUnit = 1;
		geom->setTexCoordArray( texUnit, pos.get() );
		geom->setNormalBinding( osg::Geometry::BIND_PER_VERTEX );
		geom->setColorBinding( osg::Geometry::BIND_PER_VERTEX );
		geom->addPrimitiveSet( new osg::DrawArrays( osg::PrimitiveSet::QUADS, 0, 4 * 2 * numBonds ) );
		geode->addDrawable( geom.get() );

		for( int b = 0; b != numBonds; ++b )
		{
			OBBond* bond = mol->GetBond( b );
			if( !bond ) continue;
			OBAtom* begin = bond->GetBeginAtom();
			OBAtom* end = bond->GetEndAtom();
			const bool a0 = liquorice && bonded.find( begin ) == bonded.end();
			const bool a1 = liquorice && bonded.find( end ) == bonded.end();
			osg::Vec3d p0( begin->GetX(), begin->GetY(), begin->GetZ() );
			osg::Vec3d p1( end->GetX(), end->GetY(), end->GetZ() );
			const osg::Vec3d half = ( p0 + p1 ) * 0.5;
			const osg::Vec3d bh = ( p0 + half ) * 0.5;
			const osg::Vec3d th = ( half + p1 ) * 0.5;
			const osg::Vec3d bdir = half - p0;
			const osg::Vec3d tdir = p1 - half;

			const double l = bdir.length();
			const double r = 0.2;
			osg::Matrixd R;
			osg::Vec3d d = bdir;
			d.normalize();
			R.makeRotate( osg::Vec3d( 0., 1., 0. ), d );
			osg::BoundingBox bbox;
			bbox.expandBy( osg::Vec3d( -r, -.5*l, -r ) * R );
			bbox.expandBy( osg::Vec3d(  r, -.5*l, -r ) * R );
			bbox.expandBy( osg::Vec3d(  r,  .5*l, -r ) * R );
			bbox.expandBy( osg::Vec3d( -r,  .5*l, -r ) * R );
			bbox.expandBy( osg::Vec3d( -r, -.5*l,  r ) * R );
			bbox.expandBy( osg::Vec3d(  r, -.5*l,  r ) * R );
			bbox.expandBy( osg::Vec3d(  r,  .5*l, r ) * R );
			bbox.expandBy( osg::Vec3d( -r,  .5*l, r ) * R );
			const double br = bbox.radius();


			quads->push_back( bh + osg::Vec3d( -br , -br, br ) );
			quads->push_back( bh + osg::Vec3d(  br , -br, br ) );
			quads->push_back( bh + osg::Vec3d(  br ,  br, br ) );
			quads->push_back( bh + osg::Vec3d( -br ,  br, br ) );
			dir->push_back( bdir );
			dir->push_back( bdir );
			dir->push_back( bdir );
			dir->push_back( bdir );
			colors->push_back( gAtomColors[ begin->GetAtomicNum() ] );
			colors->push_back( gAtomColors[ begin->GetAtomicNum() ] );
			colors->push_back( gAtomColors[ begin->GetAtomicNum() ] );
			colors->push_back( gAtomColors[ begin->GetAtomicNum() ] );
			pos->push_back( bh );
			pos->push_back( bh );
			pos->push_back( bh );
			pos->push_back( bh );

			quads->push_back( th + osg::Vec3d( -br , -br, br ) );
			quads->push_back( th + osg::Vec3d(  br , -br, br ) );
			quads->push_back( th + osg::Vec3d(  br ,  br, br ) );
			quads->push_back( th + osg::Vec3d( -br ,  br, br ) );
			dir->push_back( tdir );
			dir->push_back( tdir );
			dir->push_back( tdir );
			dir->push_back( tdir );
			colors->push_back( gAtomColors[ end->GetAtomicNum() ] );
			colors->push_back( gAtomColors[ end->GetAtomicNum() ] );
			colors->push_back( gAtomColors[ end->GetAtomicNum() ] );
			colors->push_back( gAtomColors[ end->GetAtomicNum() ] );
			pos->push_back( th );
			pos->push_back( th );
			pos->push_back( th );
			pos->push_back( th );

			if( liquorice )
			{
				bonded.insert( begin );
				bonded.insert( end   );
			}
		}
	}
	return group.release();
}

//------------------------------------------------------------------------------
// RAY TRACED BOX
void AddBox( const osg::Vec3d& p, double radius, osg::Geometry* geom, int n )
{
	osg::ref_ptr< osg::Vec4Array > colors = static_cast< osg::Vec4Array* >( geom->getColorArray()  );
	osg::ref_ptr< osg::Vec3Array > quads  = static_cast< osg::Vec3Array* >( geom->getVertexArray() );
	osg::ref_ptr< osg::Vec3Array > pos	  = static_cast< osg::Vec3Array* >( geom->getTexCoordArray( 1 ) );

	// back
	quads->push_back( p + osg::Vec3d(  radius, -radius, -radius ) );
	quads->push_back( p + osg::Vec3d(  -radius, -radius, -radius ) );
	quads->push_back( p + osg::Vec3d( -radius, radius, -radius ) );
	quads->push_back( p + osg::Vec3d(  radius, radius, -radius ) );
	pos->push_back( p );
	pos->push_back( p );
	pos->push_back( p );
	pos->push_back( p );
	colors->push_back( gAtomColors[ n ] );
	colors->push_back( gAtomColors[ n ] );
	colors->push_back( gAtomColors[ n ] );
	colors->push_back( gAtomColors[ n ] );
	// front
	quads->push_back( p + osg::Vec3d( -radius, -radius, radius ) );
	quads->push_back( p + osg::Vec3d(  radius, -radius, radius ) );
	quads->push_back( p + osg::Vec3d(  radius, radius, radius ) );
	quads->push_back( p + osg::Vec3d( -radius, radius, radius ) );
	pos->push_back( p );
	pos->push_back( p );
	pos->push_back( p );
	pos->push_back( p );
	colors->push_back( gAtomColors[ n ] );
	colors->push_back( gAtomColors[ n ] );
	colors->push_back( gAtomColors[ n ] );
	colors->push_back( gAtomColors[ n ] );
	// top
	quads->push_back( p + osg::Vec3d( -radius, radius, -radius ) );
	quads->push_back( p + osg::Vec3d( -radius, radius, radius ) );
	quads->push_back( p + osg::Vec3d(  radius, radius, radius ) );
	quads->push_back( p + osg::Vec3d(  radius, radius, -radius ) );
	pos->push_back( p );
	pos->push_back( p );
	pos->push_back( p );
	pos->push_back( p );
	colors->push_back( gAtomColors[ n ] );
	colors->push_back( gAtomColors[ n ] );
	colors->push_back( gAtomColors[ n ] );
	colors->push_back( gAtomColors[ n ] );
	// bottom
	quads->push_back( p + osg::Vec3d( -radius, -radius, -radius ) );
	quads->push_back( p + osg::Vec3d(  radius, -radius, -radius ) );
	quads->push_back( p + osg::Vec3d(  radius, -radius, radius ) );
	quads->push_back( p + osg::Vec3d( -radius, -radius, radius ) );
	pos->push_back( p );
	pos->push_back( p );
	pos->push_back( p );
	pos->push_back( p );
	colors->push_back( gAtomColors[ n ] );
	colors->push_back( gAtomColors[ n ] );
	colors->push_back( gAtomColors[ n ] );
	colors->push_back( gAtomColors[ n ] );
	// left
	quads->push_back( p + osg::Vec3d( -radius, -radius, -radius ) );
	quads->push_back( p + osg::Vec3d( -radius, -radius,  radius ) );
	quads->push_back( p + osg::Vec3d( -radius, radius,   radius ) );
	quads->push_back( p + osg::Vec3d( -radius, radius, -radius ) );
	pos->push_back( p );
	pos->push_back( p );
	pos->push_back( p );
	pos->push_back( p );
	colors->push_back( gAtomColors[ n ] );
	colors->push_back( gAtomColors[ n ] );
	colors->push_back( gAtomColors[ n ] );
	colors->push_back( gAtomColors[ n ] );
	// right
	quads->push_back( p + osg::Vec3d( radius, -radius, -radius ) );
	quads->push_back( p + osg::Vec3d( radius, radius,  -radius ) );
	quads->push_back( p + osg::Vec3d( radius, radius,   radius ) );
	quads->push_back( p + osg::Vec3d( radius, -radius, radius ) );
	pos->push_back( p );
	pos->push_back( p );
	pos->push_back( p );
	pos->push_back( p );
	colors->push_back( gAtomColors[ n ] );
	colors->push_back( gAtomColors[ n ] );
	colors->push_back( gAtomColors[ n ] );
	colors->push_back( gAtomColors[ n ] );
}
#include "boxsphere.vert.h"
#include "aosphere.frag.h"
osg::Group* CreateRayTracedMoleculeGeometry3( const std::string& fname,
						 				      const std::string& format,
											  MoleculeDisplayStyle ds,
											  float radScale )
{

	const bool bonds = ds != MOLDISPSTYLE_SPACEFILL;
	const bool atoms = ds != MOLDISPSTYLE_WIREFRAME;
	const bool liquorice = ds == MOLDISPSTYLE_LIQUORICE;
	const float liqRadius = radScale;
	OBConversion obConversion;
	if( !bonds ) obConversion.AddOption( "b", OBConversion::INOPTIONS );
    obConversion.SetInFormat( format.c_str() );
	OBMol* mol = new OBMol;
	std::ifstream in( fname.c_str() );
	if( !in ) return 0;
	bool ok = obConversion.Read( mol, &in );
	if( !ok ) return 0;

	const int numAtoms = mol->NumAtoms();
	typedef std::map< int, osg::ref_ptr< osg::Geode > > GEODES;
	GEODES geodes;
	const MolekelElement* et = GetElementTable();
	static osg::ref_ptr< osg::Program > aprogram;
	if( aprogram == 0 )
	{
		aprogram = new osg::Program;
		aprogram->setName( "sphere3" );
		aprogram->addShader( new osg::Shader( osg::Shader::FRAGMENT, AOSPHEREFRAG ) );
		aprogram->addShader( new osg::Shader( osg::Shader::VERTEX, BOXSPHEREVERT ) );
	}

	osg::ref_ptr< osg::Group > group = new osg::Group;

	for( int a = 0; a != numAtoms && atoms; ++a )
    {
		OBAtom* atom = mol->GetAtom( a + 1 );
		osg::Vec3d p( atom->GetX(), atom->GetY(), atom->GetZ() );
		float radius = liquorice ? liqRadius : radScale * et[ atom->GetAtomicNum() ].vdwRadius;
		if( geodes.find( atom->GetAtomicNum() ) == geodes.end() )
		{
			geodes[ atom->GetAtomicNum() ] = new osg::Geode;
			osg::StateSet* set = geodes[ atom->GetAtomicNum() ]->getOrCreateStateSet();
			set->setAttributeAndModes( aprogram.get(), osg::StateAttribute::ON );
			set->addUniform( new osg::Uniform( "radius", radius ) );
			set->setMode(GL_LIGHTING, osg::StateAttribute::OFF);

			osg::ref_ptr< osg::Geometry >  geom   = new osg::Geometry;
			osg::ref_ptr< osg::Vec4Array > colors = new osg::Vec4Array;
			osg::ref_ptr< osg::Vec3Array > quads  = new osg::Vec3Array;
			osg::ref_ptr< osg::Vec3Array > pos	  = new osg::Vec3Array;

			geom->setVertexArray( quads.get() );
			geom->setColorArray( colors.get() );
			const int texUnit = 1;
			geom->setTexCoordArray( texUnit, pos.get() );
			geom->setColorBinding( osg::Geometry::BIND_PER_VERTEX );
			geodes[ atom->GetAtomicNum() ]->addDrawable( geom.get() );

			AddBox( p, radius, geom.get(), atom->GetAtomicNum() );

			group->addChild( geodes[ atom->GetAtomicNum() ].get() );
		}
		else
		{
			osg::Geometry* geom = geodes[ atom->GetAtomicNum() ]->getDrawable( 0 )->asGeometry();
			AddBox( p, radius, geom, atom->GetAtomicNum() );
		}
	}

	if( atoms )
	{
		for( GEODES::iterator i = geodes.begin(); i != geodes.end(); ++i )
		{
			osg::Geode* g = i->second.get();
			osg::Geometry* geom = g->getDrawable( 0 )->asGeometry();
			geom->addPrimitiveSet( new osg::DrawArrays( osg::PrimitiveSet::QUADS, 0,
														geom->getVertexArray()->getNumElements() ) );
		}
	}

	// BONDS
	if( bonds )
	{
		std::set< OBAtom* > bonded;
		const int numBonds = mol->NumBonds();
		osg::ref_ptr< osg::Geode > geode = new osg::Geode;

		static osg::ref_ptr< osg::Program > bprogram;
		if( bprogram == 0 )
		{
			bprogram = new osg::Program;
			bprogram->setName( "quadcylinder2" );
			bprogram->addShader( new osg::Shader( osg::Shader::FRAGMENT, QUADCYLINDERFRAG2 ) );
			bprogram->addShader( new osg::Shader( osg::Shader::VERTEX, QUADCYLINDERVERT2 ) );
		}
		osg::StateSet* set = geode->getOrCreateStateSet();
		set->setMode( GL_CULL_FACE, osg::StateAttribute::OFF );
		set->setAttributeAndModes( bprogram.get(), osg::StateAttribute::ON );
		set->addUniform( new osg::Uniform( "radius", liqRadius ) );
		set->setMode(GL_LIGHTING, osg::StateAttribute::OFF);
		group->addChild( geode.get() );

		osg::ref_ptr< osg::Geometry >  geom   = new osg::Geometry;
		osg::ref_ptr< osg::Vec4Array > colors = new osg::Vec4Array;
		osg::ref_ptr< osg::Vec3Array > quads  = new osg::Vec3Array;
		osg::ref_ptr< osg::Vec3Array > dir	  = new osg::Vec3Array;
		osg::ref_ptr< osg::Vec3Array > pos	  = new osg::Vec3Array;
		geom->setVertexArray( quads.get() );
		geom->setColorArray( colors.get() );
		geom->setNormalArray( dir.get() );
		const int texUnit = 1;
		geom->setTexCoordArray( texUnit, pos.get() );
		geom->setNormalBinding( osg::Geometry::BIND_PER_VERTEX );
		geom->setColorBinding( osg::Geometry::BIND_PER_VERTEX );
		geom->addPrimitiveSet( new osg::DrawArrays( osg::PrimitiveSet::QUADS, 0, 4 * 2 * numBonds ) );
		geode->addDrawable( geom.get() );

		for( int b = 0; b != numBonds; ++b )
		{
			OBBond* bond = mol->GetBond( b );
			if( !bond ) continue;
			OBAtom* begin = bond->GetBeginAtom();
			OBAtom* end = bond->GetEndAtom();
			const bool a0 = liquorice && bonded.find( begin ) == bonded.end();
			const bool a1 = liquorice && bonded.find( end ) == bonded.end();
			osg::Vec3d p0( begin->GetX(), begin->GetY(), begin->GetZ() );
			osg::Vec3d p1( end->GetX(), end->GetY(), end->GetZ() );
			const osg::Vec3d half = ( p0 + p1 ) * 0.5;
			const osg::Vec3d bh = ( p0 + half ) * 0.5;
			const osg::Vec3d th = ( half + p1 ) * 0.5;
			const osg::Vec3d bdir = half - p0;
			const osg::Vec3d tdir = p1 - half;

			const double l = bdir.length();
			const double r = 0.2;
			osg::Matrixd R;
			osg::Vec3d d = bdir;
			d.normalize();
			R.makeRotate( osg::Vec3d( 0., 1., 0. ), d );
			osg::BoundingBox bbox;
			bbox.expandBy( osg::Vec3d( -r, -.5*l, -r ) * R );
			bbox.expandBy( osg::Vec3d(  r, -.5*l, -r ) * R );
			bbox.expandBy( osg::Vec3d(  r,  .5*l, -r ) * R );
			bbox.expandBy( osg::Vec3d( -r,  .5*l, -r ) * R );
			bbox.expandBy( osg::Vec3d( -r, -.5*l,  r ) * R );
			bbox.expandBy( osg::Vec3d(  r, -.5*l,  r ) * R );
			bbox.expandBy( osg::Vec3d(  r,  .5*l, r ) * R );
			bbox.expandBy( osg::Vec3d( -r,  .5*l, r ) * R );
			const double br = bbox.radius();


			quads->push_back( bh + osg::Vec3d( -br , -br, br ) );
			quads->push_back( bh + osg::Vec3d(  br , -br, br ) );
			quads->push_back( bh + osg::Vec3d(  br ,  br, br ) );
			quads->push_back( bh + osg::Vec3d( -br ,  br, br ) );
			dir->push_back( bdir );
			dir->push_back( bdir );
			dir->push_back( bdir );
			dir->push_back( bdir );
			colors->push_back( gAtomColors[ begin->GetAtomicNum() ] );
			colors->push_back( gAtomColors[ begin->GetAtomicNum() ] );
			colors->push_back( gAtomColors[ begin->GetAtomicNum() ] );
			colors->push_back( gAtomColors[ begin->GetAtomicNum() ] );
			pos->push_back( bh );
			pos->push_back( bh );
			pos->push_back( bh );
			pos->push_back( bh );

			quads->push_back( th + osg::Vec3d( -br , -br, br ) );
			quads->push_back( th + osg::Vec3d(  br , -br, br ) );
			quads->push_back( th + osg::Vec3d(  br ,  br, br ) );
			quads->push_back( th + osg::Vec3d( -br ,  br, br ) );
			dir->push_back( tdir );
			dir->push_back( tdir );
			dir->push_back( tdir );
			dir->push_back( tdir );
			colors->push_back( gAtomColors[ end->GetAtomicNum() ] );
			colors->push_back( gAtomColors[ end->GetAtomicNum() ] );
			colors->push_back( gAtomColors[ end->GetAtomicNum() ] );
			colors->push_back( gAtomColors[ end->GetAtomicNum() ] );
			pos->push_back( th );
			pos->push_back( th );
			pos->push_back( th );
			pos->push_back( th );

			if( liquorice )
			{
				bonded.insert( begin );
				bonded.insert( end   );
			}
		}
	}
	return group.release();
}

//------------------------------------------------------------------------------
osg::Group* CreatePointRayTracedMoleculeGeometry( const std::string& fname,
						 				          const std::string& format,
												  MoleculeDisplayStyle ds,
												  float radScale )
{
	const bool bonds = ds != MOLDISPSTYLE_SPACEFILL;
	const bool atoms = ds != MOLDISPSTYLE_WIREFRAME && ds != MOLDISPSTYLE_LIQUORICE;
	const bool liquorice = ds == MOLDISPSTYLE_LIQUORICE;
	const float liqRadius = 0.2f;
	OBConversion obConversion;
	if( !bonds ) obConversion.AddOption( "b", OBConversion::INOPTIONS );
    obConversion.SetInFormat( format.c_str() );
	OBMol* mol = new OBMol;
	std::ifstream in( fname.c_str() );
	if( !in ) return 0;
	bool ok = obConversion.Read( mol, &in );
	if( !ok ) return 0;

	const int numAtoms = mol->NumAtoms();
	typedef std::map< int, osg::ref_ptr< osg::Geode > > GEODES;
	GEODES geodes;
	const MolekelElement* et = GetElementTable();
	osg::ref_ptr< osg::Group > group = new osg::Group;

	//ATOMS
	static osg::ref_ptr< osg::Program > aprogram;
	if( aprogram == 0 && atoms )
	{
		aprogram = new osg::Program;
		aprogram->setName( "pointsphere" );
		aprogram->addShader( new osg::Shader( osg::Shader::FRAGMENT, POINTSPHEREFRAG ) );
		aprogram->addShader( new osg::Shader( osg::Shader::VERTEX, POINTSPHEREVERT ) );
	}

	const bool perVertexColor = true; //useful for AO computation

	for( int a = 0; a != numAtoms && atoms; ++a )
    {
		OBAtom* atom = mol->GetAtom( a + 1 );
		osg::Vec3d p( atom->GetX(), atom->GetY(), atom->GetZ() );
		osg::ref_ptr< osg::Geometry > geom;
		if( geodes.find( atom->GetAtomicNum() ) == geodes.end() )
		{
			geodes[ atom->GetAtomicNum() ] = new osg::Geode;
		   	osg::ref_ptr< osg::Vec4Array > colors = new osg::Vec4Array;
			colors->push_back( gAtomColors[ atom->GetAtomicNum() ] );
			osg::StateSet* set = geodes[ atom->GetAtomicNum() ]->getOrCreateStateSet();
			set->setAttributeAndModes( aprogram.get(), osg::StateAttribute::ON );
			const float radius = bonds ? radScale * et[ atom->GetAtomicNum() ].vdwRadius : et[ atom->GetAtomicNum() ].vdwRadius;
			//osg::Point *point = new osg::Point();
			//point->setSize( 32.0f );
			//set->setAttribute(point);
			set->addUniform( new osg::Uniform( "radius", radius ) );
			set->setMode(GL_LIGHTING, osg::StateAttribute::OFF);
			set->setMode( GL_VERTEX_PROGRAM_POINT_SIZE, osg::StateAttribute::ON );
			geom = new osg::Geometry;
			geodes[ atom->GetAtomicNum() ]->addDrawable( geom.get() );
			if( gAtomColors.size() )
			{
				osg::ref_ptr< osg::Vec4Array > colors = new osg::Vec4Array;
				colors->push_back( gAtomColors[ std::min( size_t( atom->GetAtomicNum() ), gAtomColors.size() - 1 ) ] );
				geom->setColorArray( colors.get() );
				if( perVertexColor ) geom->setColorBinding( osg::Geometry::BIND_PER_VERTEX );
				else geom->setColorBinding( osg::Geometry::BIND_OVERALL );
			}
			osg::ref_ptr< osg::Vec3Array > coord = new osg::Vec3Array;
			coord->push_back( p );
			geom->setVertexArray( coord.get() );
			group->addChild( geodes[ atom->GetAtomicNum() ].get() );
		}
		else
		{
			geom = geodes[ atom->GetAtomicNum() ]->getDrawable( 0 )->asGeometry();
			osg::Vec3Array* coord = static_cast< osg::Vec3Array* >( geom->getVertexArray() );
			if( !coord ) continue;
			coord->push_back( p );
			if( perVertexColor )
			{
				osg::Vec4Array* colors = static_cast< osg::Vec4Array* >( geom->getColorArray() );
				colors->push_back( gAtomColors[ atom->GetAtomicNum() ] );
			}
		}
	}

	GEODES::iterator pi = geodes.begin();
	const GEODES::const_iterator pe = geodes.end();
	for( ; pi != pe && atoms; ++pi )
	{
		osg::Geometry*  g = pi->second->getDrawable( 0 )->asGeometry();
		if( !g ) continue;
		const osg::Vec3Array* coord = static_cast< const osg::Vec3Array* >( g->getVertexArray() );
		if( !coord ) continue;
		g->addPrimitiveSet( new osg::DrawArrays( osg::PrimitiveSet::POINTS, 0, coord->size() ) );
	}


	// BONDS
	if( bonds )
	{
		std::set< OBAtom* > bonded;
		const int numBonds = mol->NumBonds();
		for( int b = 0; b != numBonds; ++b )
		{
			OBBond* bond = mol->GetBond( b );
			if( !bond ) continue;
			OBAtom* begin = bond->GetBeginAtom();
			OBAtom* end = bond->GetEndAtom();
			const bool a0 = liquorice && bonded.find( begin ) == bonded.end();
			const bool a1 = liquorice && bonded.find( end ) == bonded.end();
			osg::Vec3d p0( begin->GetX(), begin->GetY(), begin->GetZ() );
			osg::Vec3d p1( end->GetX(), end->GetY(), end->GetZ() );
			group->addChild( CreateBondGroup8( p0, begin->GetAtomicNum(), p1, end->GetAtomicNum(), a0, a1 ) );
			if( liquorice )
			{
				bonded.insert( begin );
				bonded.insert( end   );
			}
		}
	}
	return group.release();
}

//------------------------------------------------------------------------------
#include "npointcylinder.vert.h"
#include "accnpointcylinder.frag.h"
#include "coolpointsphere.frag.h"
#include "aopointsphere.frag.h"
#include "coolaopointsphere.frag.h"
#include "aopointsphere.vert.h"
osg::Group* CreatePointRayTracedMoleculeGeometry2( const std::string& fname,
						 				           const std::string& format,
												   MoleculeDisplayStyle ds,
												   float radScale )
{
	const bool bonds = ds != MOLDISPSTYLE_SPACEFILL;
	const bool atoms = ds != MOLDISPSTYLE_WIREFRAME;// && ds != MOLDISPSTYLE_LIQUORICE;
	const bool liquorice = ds == MOLDISPSTYLE_LIQUORICE;
	const float liqRadius = radScale;
	OBConversion obConversion;
	if( !bonds ) obConversion.AddOption( "b", OBConversion::INOPTIONS );
    obConversion.SetInFormat( format.c_str() );
	OBMol* mol = new OBMol;
	std::ifstream in( fname.c_str() );
	if( !in ) return 0;
	bool ok = obConversion.Read( mol, &in );
	if( !ok ) return 0;

	const int numAtoms = mol->NumAtoms();
	typedef std::map< int, osg::ref_ptr< osg::Geode > > GEODES;
	GEODES geodes;
	const MolekelElement* et = GetElementTable();
	osg::ref_ptr< osg::Group > group = new osg::Group;

	//ATOMS
	static osg::ref_ptr< osg::Program > aprogram;
	if( aprogram == 0 && atoms )
	{
		aprogram = new osg::Program;
		aprogram->setName( "pointsphere" );
		aprogram->addShader( new osg::Shader( osg::Shader::FRAGMENT, COOLAOPOINTSPHEREFRAG ) );
		aprogram->addShader( new osg::Shader( osg::Shader::VERTEX, AOPOINTSPHEREVERT ) );
	}

	const bool perVertexColor = true; //useful for AO computation

	for( int a = 0; a != numAtoms && atoms; ++a )
    {
		OBAtom* atom = mol->GetAtom( a + 1 );
		osg::Vec3d p( atom->GetX(), atom->GetY(), atom->GetZ() );
		osg::ref_ptr< osg::Geometry > geom;
		if( geodes.find( atom->GetAtomicNum() ) == geodes.end() )
		{
			geodes[ atom->GetAtomicNum() ] = new osg::Geode;
		   	osg::ref_ptr< osg::Vec4Array > colors = new osg::Vec4Array;
			colors->push_back( gAtomColors[ atom->GetAtomicNum() ] );
			osg::StateSet* set = geodes[ atom->GetAtomicNum() ]->getOrCreateStateSet();
			set->setAttributeAndModes( aprogram.get(), osg::StateAttribute::ON );
			const float radius = liquorice ? liqRadius :
									( bonds ? radScale * et[ atom->GetAtomicNum() ].vdwRadius :
									  et[ atom->GetAtomicNum() ].vdwRadius );
			//osg::Point *point = new osg::Point();
			//point->setSize( 32.0f );
			//set->setAttribute(point);
			set->addUniform( new osg::Uniform( "radius", radius ) );
			set->setMode(GL_LIGHTING, osg::StateAttribute::OFF);
			set->setMode( GL_VERTEX_PROGRAM_POINT_SIZE, osg::StateAttribute::ON );
			geom = new osg::Geometry;
			geodes[ atom->GetAtomicNum() ]->addDrawable( geom.get() );
			if( gAtomColors.size() )
			{
				osg::ref_ptr< osg::Vec4Array > colors = new osg::Vec4Array;
				colors->push_back( gAtomColors[ std::min( size_t( atom->GetAtomicNum() ), gAtomColors.size() - 1 ) ] );
				geom->setColorArray( colors.get() );
				if( perVertexColor ) geom->setColorBinding( osg::Geometry::BIND_PER_VERTEX );
				else geom->setColorBinding( osg::Geometry::BIND_OVERALL );
			}
			osg::ref_ptr< osg::Vec3Array > coord = new osg::Vec3Array;
			coord->push_back( p );
			geom->setVertexArray( coord.get() );
			group->addChild( geodes[ atom->GetAtomicNum() ].get() );
		}
		else
		{
			geom = geodes[ atom->GetAtomicNum() ]->getDrawable( 0 )->asGeometry();
			osg::Vec3Array* coord = static_cast< osg::Vec3Array* >( geom->getVertexArray() );
			if( !coord ) continue;
			coord->push_back( p );
			if( perVertexColor )
			{
				osg::Vec4Array* colors = static_cast< osg::Vec4Array* >( geom->getColorArray() );
				colors->push_back( gAtomColors[ atom->GetAtomicNum() ] );
			}
		}
	}

	GEODES::iterator pi = geodes.begin();
	const GEODES::const_iterator pe = geodes.end();
	for( ; pi != pe && atoms; ++pi )
	{
		osg::Geometry*  g = pi->second->getDrawable( 0 )->asGeometry();
		if( !g ) continue;
		const osg::Vec3Array* coord = static_cast< const osg::Vec3Array* >( g->getVertexArray() );
		if( !coord ) continue;
		g->addPrimitiveSet( new osg::DrawArrays( osg::PrimitiveSet::POINTS, 0, coord->size() ) );
	}


	// BONDS
	if( bonds )
	{
		std::set< OBAtom* > bonded;
		const int numBonds = mol->NumBonds();
		osg::ref_ptr< osg::Geode > geode = new osg::Geode;
		osg::ref_ptr< osg::Geometry > geom = new osg::Geometry;
		geode->addDrawable( geom.get() );
		osg::ref_ptr< osg::Vec4Array > colors = new osg::Vec4Array;
		osg::ref_ptr< osg::Vec3Array > points = new osg::Vec3Array;
		osg::ref_ptr< osg::Vec3Array > dir = new osg::Vec3Array;
		geom->setVertexArray( points.get() );
		geom->setColorArray( colors.get() );
		geom->setNormalArray( dir.get() );
		geom->setNormalBinding( osg::Geometry::BIND_PER_VERTEX );
		geom->setColorBinding( osg::Geometry::BIND_PER_VERTEX );
		geom->addPrimitiveSet( new osg::DrawArrays( osg::PrimitiveSet::POINTS, 0, 2 * numBonds ) );

		static osg::ref_ptr< osg::Program > bprogram;
		if( bprogram == 0 )
		{
			bprogram = new osg::Program;
			bprogram->setName( "npointcylinder" );
			bprogram->addShader( new osg::Shader( osg::Shader::FRAGMENT, ACCNPOINTCYLINDERFRAG ) );
			bprogram->addShader( new osg::Shader( osg::Shader::VERTEX, NPOINTCYLINDERVERT ) );
		}
		osg::StateSet* set = geode->getOrCreateStateSet();
		set->setAttributeAndModes( bprogram.get(), osg::StateAttribute::ON );
		set->addUniform( new osg::Uniform( "radius", liqRadius ) );
		set->setMode(GL_LIGHTING, osg::StateAttribute::OFF);
		set->setMode( GL_VERTEX_PROGRAM_POINT_SIZE, osg::StateAttribute::ON );
		group->addChild( geode.get() );

		for( int b = 0; b != numBonds; ++b )
		{
			OBBond* bond = mol->GetBond( b );
			if( !bond ) continue;
			OBAtom* begin = bond->GetBeginAtom();
			OBAtom* end = bond->GetEndAtom();
			const bool a0 = liquorice && bonded.find( begin ) == bonded.end();
			const bool a1 = liquorice && bonded.find( end ) == bonded.end();
			osg::Vec3d p0( begin->GetX(), begin->GetY(), begin->GetZ() );
			osg::Vec3d p1( end->GetX(), end->GetY(), end->GetZ() );
			const osg::Vec3d half = ( p0 + p1 ) * 0.5;
			const osg::Vec3d bh = ( p0 + half ) * 0.5;
			const osg::Vec3d th = ( half + p1 ) * 0.5;
			const osg::Vec3d bdir = half - p0;
			const osg::Vec3d tdir = p1 - half;
			points->push_back( bh );
			dir->push_back( bdir );
			colors->push_back( gAtomColors[ begin->GetAtomicNum() ] );
			points->push_back( th );
			dir->push_back( tdir );
			colors->push_back( gAtomColors[ end->GetAtomicNum() ] );
			if( liquorice )
			{
				bonded.insert( begin );
				bonded.insert( end   );
			}
		}
	}
	return group.release();
}

//------------------------------------------------------------------------------
// RAY TRACED POINTS
#include "ssaopointsphere.frag.h"
#include "ssaopointsphere_trace_frag.vert.h"
#include "ssaopointsphere_trace_frag.frag.h"
#include "ssaopointsphere_trace_vert.vert.h"
#include "ssaopointsphere_trace_vert.frag.h"
#include "ssaopointsphere.vert.h"
#include "ssaoaccnpointcylinder.frag.h"
#include "ssaoaccnpointcylinder_trace.frag.h"
#include "ssaonpointcylinder.vert.h"
osg::Group* CreateAnimatedPointRayTracedMoleculeGeometry2( OBMol* mol,
												           MoleculeDisplayStyle ds,
												           float radScale,
														   bool ssao,
                                                           SSAOMode ssaoMode )
{
	if( !mol ) return 0;
	const bool bonds = ds != MOLDISPSTYLE_SPACEFILL;
	const bool atoms = ds != MOLDISPSTYLE_WIREFRAME;// && ds != MOLDISPSTYLE_LIQUORICE;
	const bool liquorice = ds == MOLDISPSTYLE_LIQUORICE;
	const float liqRadius = radScale;

	const int numAtoms = mol->NumAtoms();
	typedef std::map< int, osg::ref_ptr< osg::Geode > > GEODES;
	GEODES geodes;
	const MolekelElement* et = GetElementTable();
	osg::ref_ptr< osg::Group > group = new osg::Group;

	//ATOMS
	static osg::ref_ptr< osg::Program > aprogram;
	if( aprogram == 0 && atoms )
	{
		aprogram = new osg::Program;
		if( ssao )
		{
			aprogram->setName( "SSAO pointsphere" );
            switch( ssaoMode )
            {
            case SSAO_TRACE_FRAG: aprogram->addShader( new osg::Shader( osg::Shader::FRAGMENT, SSAOPOINTSPHEREFRAG_TRACE_FRAG ) );
                                  aprogram->addShader( new osg::Shader( osg::Shader::VERTEX, SSAOPOINTSPHEREVERT_TRACE_FRAG ) );
                                  break;
            case SSAO_TRACE_VERT: aprogram->addShader( new osg::Shader( osg::Shader::FRAGMENT, SSAOPOINTSPHEREFRAG_TRACE_VERT ) );
                                  aprogram->addShader( new osg::Shader( osg::Shader::VERTEX, SSAOPOINTSPHEREVERT_TRACE_VERT ) );
                                  break;
            default: aprogram->addShader( new osg::Shader( osg::Shader::FRAGMENT, SSAOPOINTSPHEREFRAG ) );
                     aprogram->addShader( new osg::Shader( osg::Shader::VERTEX, SSAOPOINTSPHEREVERT ) );
                     break;
            }

		}
		else
		{
			aprogram->setName( "pointsphere" );
			aprogram->addShader( new osg::Shader( osg::Shader::FRAGMENT, COOLAOPOINTSPHEREFRAG ) );
			aprogram->addShader( new osg::Shader( osg::Shader::VERTEX, AOPOINTSPHEREVERT ) );
		}
	}

	const bool perVertexColor = true; //useful for AO computation

	for( int a = 0; a != numAtoms && atoms; ++a )
    {
		OBAtom* atom = mol->GetAtom( a + 1 );
		osg::Vec3d p( atom->GetX(), atom->GetY(), atom->GetZ() );
		osg::ref_ptr< osg::Geometry > geom;
		if( geodes.find( atom->GetAtomicNum() ) == geodes.end() )
		{
			geodes[ atom->GetAtomicNum() ] = new osg::Geode;
		   	osg::ref_ptr< osg::Vec4Array > colors = new osg::Vec4Array;
			colors->push_back( gAtomColors[ atom->GetAtomicNum() ] );
			osg::StateSet* set = geodes[ atom->GetAtomicNum() ]->getOrCreateStateSet();
			set->setAttributeAndModes( aprogram.get(), osg::StateAttribute::ON );
			const float radius = liquorice ? liqRadius :
									( bonds ? radScale * et[ atom->GetAtomicNum() ].vdwRadius :
									  et[ atom->GetAtomicNum() ].vdwRadius );
			//osg::Point *point = new osg::Point();
			//point->setSize( 32.0f );
			//set->setAttribute(point);
			set->addUniform( new osg::Uniform( "radius", radius ) );
			set->setMode(GL_LIGHTING, osg::StateAttribute::OFF);
			set->setMode( GL_VERTEX_PROGRAM_POINT_SIZE, osg::StateAttribute::ON );
			geom = new osg::Geometry;

			// super simple approximation of bounding box: simply add 2.0 to radius;
			// this is required because atom are added as points i.e. spheres with zero radius
			struct BBCBack : public osg::Drawable::ComputeBoundingBoxCallback
			{
				osg::BoundingBox computeBound( const osg::Drawable& d ) const
				{
					const float dr = 2.0;
					const osg::BoundingBox& bb = d.computeBoundingBox();
					return osg::BoundingBox( bb.xMin() - dr, bb.yMin() - dr, bb.zMin() - dr,
											 bb.xMax() + dr, bb.yMax() + dr, bb.zMax() + dr );
				}
			};
			geom->setComputeBoundingBoxCallback( new BBCBack );

			geodes[ atom->GetAtomicNum() ]->addDrawable( geom.get() );
			if( gAtomColors.size() )
			{
				osg::ref_ptr< osg::Vec4Array > colors = new osg::Vec4Array;
				colors->push_back( gAtomColors[ std::min( size_t( atom->GetAtomicNum() ), gAtomColors.size() - 1 ) ] );
				geom->setColorArray( colors.get() );
				if( perVertexColor ) geom->setColorBinding( osg::Geometry::BIND_PER_VERTEX );
				else geom->setColorBinding( osg::Geometry::BIND_OVERALL );
			}
			osg::ref_ptr< osg::Vec3Array > coord = new osg::Vec3Array;
			coord->push_back( p );
			geom->setVertexArray( coord.get() );
			group->addChild( geodes[ atom->GetAtomicNum() ].get() );
		}
		else
		{
			geom = geodes[ atom->GetAtomicNum() ]->getDrawable( 0 )->asGeometry();
			osg::Vec3Array* coord = static_cast< osg::Vec3Array* >( geom->getVertexArray() );
			if( !coord ) continue;
			coord->push_back( p );
			if( perVertexColor )
			{
				osg::Vec4Array* colors = static_cast< osg::Vec4Array* >( geom->getColorArray() );
				colors->push_back( gAtomColors[ atom->GetAtomicNum() ] );
			}
		}
	}

	GEODES::iterator pi = geodes.begin();
	const GEODES::const_iterator pe = geodes.end();
	for( ; pi != pe && atoms; ++pi )
	{
		osg::Geometry*  g = pi->second->getDrawable( 0 )->asGeometry();
		if( !g ) continue;
		const osg::Vec3Array* coord = static_cast< const osg::Vec3Array* >( g->getVertexArray() );
		if( !coord ) continue;
		g->addPrimitiveSet( new osg::DrawArrays( osg::PrimitiveSet::POINTS, 0, coord->size() ) );
	}


	// BONDS
	if( bonds )
	{
		std::set< OBAtom* > bonded;
		const int numBonds = mol->NumBonds();
		osg::ref_ptr< osg::Geode > geode = new osg::Geode;
		osg::ref_ptr< osg::Geometry > geom = new osg::Geometry;
		geode->addDrawable( geom.get() );
		osg::ref_ptr< osg::Vec4Array > colors = new osg::Vec4Array;
		osg::ref_ptr< osg::Vec3Array > points = new osg::Vec3Array;
		osg::ref_ptr< osg::Vec3Array > dir = new osg::Vec3Array;
		geom->setVertexArray( points.get() );
		geom->setColorArray( colors.get() );
		geom->setNormalArray( dir.get() );
		geom->setNormalBinding( osg::Geometry::BIND_PER_VERTEX );
		geom->setColorBinding( osg::Geometry::BIND_PER_VERTEX );
		geom->addPrimitiveSet( new osg::DrawArrays( osg::PrimitiveSet::POINTS, 0, 2 * numBonds ) );

		static osg::ref_ptr< osg::Program > bprogram;
		if( bprogram == 0 )
		{
			bprogram = new osg::Program;
			if( ssao )
			{
				bprogram->setName( "SSAO npointcylinder" );
                switch( ssaoMode )
                {
                case SSAO_TRACE_FRAG: bprogram->addShader( new osg::Shader( osg::Shader::FRAGMENT, SSAOACCNPOINTCYLINDERFRAG_TRACE ) );
                                      break;
                default: bprogram->addShader( new osg::Shader( osg::Shader::FRAGMENT, SSAOACCNPOINTCYLINDERFRAG ) );
                         break;
                }
				bprogram->addShader( new osg::Shader( osg::Shader::VERTEX, SSAONPOINTCYLINDERVERT ) );
			}
			else
			{
				bprogram->setName( "npointcylinder" );
				bprogram->addShader( new osg::Shader( osg::Shader::FRAGMENT, ACCNPOINTCYLINDERFRAG ) );
				bprogram->addShader( new osg::Shader( osg::Shader::VERTEX, NPOINTCYLINDERVERT ) );
			}
		}
		osg::StateSet* set = geode->getOrCreateStateSet();
		set->setAttributeAndModes( bprogram.get(), osg::StateAttribute::ON );
		set->addUniform( new osg::Uniform( "radius", liqRadius ) );
		set->setMode(GL_LIGHTING, osg::StateAttribute::OFF);
		set->setMode( GL_VERTEX_PROGRAM_POINT_SIZE, osg::StateAttribute::ON );
		group->addChild( geode.get() );

		for( int b = 0; b != numBonds; ++b )
		{
			OBBond* bond = mol->GetBond( b );
			if( !bond ) continue;
			OBAtom* begin = bond->GetBeginAtom();
			OBAtom* end = bond->GetEndAtom();
			const bool a0 = liquorice && bonded.find( begin ) == bonded.end();
			const bool a1 = liquorice && bonded.find( end ) == bonded.end();
			osg::Vec3d p0( begin->GetX(), begin->GetY(), begin->GetZ() );
			osg::Vec3d p1( end->GetX(), end->GetY(), end->GetZ() );
			const osg::Vec3d half = ( p0 + p1 ) * 0.5;
			const osg::Vec3d bh = ( p0 + half ) * 0.5;
			const osg::Vec3d th = ( half + p1 ) * 0.5;
			const osg::Vec3d bdir = half - p0;
			const osg::Vec3d tdir = p1 - half;
			points->push_back( bh );
			dir->push_back( bdir );
			colors->push_back( gAtomColors[ begin->GetAtomicNum() ] );
			points->push_back( th );
			dir->push_back( tdir );
			colors->push_back( gAtomColors[ end->GetAtomicNum() ] );
			if( liquorice )
			{
				bonded.insert( begin );
				bonded.insert( end   );
			}
		}
	}
	return group.release();
}


//------------------------------------------------------------------------------
osg::Group* CreateImpostorRayTracedMoleculeGeometry( const std::string& fname,
						 				             const std::string& format )
{

	OBConversion obConversion;
	obConversion.AddOption( "b", OBConversion::INOPTIONS );
    obConversion.SetInFormat( format.c_str() );
	OBMol* mol = new OBMol;
	std::ifstream in( fname.c_str() );
	if( !in ) return 0;
	bool ok = obConversion.Read( mol, &in );
	if( !ok ) return 0;

	const int numAtoms = mol->NumAtoms();
	typedef std::map< int, osg::ref_ptr< osg::Billboard > > BILLBOARDS;
	BILLBOARDS billboards;
	const MolekelElement* et = GetElementTable();

	static osg::ref_ptr< osg::Program > aprogram;
	if( aprogram == 0 )
	{
		aprogram = new osg::Program;
		aprogram->setName( "sphere" );
		aprogram->addShader( osg::Shader::readShaderFile( osg::Shader::FRAGMENT, "sphere3.frag" ) );
		aprogram->addShader( osg::Shader::readShaderFile( osg::Shader::VERTEX, "sphere3.vert" ) );
	}

	osg::ref_ptr< osg::Group > group = new osg::Group;



	struct VPCback : osg::Uniform::Callback
	{
		void operator()( osg::Uniform* u, osg::NodeVisitor* nv )
		{
			if( !u || !nv ) return;
			osgGA::EventVisitor* ev = dynamic_cast< osgGA::EventVisitor* >( nv );
			if( ev )
			{
				std::cout << "osgGA::EventVisitor" << std::endl;
				osg::ref_ptr< osgGA::Event > e = ev->getEvents().back();
                osgGA::GUIEventAdapter* ev = static_cast< osgGA::GUIEventAdapter* >(e.get());
                osg::Vec2 vp( ev->getWindowWidth(), ev->getWindowHeight() );
				u->set( vp );
				return;
			}
			osgUtil::CullVisitor* cv = dynamic_cast< osgUtil::CullVisitor* >( nv );
			if( cv )
			{
				std::cout << "osgGA::CullVisitor" << std::endl;
				osg::Vec2 vp( cv->getState()->getCurrentViewport()->width(), cv->getState()->getCurrentViewport()->height() );
				u->set( vp );
				return;
			}
		}
	};


	for( int a = 0; a != numAtoms; ++a )
    {
		OBAtom* atom = mol->GetAtom( a + 1 );
		osg::Vec3d p( atom->GetX(), atom->GetY(), atom->GetZ() );
		osg::ref_ptr< osg::Geometry > geom;
		if( billboards.find( atom->GetAtomicNum() ) == billboards.end() )
		{
			billboards[ atom->GetAtomicNum() ] = new osg::Billboard;
		    billboards[ atom->GetAtomicNum() ]->setMode( osg::Billboard::POINT_ROT_EYE );
			osg::ref_ptr< osg::Vec4Array > colors = new osg::Vec4Array;
			colors->push_back( gAtomColors[ atom->GetAtomicNum() ] );

			osg::StateSet* set = billboards[ atom->GetAtomicNum() ]->getOrCreateStateSet();
			set->setMode( GL_BLEND, osg::StateAttribute::ON );
			set->setAttributeAndModes( aprogram.get(), osg::StateAttribute::ON );
			float radius = et[ atom->GetAtomicNum() ].vdwRadius;
			// POINT
			//osg::ref_ptr< osg::Uniform > vp = new osg::Uniform( "vport",  osg::Vec2( 512, 512 ) );
			//vp->setUpdateCallback( new VPCback );
			//set->addUniform( vp.get() );
			//osg::Point *point = new osg::Point();
			//point->setSize( 63.0f );
			//set->setAttribute(point);
			set->addUniform( new osg::Uniform( "radius", radius ) );
			set->setMode(GL_LIGHTING, osg::StateAttribute::OFF);

			osg::Geometry* g = gAtoms[ atom->GetAtomicNum() ]->getDrawable( 0 )->asGeometry();
			geom = new osg::Geometry;
			geom->setColorArray( g->getColorArray() );
			geom->setColorBinding( osg::Geometry::BIND_OVERALL );
			osg::Vec3Array* coord = new osg::Vec3Array;
			// put quad in front: billboards are automatically rotated such that -y is parallel to z in eye
			// coordinates
			// With point
			//coord->push_back( osg::Vec3d( 0, 0, 0 ) );
			// With billboard
			coord->push_back( osg::Vec3d( -radius, -radius, -radius ) );
			coord->push_back( osg::Vec3d(  radius, -radius, -radius ) );
			coord->push_back( osg::Vec3d(  radius, -radius,  radius ) );
			coord->push_back( osg::Vec3d( -radius, -radius,  radius ) );
			geom->setVertexArray( coord );
			geom->addPrimitiveSet( new osg::DrawArrays( osg::PrimitiveSet::QUADS, 0, 4 ) );
			group->addChild( billboards[ atom->GetAtomicNum() ].get() );
		}
		else geom = billboards[ atom->GetAtomicNum() ]->getDrawable( 0 )->asGeometry();
		billboards[ atom->GetAtomicNum() ]->addDrawable( geom.get(), p );
	}


	osg::ref_ptr< osgSim::Impostor > impostor = new osgSim::Impostor();
    impostor->setImpostorThresholdToBound( static_cast< float >( 1.f ) );
    impostor->addChild( group.get() );
    //impostor->setRange(0, 0.0f, 1e7f);
    impostor->setCenter( group->getBound().center());
    return impostor.release();

	return group.release();
}

//------------------------------------------------------------------------------
osg::Vec3d CRomSpline( const osg::Vec3d& P0, const osg::Vec3d& P1,
					   const osg::Vec3d& P2, const osg::Vec3d& P3,
					   double t )
{
	return ( ( P1 * 2 ) + ( -P0 + P2 ) * t + ( P0 * 2 - P1 * 5 + P2 * 4 - P3 ) * t * t +
		      ( -P0 + P1 * 3 - P2 * 3 + P3 ) * t * t * t ) * 0.5;
}

//------------------------------------------------------------------------------
osg::Vec3Array* CRomInterpolation( const osg::Vec3Array* in, int subdivisions )
{
	assert( in );
	assert( in->size() >= 4 );
	assert( subdivisions > in->size() );
	osg::ref_ptr< osg::Vec3Array > out = new osg::Vec3Array;
	const double dt = ( in->size() - 1 ) / double( subdivisions );
	const osg::Vec3Array& p = *in;
	out->push_back( in->front() );
	for( int step = 0; step != subdivisions; ++step )
	{
		const double t = dt * step;
		const int i = int( t ); //floor t
		const double s = t - int( t );
		if( i == 0 )
		{
			out->push_back( CRomSpline( p[ 0 ], p[ 0 ], p[ 1 ], p[ 2 ], s ) );
		}
		else
		{
			out->push_back( CRomSpline( p[ i - 1 ], p[ i ], p[ std::min( i + 1, int( in->size() - 1 ) ) ], p[ std::min( i + 2, int( in->size() - 1 ) ) ], s ) );
		}
	}
	out->push_back( in->back() );
	return out.release();
}

//------------------------------------------------------------------------------
inline double B( double t )
{
	static const double c = 1.0 / 6.0;
	if( t >= -3 && t < -2 ) return c * ( ( t + 3 ) * ( t + 3 ) * ( t + 3 ) );
	if( t >= -2 && t < -1 ) return c * ( -3 * t * t * t - 15 * t * t - 21 * t - 5 );
	if( t >= -1 && t < 0  ) return c * ( 3 * t * t * t + 3 * t * t - 3 * t + 1 );
	if( t >= 0  && t < 1  ) return c * ( ( 1 - t ) * ( 1 - t  ) * ( 1 - t ) );
	return 0.0;
}

//------------------------------------------------------------------------------
inline osg::Vec3d UBSpline( const osg::Vec3d& P0, const osg::Vec3d& P1,
							const osg::Vec3d& P2, const osg::Vec3d& P3,
							double t )
{
	return P0 * B( t ) + P1 * B( t - 1 ) + P2 * B( t - 2 ) + P3 * B( t - 3 );
}

//------------------------------------------------------------------------------

osg::Vec3Array* BSplineApproximation( const osg::Vec3Array* in, int subdivisions )
{
	assert( in );
	assert( in->size() >= 4 );
	assert( subdivisions > in->size() );
	BSpline b( in );
	return b.BuildSpline( 0.1 );
}

//------------------------------------------------------------------------------
osg::Group* CreateMoleculeBackBoneGeometry( const std::string& fname, const std::string& format )
{


	OBConversion obConversion;
	obConversion.AddOption( "b", OBConversion::INOPTIONS );
    obConversion.SetInFormat( format.c_str() );
	OBMol* mol = new OBMol;
	std::ifstream in( fname.c_str() );
	if( !in ) return 0;
	bool ok = obConversion.Read( mol, &in );
	if( !ok ) return 0;

	typedef std::map< int, std::vector< OBResidue* > > Chains;
	Chains chains;
	for( OBResidueIterator ri = mol->BeginResidues(); ri != mol->EndResidues(); ++ri )
	{
		OBResidue* r = *ri;
		if( !r || r->GetNumAtoms() < 3 ) continue;
		chains[ r->GetChainNum() ].push_back( r );
	}
	if( chains.empty() )
	{
		std::cerr << "EMPTY CHAINS" << std::endl;
		return 0;
	}

	osg::ref_ptr< osg::Group > group = new osg::Group;
	//Chain colors
	std::vector< osg::Vec4 > colors( 10 );
	colors[ 0 ] = osg::Vec4( 1, 0, 0, 1 );
	colors[ 1 ] = osg::Vec4( 0, 1, 0, 1 );
	colors[ 2 ] = osg::Vec4( 0, 0, 1, 1 );
	colors[ 3 ] = osg::Vec4( 1, 1, 1, 1 );
	colors[ 4 ] = osg::Vec4( 1, 0, 1, 1 );
	colors[ 5 ] = osg::Vec4( 0, 1, 1, 1 );
	colors[ 6 ] = osg::Vec4( 1, 1, 0, 1 );
	colors[ 7 ] = osg::Vec4( 0.5, 1, 0.5, 1 );
	colors[ 8 ] = osg::Vec4( 1, 0.5, 1, 1 );
	colors[ 9 ] = osg::Vec4( 0.4, 0.5, 0.6, 1 );
	int col = 10;
	for( Chains::iterator i = chains.begin(); i != chains.end(); ++i )
	{
		std::cout << "CHAIN" << std::endl;
		std::vector< OBResidue* > r = i->second;
		osg::ref_ptr< osg::Geode > g = new osg::Geode;
		osg::ref_ptr< osg::Geometry > geom = new osg::Geometry;
		osg::ref_ptr< osg::Vec3Array > parray = new osg::Vec3Array;
		
		std::vector< OBResidue* >::iterator ri;
        for(  ri = r.begin() ; ri != r.end() ; ++ri )
		{

			typedef std::vector< OBAtom* > Atoms;
			Atoms atoms = ( *ri )->GetAtoms();

			//std::cout << "Atoms: " << atoms.size() << " Name: " << res->GetName() << " Num: " << res->GetNum()
			//    << " Chain: " << res->GetChain() << " Chain num: " << res->GetChainNum() << " Idx: " << res->GetIdx()
			//	<< std::endl;

			for( Atoms::iterator ai = atoms.begin(); ai != atoms.end(); ++ai )
			{
				OBAtom* atom = *ai;
				if( atom && /*res->GetAtomID*/ atom->GetAtomicNum() == 6 && !( *ri )->IsHetAtom( atom ) ) // could use "CA" atom instead
				{
					//std::cout << "CA " << atom->GetX() << ' ' << atom->GetY() << ' ' << atom->GetZ() << std::endl;
					parray->push_back( osg::Vec3d( atom->GetX(), atom->GetY(), atom->GetZ() ) );
					break;
				}
			}
		}


		if( !parray->size() ) continue;
/*		geom->setVertexArray( parray.get() );
		g->addDrawable( geom.get() );
		osg::ref_ptr< osg::DrawArrays > prs = new osg::DrawArrays( osg::PrimitiveSet::LINE_STRIP, 0, parray->size() );
		geom->addPrimitiveSet( prs.get() );
		osg::ref_ptr< osg::StateSet > ss = geom->getOrCreateStateSet();
		osg::Material* m = new osg::Material;
		m->setDiffuse( osg::Material::FRONT, colors[ col   % 10 ] );
		m->setAmbient( osg::Material::FRONT, colors[ col++ % 10 ] );
		ss->setAttribute( m, osg::StateAttribute::ON );
*/
		// ribbon cross section
		osg::ref_ptr< osg::Vec3Array > section = new osg::Vec3Array;
		section->push_back( osg::Vec3d( -.1, -1, 0 ) );
		section->push_back( osg::Vec3d( -.1, 1, 0 ) );
		section->push_back( osg::Vec3d( .1, 1, 0 ) );
		section->push_back( osg::Vec3d( .1, -1, 0 ) );
		section->push_back( osg::Vec3d( -.1, -1, 0 ) );

		// sweep cross section along cubic spline
		geom = CreateExtrudedSurfaceGeometry( section.get(), BSplineApproximation( parray.get(), 4 * parray->size() ) );
		g->addDrawable( geom.get() );
		osg::ref_ptr< osg::StateSet > ss = geom->getOrCreateStateSet();
		osg::Material* m = new osg::Material;
		m->setDiffuse( osg::Material::FRONT, colors[ col   % colors.size() ] );
		m->setAmbient( osg::Material::FRONT, colors[ col++ % colors.size() ] );
		ss->setAttribute( m, osg::StateAttribute::ON );
		group->addChild( g.get() );
	}

	return  group.release();
}


#ifdef MAIN_
//------------------------------------------------------------------------------
int main( int argc, char **argv )
{
    // use an ArgumentParser object to manage the program arguments.
    osg::ArgumentParser arguments(&argc,argv);
   
    //osg::DisplaySettings::instance()->setMinimumNumAlphaBits(8);
   
    // construct the viewer.
    osgViewer::Viewer viewer;

    // if user request help write it out to cout.
    if (arguments.read("-h") || arguments.read("--help"))
    {
        arguments.getApplicationUsage()->write( std::cout );
        return 1;
    }
    
	std::string filepath;
	std::string cfilepath;
	std::string format;
	if( !arguments.read( "--molecule", filepath ) ) return 1;
	if( !arguments.read( "--colors", cfilepath ) ) return 1;
	if( !arguments.read( "--format", format ) ) return 1; 
	InitAtomColors( cfilepath.c_str() );
	InitAtomGeometries();
    osg::ref_ptr<osg::Group> node = CreateMoleculeGeometry( filepath, format );
	if( !node )
	{
		std::cerr << "Error reading file " << filepath << std::endl;
		return 1;
	}

	osgUtil::Optimizer optimizer;
	optimizer.optimize( node.get() );

	// add model to viewer.
    viewer.setSceneData( node.get() );

	 // add the state manipulator
    viewer.addEventHandler( new osgGA::StateSetManipulator(viewer.getCamera()->getOrCreateStateSet()) );
    
    // add the thread model handler
    viewer.addEventHandler(new osgViewer::ThreadingHandler);

    // add the window size toggle handler
    viewer.addEventHandler(new osgViewer::WindowSizeHandler);
        
    // add the stats handler
    viewer.addEventHandler(new osgViewer::StatsHandler);

    // add the help handler
    viewer.addEventHandler(new osgViewer::HelpHandler(arguments.getApplicationUsage()));

    return viewer.run();
}
#endif

