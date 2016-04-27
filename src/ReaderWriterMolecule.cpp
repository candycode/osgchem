
#if defined(_MSC_VER)
    #pragma warning( disable : 4786 )
#endif

#include <string>
#include <sstream>

#include <osg/Notify>
#include <osg/Node>
#include <osg/Group>
#include <osg/MatrixTransform>
#include <osg/Geode>

#include <osg/Geometry>
#include <osg/StateSet>
#include <osg/Material>
#include <osg/Texture2D>
#include <osg/TexGen>
#include <osg/TexMat>

#include <osgDB/Registry>
#include <osgDB/ReadFile>
#include <osgDB/FileUtils>
#include <osgDB/FileNameUtils>

#include <osgUtil/TriStripVisitor>
#include <osgUtil/SmoothingVisitor>
#include <osgUtil/Tessellator>

#include <osg/Sequence>


#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/obmolecformat.h>

#include "osgob.h"

using namespace OpenBabel;

class ReaderWriterMolecule : public osgDB::ReaderWriter
{
	mutable bool ribbons_;
	mutable bool atoms_;
	mutable std::string atomColorsFilePath_;
	mutable bool billboard_;
	mutable int slices_;
	mutable int stacks_;
	mutable float scale_;
	mutable MoleculeDisplayStyle molDisplayStyle_;
	mutable double dt_;
	mutable bool animate_;
	mutable std::string renderMethod_;
	mutable bool tex_;
	mutable bool ssao_;
        mutable SSAOMode ssaoMode_;

	void ParseOptions( const std::string& o ) const
	{
		std::istringstream is( o );
		std::string s;
		while( is )
		{
			is >> s;
			if     ( s == "-atoms"     ) is >> atoms_;
			else if( s == "-ribbons"   ) is >> ribbons_;
			else if( s == "-colorfile" ) is >> atomColorsFilePath_;
			else if( s == "-slices"    ) is >> slices_;
			else if( s == "-stacks"    ) is >> stacks_;
			else if( s == "-rep"       )
			{
				std::string r;
				is >> r;
				if( r == "bs" ) molDisplayStyle_ = MOLDISPSTYLE_BALLANDSTICK;
				else if( r == "liq" ) molDisplayStyle_ = MOLDISPSTYLE_LIQUORICE;
				else if( r == "sf"  ) molDisplayStyle_ = MOLDISPSTYLE_SPACEFILL;
			}
			else if( s == "-scale"  ) is >> scale_;
			else if( s == "-dt"     ) is >> dt_;
			else if( s == "-anim"   ) animate_ = true;
			else if( s == "-method" ) is >> renderMethod_;
			else if( s == "-bb"     ) is >> billboard_;
			else if( s == "-tex"    ) is >> tex_;
			else if( s == "-ssao"   ) 
            {
                ssao_ = true;
                std::string ssaoMode;
                is >> ssaoMode;
                if( ssaoMode == "trace_frag" ) ssaoMode_ = SSAO_TRACE_FRAG;
                else if( ssaoMode == "trace_frag_mrt" ) ssaoMode_ = SSAO_TRACE_FRAG_MRT;
                else if( ssaoMode == "trace_vert" ) ssaoMode_ = SSAO_TRACE_VERT;
                else if( ssaoMode == "trace_vert_mrt" ) ssaoMode_ = SSAO_TRACE_VERT_MRT;
            }				
		}
	}

	osg::Group* ReadAnimatedMolecule( const char* fname, const char* format, MoleculeDisplayStyle ds, float scale ) const 
	{
		const bool bonds = ( ds != MOLDISPSTYLE_SPACEFILL );
		OBConversion obConversion;
		if( !bonds ) obConversion.AddOption( "b", OBConversion::INOPTIONS );
		obConversion.SetInFormat( format );
		OBMol mol;
		std::ifstream in( fname );
		if( !in ) return 0;
		osg::ref_ptr< osg::Sequence > seq = new osg::Sequence;
		bool ok = false;
		double t = 0.0;
		double dt = 0.01;
		osg::ref_ptr< osg::Group > node;
		do
		{
			ok = obConversion.Read( &mol, &in );
			if( ok )
			{
				node = CreateAnimatedPointRayTracedMoleculeGeometry2( &mol, ds, scale, ssao_, ssaoMode_ );
				if( node != 0 ) 
				{
					seq->addChild( node.get() );
					seq->setTime( seq->getNumChildren() - 1, dt_ );
				}
			}
			if( !animate_ ) break;
			t += dt;
		} while( ok );
		if( seq->getNumChildren() == 0 ) return 0;

		if( seq->getNumChildren() == 1 )
		{
			seq->setInterval( osg::Sequence::LOOP, 0, 0 );
			seq->setDuration( 0.0f, 1 );
			seq->setMode( osg::Sequence::START );
		}
		else
		{
			// loop through all children
			seq->setInterval( osg::Sequence::LOOP, 0, -1 );	
			// real-time playback, repeat indefinitively
			seq->setDuration( 1.0f, -1 );
			seq->setMode( osg::Sequence::START );
		}
		
		return seq.release();
		//seq->setSync( true );
		
	}

public:

	ReaderWriterMolecule() : ribbons_( false ), atoms_( true ), slices_( 16 ), stacks_( 8 ),
							 scale_( 1.0f ),  molDisplayStyle_( MOLDISPSTYLE_SPACEFILL ), dt_( 0.03 ),
							 animate_( false ), renderMethod_( "triangles" ), billboard_( false ),
							 tex_( false ), ssao_( false ), ssaoMode_( SSAO_SIMPLE ) {}

    virtual const char* className() const { return "Molecule Reader"; }

    virtual bool acceptsExtension( const std::string& extension) const {
        return osgDB::equalCaseInsensitive( extension, "pdb"  ) ||
			   osgDB::equalCaseInsensitive( extension, "xyz"  ) ||
			   osgDB::equalCaseInsensitive( extension, "mol"  ) ||
			   osgDB::equalCaseInsensitive( extension, "mol2" ) ||
			   osgDB::equalCaseInsensitive( extension, "cube" ) ||
			   osgDB::equalCaseInsensitive( extension, "hin"  ) ||
			   osgDB::equalCaseInsensitive( extension, "ent"  ) ||
			   osgDB::equalCaseInsensitive( extension, "g98"  ) ||
			   osgDB::equalCaseInsensitive( extension, "g03"  ) ||
			   osgDB::equalCaseInsensitive( extension, "gam"  ) ||
			   osgDB::equalCaseInsensitive( extension, "coor" ) ||
			   osgDB::equalCaseInsensitive( extension, "ref"  );
    }

    virtual ReadResult readNode(const std::string& fileName, const osgDB::ReaderWriter::Options* options) const;

    virtual ReadResult readNode(std::istream& fin, const Options* options) const;

};


// register with Registry to instantiate the above reader/writer.
REGISTER_OSGPLUGIN(molecule, ReaderWriterMolecule)

// read file and convert to OSG.
osgDB::ReaderWriter::ReadResult ReaderWriterMolecule::readNode(const std::string& file, const osgDB::ReaderWriter::Options* options) const
{
    std::string ext = osgDB::getLowerCaseFileExtension(file);
    if (!acceptsExtension(ext)) return osgDB::ReaderWriter::ReadResult::FILE_NOT_HANDLED;

    std::string fileName = osgDB::findDataFile( file, options );
    if (fileName.empty()) return osgDB::ReaderWriter::ReadResult::FILE_NOT_FOUND;

	if( ext == "coor" || ext == "ref" ) ext = "pdb";	

	if( options ) ParseOptions( options->getOptionString() );

	osg::ref_ptr< osg::Group > atoms;
	if( atoms_ )
	{
		InitAtomColors( atomColorsFilePath_.c_str() );
		if( billboard_ ) InitAtomTextureGeometries();
		else InitAtomGeometries( slices_, stacks_ );

		if( renderMethod_ == "RayTracedPoint" )
		{
			atoms = ReadAnimatedMolecule( fileName.c_str(), ext.c_str(), molDisplayStyle_, scale_ );
		}
		else if( renderMethod_ == "RayTracedPointQ" )
		{
			atoms = CreatePointRayTracedMoleculeGeometry( fileName.c_str(), ext.c_str(), molDisplayStyle_, scale_ );
		}
		else if( renderMethod_ == "RayTracedQuad" )
		{
			atoms = CreateRayTracedMoleculeGeometry2( fileName.c_str(), ext.c_str(), molDisplayStyle_, scale_ );
		}
		else if( renderMethod_ == "RayTracedBox" )
		{
			atoms = CreateRayTracedMoleculeGeometry3( fileName.c_str(), ext.c_str(), molDisplayStyle_, scale_ );
		}
		else if( renderMethod_ == "PointSprite" )
		{
			atoms = CreatePointSpriteMoleculeGeometry( fileName, ext.c_str(), molDisplayStyle_ );
		}
		else if( renderMethod_ == "BillBoard" )
		{
			atoms = CreateBillboardMoleculeGeometry( fileName, ext.c_str(), molDisplayStyle_ ); 
		}
		else
		{
			atoms = CreateMoleculeGeometry( fileName, ext.c_str(), molDisplayStyle_, billboard_, tex_ );
		}
		if( !atoms )
		{
			return osgDB::ReaderWriter::ReadResult::FILE_NOT_HANDLED;
		}
	}
	osg::ref_ptr< osg::Group > ribbons = ribbons_ ? CreateMoleculeBackBoneGeometry( fileName, ext.c_str() ) : 0;
	if( atoms != 0 && ribbons != 0 )
	{
		osg::ref_ptr< osg::Group > group = new osg::Group;
		group->addChild( atoms.get()   );
		group->addChild( ribbons.get() );
		return group.release();
	}

	if( atoms != 0 ) return atoms.release();
	else if( ribbons != 0 ) return ribbons.release();

	return 0;
}

osgDB::ReaderWriter::ReadResult ReaderWriterMolecule::readNode(std::istream& fin, const Options* options) const
{

    return osgDB::ReaderWriter::ReadResult::FILE_NOT_HANDLED;
}


