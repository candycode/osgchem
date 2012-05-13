#include <osg/Array>
#include <osg/ref_ptr>
#include <vector>

class BSpline
{
	int n;
	int m;
	std::vector< int > nV;
	
	osg::ref_ptr< const osg::Vec3Array > cP;
	
	osg::Vec3d deBoor( int k, int i, double u )
	{
		if( k == 0 ) return ( *cP )[ i ];
		const double pre = ( u - nV[ i + k ] ) / ( nV[ i + n + 1 ] - nV[ i + k ] );
		return deBoor( k - 1, i, u ) * ( 1. - pre ) + deBoor( k - 1, i + 1, u ) * pre;
	}
	
	void createNodeVector()
	{
		int knoten = 0;
		nV.clear();
		for( int j = 0; j < ( n + m + 1 ); ++j ) nV.push_back( 0 );
		for( int i = 0; i < ( n + m + 1 ); ++i )
		{
			if( i > n )
			{
				if( i <= m ) nV[ i ] = ++knoten;
				else nV[ i ] = knoten;
			}
			else nV[ i ] = knoten;
		}
	}
public:

	BSpline( const osg::Vec3Array* cp ) : cP( cp ) { m = cp->size(); n = 3; createNodeVector(); std::cout << "!!!!!!!!!!!!" << std::endl; }
	
	osg::Vec3Array* BuildSpline( double du )
	{
		std::cout << "!!!!!!!!!!!!" << std::endl;
		osg::ref_ptr< osg::Vec3Array > out = new osg::Vec3Array;
		const double mu = nV[ m + n ];
		osg::Vec3d temp;
		for( double u=0; u < mu; u += du )
		{		
			for( int j = 0; j < m; ++j )							// until  every cP is passed
			{ 
				if( u >= double( j ) ) temp = deBoor( n, j, u );
			}
			out->push_back( temp ); 
		}
		return out.release();
	}
};
