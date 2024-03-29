static const char SSAOACCNPOINTCYLINDERFRAG_TRACE[] =
"vec3 raydir;\n"
"varying vec3 rayorigin;\n"
"varying vec4 color;\n"
"uniform vec2 osg_Viewport;\n"
"uniform int ssao;\n"
"uniform int shade;\n"
"uniform float occlusionFactor;\n"
"varying float pointSize;\n"
"varying float L;\n"
"varying vec3 D;\n"
"varying vec3 Dn;\n"
"varying vec3 P0;\n"
"varying vec3 P1;\n"
"varying float sqradius;\n"
"uniform float radius;\n"
"#extension GL_ARB_texture_rectangle : enable\n"
"uniform sampler2DRect depthMap;\n"
"uniform float samplingStep;\n"
"uniform float halfSamples;\n"
"struct I\n"
"{\n"
"  vec3 P;\n"
"  vec3 N;\n"
"  float t;\n"
"};\n"
"// returns occlusion at pixel x, y;\n"
"float ComputeOcclusion( vec4 p )\n"
"{\n"
"	float zt = 0.0;\n"
"	// attempt to correlate the radius and step with disatance from viewer; it really\n"
"	// depends on the object topology\n"
"	float step = samplingStep; //int( max( 1.0, samplingStep * ( 1.0 - .3 * z ) ) );\n"
"	int hs = int( halfSamples ); //int( max( 1.0, float( halfSamples ) * ( 1.0 - .5 * z ) ) );\n"
"	for( int i = -hs; i != hs + 1; ++i )\n"
"	{\n"
"	  for( int j = -hs; j != hs + 1; ++j )\n"
"	  {\n"
"	    //if( i == 0 && j == 0 ) continue;\n"
"	    float zz = texture2DRect( depthMap, vec2( p.x + step * float( i ), p.y + step * float( j ) ) ).x;\n"
"	    //if( zz < MAX_Z )\n"
"	    {\n"
"			zt += zz;\n"
"		}\n"
"	  }\n"
"	}\n"
"   float P = 2.0 * halfSamples + 1.0;\n"
"   P *= P;\n"
"   P = 1. / P;\n"
"	//If z > 0 it means the fragment is behind (further away from the viewpoint) the neighboring fragments.\n"
"	//The distance between the fragment and the average of its neighbors is considered an occlusion value to be\n"
"	//subtracted from or multiplied by the pixel luminance/color.\n"  
"	return ( sqrt( clamp( occlusionFactor * ( p.z - zt * P ), 0.0, 1.0 ) ) );\n"
"}\n"
"\n"
"vec3 ComputeNormal( in vec3 P )\n"
"{\n"
"  return normalize( P - (P0 + dot( P - P0, Dn ) * Dn ) );\n" 
"}\n"
"I ComputeRayCylinderIntersection()\n"
"{\n"
"  I ip;\n"
"  ip.t = -1.0;\n"
"  float a = dot( D, D );\n"
"  float b = dot( D, raydir );\n"
"  float c = dot( raydir, raydir );\n"
"  float d = dot( D, P0 - rayorigin );\n"
"  float e = dot( raydir, P0 - rayorigin );\n"
"  float delta = a * c - b * b;\n"
"  if( delta == 0.0 ) return ip;\n"
"  float tc = ( a * e - b * d ) / delta;\n"  
"  if( tc < 0. ) return ip;\n"
"  float sc = ( b * e - c * d ) / delta;\n"
"  vec3 pc = P0 + sc * D;\n"
"  vec3 pn = rayorigin + tc * raydir;\n"
"  float l = length( pn - pc );\n"
"  if( l > radius ) return ip;\n"
"  vec3 O = cross( normalize( pn - pc ), Dn );\n"
"  float dt = sqrt( sqradius - l * l ) / dot( normalize( raydir ), O );\n"	  
"  vec3 Dt = dt * normalize( raydir );\n"
"  vec3 Pc = pn - Dt;\n"
"  vec3 Pf = pn + Dt;\n" 
"  float l1 = dot( Pc - P0, Dn );\n"
"  float l2 = dot( Pf - P0, Dn );\n"
"  float d1 = length( Pc - rayorigin );\n"
"  float d2 = length( Pf - rayorigin );\n"
"  if( d1 < d2 )\n"
"  {\n"
"    if( l1 <= L && l1 > 0. )\n"
"    {\n"
"      ip.P = Pc;\n"
"      ip.N = ComputeNormal( ip.P );\n"
"      ip.t = 0.;\n"
"    }\n"
"    else if( l2 <= L && l2 > 0. )\n"
"    {\n"
"      ip.P = Pf;\n"
"      ip.N = ComputeNormal( ip.P );\n"
"      ip.t = 0.;\n"
"    }\n"
"  }\n"
"  else\n"
"  {\n"
"    if( l2 <= L && l2 > 0. )\n"
"    {\n"
"      ip.P = Pf;\n"
"      ip.N = ComputeNormal( ip.P );\n"
"      ip.t = 0.;\n"
"    }\n"
"    else if( l1 <= L && l1 > 0. )\n"
"    {\n"
"      ip.P = Pc;\n"
"      ip.N = ComputeNormal( ip.P );\n"
"      ip.t = 0.;\n"
"    }\n"
"  }\n"
"  return ip;\n"
"}\n"
"vec3 lightDir = vec3( 0, 0, -1 );\n"
"float kd = 1.0;\n"
"float ka = 0.02;\n"
"float ks = 0.0;\n"
"float sh = 90.0;\n"
"vec3 refcolor = vec3( 1, 1, 1 );\n"
"vec4 ComputeColor( vec4 color, vec3 n )\n"
"{\n"
"  if( pointSize < 4.0 ) return color;\n"
"  vec3 N = faceforward( n, lightDir, n );\n"
"  float d = dot( N, -normalize( lightDir ) );\n"
"  float s = pow( max( 0.0, dot( vec3( 0, 0, 1 ), reflect( lightDir, N ) ) ), sh );\n"
"  return vec4(  ks * s * refcolor + kd * d * color.rgb + ka * color.rgb, color.a );\n"
"}\n"
"// entry point \n"
"void main(void)\n"
"{\n"
"  //if( pointSize > 255.0 || pointSize < 1.0 ) discard;\n"
"  vec3 fc = vec3( gl_FragCoord );// * gl_FragCoord.w;\n"
"  fc.x = gl_FragCoord.x / osg_Viewport.x;\n"
"  fc.y = gl_FragCoord.y / osg_Viewport.y;\n"
"  fc.z = gl_FragCoord.z;\n"
"  fc *= 2.0;\n"
"  fc -= 1.0;\n"
"  vec4 p = gl_ProjectionMatrixInverse * vec4( fc, 1. );\n"
"  bool perspective = gl_ProjectionMatrix[ 3 ][ 3 ] < 0.001 && gl_ProjectionMatrix[ 2 ][ 3 ] != 0.0;\n"
"  if( perspective )\n"
"  {\n"
"    raydir = vec3( p ) / p.w;\n"
"  }  \n"
"  else\n"
"  {\n"
"    raydir = vec3( 0, 0, -1 );\n"
"  }\n"
"  I i = ComputeRayCylinderIntersection();\n"
"  if( i.t < 0.0 ) discard;\n"
"  if( bool( shade ) ) gl_FragColor = ComputeColor( color, i.N );\n"
"  if( bool( ssao ) && bool( shade ) )\n"
"  {\n"
"    float l = ComputeOcclusion( gl_FragCoord );\n"
"    gl_FragColor.rgb -= l;\n"
"  }\n"
"  float z = dot( vec4( i.P, 1 ), gl_ProjectionMatrixTranspose[ 2 ] );\n"
"  float w = dot( vec4( i.P, 1 ), gl_ProjectionMatrixTranspose[ 3 ] );\n"
"  gl_FragDepth = 0.5 * ( z / w + 1.0 );\n"
"}\n";
