static const char NPOINTCYLINDERFRAG[] =
"vec3 raydir;\n"
"varying vec3 rayorigin;\n"
"varying vec4 color;\n"
"uniform vec2 osg_Viewport;\n"
"varying float pointSize;\n"
"varying float L;\n"
"varying vec3 D;\n"
"varying vec3 Dn;\n"
"varying vec3 P0;\n"
"varying vec3 P1;\n"
"varying float sqradius;\n"
"uniform float radius;\n"
"struct I\n"
"{\n"
"  vec3 P;\n"
"  vec3 N;\n"
"  float t;\n"
"};\n"
"vec3 ComputeNormal( vec3 P )\n"
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
"  float sc = ( b * e - c * d ) / delta;\n"
"  if( sc > 1.0 || sc < 0. ) return ip;\n" 
"  float tc = ( a * e - b * d ) / delta;\n"  
"  if( tc < 0. ) return ip;\n"
"  vec3 pc = P0 + sc * D;\n"
"  vec3 pn = rayorigin + tc * raydir;\n"
"  float l = length( pn - pc );\n"
"  if( l > radius ) return ip;\n"
"  float dt = sqrt( sqradius - l * l );\n"
"  ip.P = pn - dt * normalize( raydir );\n"
"  ip.N = ComputeNormal( ip.P );\n"
"  ip.t = 0.;\n"
"  return ip;\n"
"}\n"
"vec3 lightDir = vec3( 0, 0, -1 );\n"
"float kd = 1.0;\n"
"float ka = 0.01;\n"
"float ks = .5;\n"
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
"  if( pointSize > 256.0 || pointSize < 2.0 ) discard;\n"
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
"  gl_FragColor = ComputeColor( color, i.N );\n"
"  float z = dot( vec4( i.P, 1 ), gl_ProjectionMatrixTranspose[ 2 ] );\n"
"  float w = dot( vec4( i.P, 1 ), gl_ProjectionMatrixTranspose[ 3 ] );\n"
"  gl_FragDepth = 0.5 * ( z / w + 1.0 );\n"
"}\n";
