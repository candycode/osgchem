static const char COOLPOINTSPHEREFRAG[] =
"uniform float radius;\n"
"uniform vec2 osg_Viewport;\n"
"varying vec4 color;\n"
"varying vec3 rayorigin;\n"
"vec3 raydir;\n"
"varying vec3 spherepos;\n"
"float sphereradsq;\n"
"varying float pointSize;\n"
"struct I\n"
"{\n"
"  vec3 P;\n"
"  vec3 N;\n"
"  float t;  \n"
"};\n"
"// computes intersection of ray with sphere\n"
"I ComputeRaySphereIntersection( vec3 R, vec3 D )\n"
"{\n"
"  I i;\n"
"  i.t = -1.0;\n"
"  vec3 P = R - spherepos;\n"
"  float a = D.x * D.x + D.y * D.y + D.z * D.z;\n"
"  float b = 2.0 * ( D.x * P.x + D.y * P.y + D.z * P.z );\n"
"  float c = P.x * P.x + P.y * P.y + P.z * P.z - sphereradsq;\n"
"  float delta = ( b * b - 4. * a * c );\n"
"  if( delta < 0.0 ) return i;\n"
"  float d = sqrt( delta );\n"
"  a = .5 / a;\n"
"  float t2 = ( -b + d ) * a;\n"
"  float t1 = ( -b - d ) * a;\n"
"  float t = min( t1, t2 );\n"
"  if( t < 0.0 ) return i;\n"
"  i.P = R + t * D;\n"
"  i.N = normalize( i.P - spherepos );\n"
"  i.t = t;\n"
"  return i;\n"
"}\n"
"vec3 lightDir = vec3( 0., 0., -1. );\n"
"float kd = 1.0;\n"
"float ka = 0.01;\n"
"float ks = 0.5;\n"
"float sh = 90.0;\n"
"vec3 refcolor = vec3( .1, .1, .1 );\n"
"\n"
"const float C1 = 0.429043;\n"
"const float C2 = 0.511664;\n"
"const float C3 = 0.743125;\n"
"const float C4 = 0.886227;\n"
"const float C5 = 0.247708;\n"
"// Constants for Grace Cathedral lighting\n"
"const vec3 L00  = vec3( 0.78908,  0.43710,  0.54161);\n"
"const vec3 L1m1 = vec3( 0.39499,  0.34989,  0.60488);\n"
"const vec3 L10  = vec3(-0.33974, -0.18236, -0.26940);\n"
"const vec3 L11  = vec3(-0.29213, -0.05562,  0.00944);\n"
"const vec3 L2m2 = vec3(-0.11141, -0.05090, -0.12231);\n"
"const vec3 L2m1 = vec3(-0.26240, -0.22401, -0.47479);\n"
"const vec3 L20  = vec3(-0.15570, -0.09471, -0.14733);\n"
"const vec3 L21  = vec3( 0.56014,  0.21444,  0.13915);\n"
"const vec3 L22  = vec3( 0.21205, -0.05432, -0.30374);\n"
"\n"
"vec4 ComputeColor( in vec3 n )\n"
"{\n"
"  float ScaleFactor = 0.9;\n" 
"  vec3 tnorm = faceforward( n, lightDir, n );\n"
"  vec3 DiffuseColor    = C1 * L22 * (tnorm.x * tnorm.x - tnorm.y * tnorm.y) +\n"
"                     C3 * L20 * tnorm.z * tnorm.z +\n"
"                      C4 * L00 -\n"
"                      C5 * L20 +\n"
"                      2.0 * C1 * L2m2 * tnorm.x * tnorm.y +\n"
"                      2.0 * C1 * L21  * tnorm.x * tnorm.z +\n"
"                      2.0 * C1 * L2m1 * tnorm.y * tnorm.z +\n"
"                      2.0 * C2 * L11  * tnorm.x +\n"
"                      2.0 * C2 * L1m1 * tnorm.y +\n"
"                      2.0 * C2 * L10  * tnorm.z;\n"
"\n"
"  DiffuseColor   *= vec3( color ) * ScaleFactor + refcolor;\n"
"  return vec4( DiffuseColor, color.a );\n" 
"}\n"
"void main(void)\n"
"{\n"
"  if( pointSize > 255.0 || pointSize < 2.0 ) discard;\n"
"  sphereradsq = radius * radius;\n"
"  bool perspective = gl_ProjectionMatrix[ 3 ][ 3 ] < 0.001 && gl_ProjectionMatrix[ 2 ][ 3 ] != 0.0;\n"
"  vec3 fc = vec3( gl_FragCoord );// * gl_FragCoord.w;\n"
"  fc.x = gl_FragCoord.x / osg_Viewport.x;\n"
"  fc.y = gl_FragCoord.y / osg_Viewport.y;\n"
"  fc.z = gl_FragCoord.z;\n"
"  fc *= 2.0;\n"
"  fc -= 1.0;\n"
"  vec4 p = gl_ProjectionMatrixInverse * vec4( fc, 1. );\n"
"  if( perspective )\n"
"  {\n"
"    raydir = vec3( p ) / p.w;\n"
"  }  \n"
"  else\n"
"  {\n"
"    raydir = vec3( 0, 0, -1 );\n"
"  }\n"
"  I i = ComputeRaySphereIntersection( rayorigin, raydir );\n"
"  if( i.t < 0.0 ) discard;\n"
"  gl_FragColor = ComputeColor( i.N );	  \n"
"  float z = dot( vec4( i.P, 1 ), gl_ProjectionMatrixTranspose[ 2 ] );\n"
"  float w = dot( vec4( i.P, 1 ), gl_ProjectionMatrixTranspose[ 3 ] );\n"
"  gl_FragDepth = 0.5 * ( z / w + 1.0 );\n"
"}\n";
