static const char SSAONPOINTCYLINDERVERT[] =
"uniform vec2 osg_Viewport;\n"
"varying vec3 rayorigin;\n"
"varying vec4 color;\n"
"varying float pointSize;\n"
"varying float L;\n"
"varying vec3 D;\n"
"varying vec3 Dn;\n"
"varying vec3 P0;\n"
"varying vec3 P1;\n"
"varying float sqradius;\n"
"uniform float radius;\n"
"uniform int shade;\n"
"/*uniform*/ float saturation = 0.5;\n"
"void main()\n"
"{\n"
"  sqradius = radius * radius;\n"
"  if( bool( shade ) ) color = vec4( saturation * gl_Color.rgb + ( 1.0 - saturation ) * vec3( 1., 1., 1. ), gl_Color.a );\n" 
"  D = gl_NormalMatrix * gl_Normal;\n"
"  vec3 taxis = normalize( D );\n"
"  Dn = taxis;\n"
"  L = length( D );\n"
"  vec4 p = gl_ModelViewMatrix * gl_Vertex;\n"
"  vec3 ep = p.xyz / p.w;\n"
"  P0 = ep - 0.5 * L * taxis;\n"
"  P1 = ep + 0.5 * L * taxis;\n"
"  bool perspective = gl_ProjectionMatrix[ 3 ][ 3 ] < 0.001 && gl_ProjectionMatrix[ 2 ][ 3 ] != 0.0;\n"
"  gl_ClipVertex = p;\n"
"  if( perspective )\n"
"  {\n"
"    rayorigin = vec3( 0., 0., 0. );\n"
"  }\n"
"  else\n"
"  {\n"
"    rayorigin = vec3( ep.x, ep.y, 0. );\n"
"  }\n"
"  gl_Position = gl_ProjectionMatrix * p;\n"
"  // compute pixel size\n"
"  vec3 axis = cross( normalize( ep - rayorigin ), vec3( 0., 1., 0. ) );\n"
"  vec4 xp1 = gl_ProjectionMatrix * vec4( ep - 0.5 * L * axis, 1.0 );\n"
"  vec4 xp2 = gl_ProjectionMatrix * vec4( ep + 0.5 * L * axis, 1.0 );\n"
"  xp1.xy = osg_Viewport * 0.5 * ( xp1.xy / xp1.w + 1.0 );\n"
"  xp2.xy = osg_Viewport * 0.5 * ( xp2.xy / xp2.w + 1.0 );\n"
"  gl_PointSize = 1.41 * length( xp1.xy - xp2.xy );\n"
"  pointSize = gl_PointSize;\n"
"}\n";