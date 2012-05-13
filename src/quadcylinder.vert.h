static const char QUADCYLINDERVERT[] =
"#version 120\n"
"varying vec3 rayorigin;\n"
"varying vec4 color;\n"
"varying float pointSize;\n"
"varying float L;\n"
"varying vec3 D;\n"
"varying vec3 Dn;\n"
"varying vec3 P0;\n"
"varying vec3 P1;\n"
"varying float sqradius;\n"
"varying vec3 raydir;\n"
"uniform float radius;\n"
"uniform mat4 osg_ViewMatrix;\n"
"void main()\n"
"{\n"
"  sqradius = radius * radius;\n"
"  color = gl_Color;\n"
"  D = mat3( osg_ViewMatrix ) * gl_Normal;\n"
"  vec3 taxis = normalize( D );\n"
"  Dn = taxis;\n"
"  L = length( D );\n"
"  bool perspective = gl_ProjectionMatrix[ 3 ][ 3 ] < 0.001 && gl_ProjectionMatrix[ 2 ][ 3 ] != 0.0;\n"
"  vec4 p = gl_ModelViewMatrix * gl_Vertex;\n"
"  vec4 center = gl_ModelViewMatrix * vec4( 0., 0., 0., 1. );\n"
"  vec3 ep = center.xyz / center.w;\n"
"  P0 = ep - 0.5 * L * Dn;\n"
"  P1 = ep + 0.5 * L * Dn;\n"
"  gl_ClipVertex = p;\n"
"  if( perspective )\n"
"  {\n"
"    raydir = vec3( p ) / p.w;\n"
"    rayorigin = vec3( 0., 0., 0. );\n"
"  }\n"
"  else\n"
"  {\n"
"    raydir = vec3( 0., 0., -1. );\n"
"    rayorigin = vec3( p.x / p.w, p.y / p.w, 0 );\n"
"  }\n"
"  gl_Position = gl_ProjectionMatrix * p;\n"
"}\n";