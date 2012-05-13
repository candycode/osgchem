static const char POINTSPHEREVERT[] =
"varying vec3 spherepos;\n"
"varying vec4 color;\n"
"uniform vec2 osg_Viewport;\n"
"uniform float radius;\n"
"varying vec3 rayorigin;\n"
"varying float pointSize;\n"
"void main(void)\n"
"{\n"
"  color = gl_Color;\n"
"  vec4 P = gl_ModelViewMatrix * gl_Vertex;\n"
"  spherepos = vec3( P ) / P.w;\n"
"  gl_Position = gl_ProjectionMatrix * P;\n"
"  bool perspective = gl_ProjectionMatrix[ 3 ][ 3 ] < 0.001 && gl_ProjectionMatrix[ 2 ][ 3 ] != 0.0;\n"
"  if( perspective )\n"
"  {\n"
"    rayorigin = vec3( 0, 0, 0 );\n"
"  }\n"
"  else\n"
"  {\n"
"    rayorigin = vec3( P.x / P.w, P.y / P.w, 0 );\n"
"  }\n"
"  // compute pixel size\n"
"  vec3 axis = cross( normalize( spherepos - rayorigin ), vec3( 0., 1., 0. ) );\n"
"  vec4 xp1 = gl_ProjectionMatrix * vec4( spherepos - radius * axis, 1.0 );\n"
"  vec4 xp2 = gl_ProjectionMatrix * vec4( spherepos + radius * axis, 1.0 );\n"
"  float xv1 = osg_Viewport.x * .5 * ( xp1.x / xp1.w + 1.0 );\n" 
"  float xv2 = osg_Viewport.x * .5 * ( xp2.x / xp2.w + 1.0 );\n"
"  gl_PointSize = abs( xv1 - xv2 );\n"
"  pointSize = gl_PointSize;\n"
"}\n";
