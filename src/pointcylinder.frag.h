static const char POINTCYLINDERFRAG[] =
"vec3 raydir;\n"
"varying vec3 rayorigin;\n"
"varying vec3 ccenter;\n"
"varying vec3 tyaxis;\n"
"uniform float hlength;\n"
"varying vec4 color;\n"
"uniform vec2 osg_Viewport;\n"
"varying float a;\n"
"varying float b;\n"
"varying float c;\n"
"varying float d;\n"
"varying float e;\n"
"varying float f;\n"
"varying float g;\n"
"varying float h;\n"
"varying float i;\n"
"varying float j;\n"
"varying float pointSize;\n"
"struct I\n"
"{\n"
"  vec3 P;\n"
"  vec3 N;\n"
"  float t;\n"
"};\n"
"vec3 ComputeNormal( vec3 P )\n"
"{\n"
"  vec3 N;\n"
"  N.x = a*P.x + d*P.y + e*P.z + g;\n"  // should multiply by 2
"  N.y = b*P.y + d*P.x + f*P.z + h;\n"  // should multiply by 2   
"  N.z = c*P.z + e*P.x + f*P.y + i;\n"  // should multiply by 2
"  return normalize( N );\n"
"}\n"
"I ComputeRayQuadricIntersection()\n"
"{\n"
"  I ip;\n"
"  ip.t = -1.0;\n"
"  vec3 P = rayorigin - ccenter;\n"
"  vec3 D = raydir;\n"
"  float A = a*D.x*D.x + b*D.y*D.y + c*D.z*D.z + 2.*( d*D.x*D.y + e*D.x*D.z + f*D.y*D.z );\n"
"  float B = 2.*( a*P.x*D.x + b*P.y*D.y + c*P.z*D.z + d*D.x*P.y + d*D.y*P.x + e*P.x*D.z +\n"
"            e*P.z*D.x + f*P.y*D.z + f*P.z*D.y + g*D.x + h*D.y + i*D.z );\n"
"  float C = a*P.x*P.x + b*P.y*P.y + c*P.z*P.z +\n"
"            2.*( d*P.x*P.y + e*P.x*P.z + f*P.y*P.z + g*P.x + h*P.y +i*P.z ) + j;\n"
"  float delta = B * B - 4. * A * C;\n"
"  if( delta < 0.0 ) return ip;\n"
"  float d = sqrt( delta );\n"
"  A = 1. / A;\n"
"  A *= 0.5;\n"
"  float t2 = A * ( -B + d );\n"
"  float t1 = A * ( -B - d );\n"
"  vec3 P1 = rayorigin + D * min( t1, t2 );\n"
"  vec3 P2 = rayorigin + D * max( t1, t2 );\n"
"  if( abs( dot( P1 - ccenter, tyaxis ) ) <= hlength )\n"
"  {\n"
"    ip.P = P1;\n"
"    ip.N = ComputeNormal( P1 - ccenter );\n"
"    ip.t = 0.;\n"
"    return ip; \n"
"  }\n"
"  if( abs( dot( P2 - ccenter, tyaxis ) ) <= hlength )\n"
"  {\n"
"    ip.P = P2;\n"
"    ip.N = ComputeNormal( P2 - ccenter );\n"
"    ip.t = 0.;\n"
"    return ip; \n"
"  }\n"
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
"  if( pointSize < 8.0 ) return color;\n"
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
"  I i = ComputeRayQuadricIntersection();\n"
"  if( i.t < 0.0 ) discard;\n"
"  gl_FragColor = ComputeColor( color, i.N );\n"
"  float z = dot( vec4( i.P, 1 ), gl_ProjectionMatrixTranspose[ 2 ] );\n"
"  float w = dot( vec4( i.P, 1 ), gl_ProjectionMatrixTranspose[ 3 ] );\n"
"  gl_FragDepth = 0.5 * ( z / w + 1.0 );\n"
"}\n";
