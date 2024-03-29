static const char SSAOPOINTSPHEREFRAG_TRACE_FRAG[] =
"uniform float radius;\n"
"uniform vec2 osg_Viewport;\n"
"uniform int ssao;\n"
"uniform int shade;\n"
"uniform float occlusionFactor;\n"
//"float occlusionFactor = 0.5;\n"
"varying vec4 color;\n"
"varying vec3 rayorigin;\n"
"vec3 raydir;\n"
"varying vec3 spherepos;\n"
"float sphereradsq;\n"
"varying float pointSize;\n"
"\n"
"#extension GL_ARB_texture_rectangle : enable\n"
"uniform sampler2DRect depthMap;\n"
"uniform float numSamples; //number of rays\n"
"uniform float hwMax; //max pixels\n"
"uniform float dstep; //step multiplier\n"
//"float numSamples = 16.0; //number of rays\n"
//"float hwMax = 20.0; //max pixels\n"
//"float dstep = 1.0; //step multiplier\n"
"varying float R; // radius in world coordinates\n"
"varying float pixelRadius; //radius in screen coordinates\n"
//"float R=10.0; // radius in world coordinates\n"
//"float pixelRadius=20.0; //radius in screen coordinates\n"
"float width = osg_Viewport.x;\n"
"float height = osg_Viewport.y;\n"
"vec3 screenPosition;\n"
"vec3 normal;\n"
"vec3 worldPosition;\n"
"//-----------------------------------------------------------------------------\n"
"vec3 ssUnproject( vec3 v )\n"
"{\n"
"  vec4 p = vec4( v, 1.0 );\n"
"  p.x /= width;\n"
"  p.y /= height;\n"
"  p.xyz -= 0.5;\n"
"  p.xyz *= 2.0;\n"
"  p = gl_ProjectionMatrixInverse * p;\n"
"  p.xyz /= p.w;\n"
"  return p.xyz;\n"
"}\n"
"\n"
"//------------------------------------------------------------------------------\n"
"//cosine of mininum angle used for angle occlusion computation (~30 deg. best)\n"
"const float minCosAngle = 0.2; // ~78 deg. from normal, ~22 deg from tangent plane\n"
"// adjusted pixel radius\n"
"float PR = 1.0;\n"
"// adjusted world space radius\n"
"float PRP = 1.0;\n"
"// attenuation coefficient\n"
"float B = 0.0;\n"
"\n"
"void ComputeRadiusAndOcclusionAttenuationCoeff()\n"
"{\n"
"  // clamp pixel radius to [1.0, max radius]  \n"
"  PR = clamp( pixelRadius, 1.0, hwMax );\n"
"  // set coefficient for occlusion attenuation function 1 / ( 1 + B * dist^2 )\n"
"  // Rw : Rs = NRw : NRs => NRw = Rw * NRs / Rs \n"
"  PRP = R * PR / pixelRadius; // clamp world space radius\n"
"  // the coefficient is set such that the occlusion at maximum distance\n"
"  // is ~ 1/2 of the peak value\n"
"  B = ( 1.0 - 0.1 ) / ( PRP * PRP );\n"
"}\n"
"struct I\n"
"{\n"
"  vec3 P;\n"
"  vec3 N;\n"
"  float t;  \n"
"};\n"
"\n"
"vec3 refcolor = vec3( .1, .1, .1 );\n"
"vec3 lightDir = vec3( 0., 0., -1. );\n"
"// returns occlusion at pixel x, y;\n"
"//------------------------------------------------------------------------------\n"
"// occlusion function for horizontal (y=y0) lines\n"
"float hocclusion( float ds )\n"
"{\n"
" \n"
"  // follow line y = y0 in screen space along negative or positive direction according to ds\n"
"  // where:\n"
"  //    x0 = screenPosition.x\n"
"  //    y0 = screenPosition.y  \n"
"  // each point visible from current fragment (i.e. with z < current z)\n"
"  // will contribute to the occlusion factor\n"
"  int occSteps = 0;  // number of occlusion rays \n"
"  vec3 p = screenPosition;\n"
"  float occl = 0.0; // occlusion\n"
"  float z = 1.0; // z in depth map\n"
"  float dz = 0.; //\n"
"  float dist = 1.; // distance between current point and shaded point \n"
"  float prev = 0.; // previous angular coefficient \n"
"  vec3 I; // vector from point in depth map to shaded point\n"
"  int upperI = int( PR / abs( ds ) );\n"
"  for( int i = 0; i != upperI; ++i )\n"
"  {\n"
"      p.x += ds;\n"
"    z = texture2DRect( depthMap, p.xy ).x;\n"
"    // compute angular coefficient: if angular coefficient\n"
"    // is greater than last computed coefficient it means the point is \n"
"    // visible from screenPosition and therefore occludes it \n"
"    dz = screenPosition.z - z;\n"
"    dist = distance( p.xy, screenPosition.xy );\n"
"    float angCoeff = dz / dist;\n"
"    if( angCoeff > prev )    \n"
"    {\n"
"      p.z = z;\n"
"      prev = angCoeff;\n"
"      // compute vector from point in world coordinates to point whose\n"
"      // occlusion is being computed\n"
"      vec3 X = ssUnproject( p ) - worldPosition.xyz;\n"
"      // ADD AO contribution: function of angle between normal and ray\n"
"      float k = dot( normal, normalize( X ) );\n"
"      // if occluded compute occlusion and increment number of occlusion\n"
"      // contributions\n"
"      if( k > minCosAngle )\n"
"      {\n"
"           occl += ( k / ( 1. + B * dot( X, X ) ) );\n"
"           ++occSteps;\n"
"      }\n"
"    }\n"     
"  }\n"
"  // return average occlusion along ray: divide by number of occlusions found \n"
"  return occl / max( 1.0, float( occSteps ) );\n"
"}\n"
"//------------------------------------------------------------------------------\n"
"// occlusion function for vertical (x=x0) lines\n"
"float vocclusion( float ds )\n"
"{\n"
" \n"
"  // follow line x = x0 in screen space along negative or positive direction according to ds\n"
"  // where:\n"
"  //    x0 = screenPosition.x\n"
"  //    y0 = screenPosition.y\n"  
"  // each point visible from current fragment (i.e. with z < current z)\n"
"  // will contribute to the occlusion factor\n"
"  int occSteps = 0;  // number of occlusion rays \n"
"  vec3 p = screenPosition;\n"
"  float occl = 0.0; // occlusion\n"
"  float z = 1.0; // z in depth map\n"
"  float dz = 0.; //\n"
"  float dist = 1.; // distance between current point and shaded point \n"
"  float prev = 0.; // previous angular coefficient \n"
"  vec3 I; // vector from point in depth map to shaded point\n"
"  int upperI = int( PR / abs( ds ) );\n"
"  for( int i = 0; i != upperI; ++i )\n"
"  {\n"
"    p.y += ds; \n"
"    z = texture2DRect( depthMap, p.xy ).x;\n"
"    // compute angular coefficient: if angular coefficient\n"
"    // is greater than last computed coefficient it means the point is \n"
"    // visible from screenPosition and therefore occludes it \n"
"    dz = screenPosition.z - z;\n"
"    dist = distance( p.xy, screenPosition.xy );\n"
"    float angCoeff = dz / dist;\n"
"    if( angCoeff > prev )  \n"    
"    {\n"
"      p.z = z;\n"
"      prev = angCoeff;\n"
"      // compute vector from point in world coordinates to point whose\n"
"      // occlusion is being computed\n"
"      vec3 X = ssUnproject( p ) - worldPosition.xyz;\n"
"      // ADD AO contribution: function of angle between normal and ray\n"
"      float k = dot( normal, normalize( X ) );\n"
"      // if occluded compute occlusion and increment number of occlusion\n"
"      // contributions\n"
"      if( k > minCosAngle )\n"
"      {\n"
"        occl += ( k / ( 1. + B * dot( X, X ) ) );\n"
"        ++occSteps;\n"
"      }\n"
"    } \n"     
"  }\n"
"  // return average occlusion along ray: divide by number of occlusions found \n"
"  return occl / max( 1.0, float( occSteps ) );\n"
"}\n"
"//------------------------------------------------------------------------------\n"
"// occlusion function for non degenerate (i.e. lines not paralles to x or y axis) lines\n"
"float occlusion( vec2 dir )\n"
"{\n"
"  // compute angular coefficient and steps\n"
"  float m = dir.y / dir.x;\n"
"  \n"
"  // follow line y = m * ( x - x0 ) + y0 in screen space\n"
"  // where:\n"
"  //    m = dir.y / dir.x\n"
"  //    x0 = screenPosition.x\n"
"  //    y0 = screenPosition.y  \n"
"  // each point visible from current fragment (i.e. with z < current z)\n"
"  // will contribute to the occlusion factor\n"
"  int occSteps = 0;  // number of occlusion rays \n"
"  vec3 p = screenPosition;\n"
"  // compute number of i (x) steps\n"
"  // size of radius = sqrt( (num x steps)^2 + (ang. coeff. * num x steps)^2 )\n"
"  int upperI = int( PR * inversesqrt( 1.0 + m * m ) / abs( dstep ) );\n"
"  float occl = 0.0; // occlusion\n"
"  float z = 1.0; // z in depth map\n"
"  float dz = 0.; //\n"
"  float dist = 1.; // distance between current point and shaded point \n"
"  float prev = 0.; // previous angular coefficient \n"
"  //vec3 I; // vector from point in depth map to shaded point\n"
"  float ds = sign( dir.x ) * dstep;\n"
"  for( int i = 0; i != upperI; ++i )\n"
"  {\n"
"    p.x += ds;\n"
"    p.y += ds * m;\n"
"    z = texture2DRect( depthMap, p.xy ).x;\n"
"    // compute angular coefficient: if angular coefficient\n"
"    // is greater than last computed coefficient it means the point is \n"
"    // visible from screenPosition and therefore occludes it \n"
"    dz = screenPosition.z - z;\n"
"    dist = distance( p.xy, screenPosition.xy );\n"
"    float angCoeff = dz / dist;\n"
"    if( angCoeff > prev )\n"      
"    {\n"
"      p.z = z;\n"
"      prev = angCoeff;\n"
"      // compute vector from point in world coordinates to point whose\n"
"      // occlusion is being computed\n"
"      vec3 X = ssUnproject( p ) - worldPosition.xyz;\n"
"      // ADD AO contribution: function of angle between normal and ray\n"
"      float k = dot( normal, normalize( X ) );\n"
"      // if occluded compute occlusion and increment number of occlusion\n"
"      // contributions\n"
"      if( k > minCosAngle )\n"
"      {\n"
"        occl += ( k / ( 1. + B * dot( X, X ) ) );\n"
"        ++occSteps;\n"
"      }\n"
"    }\n"      
"  }\n"
"\n"
"  // return average occlusion along ray: divide by number of occlusions found \n"
"  return occl / max( 1.0, float( occSteps ) );\n"
"}\n"
"//------------------------------------------------------------------------------\n"
"float ComputeOcclusion()\n"
"{\n"
"    // numSamples is an upper limit fr the number of directions == number of rays\n"
"    // if the 2D (i,j) sample points are laid out on the edges of a square\n"
"    // and numSamples is the total number of points then the number of samples on one edge\n"
"    // is numSamples / 4;\n"
"    // the i and j indices are assumed to be in the range:\n"
"    //  [-(numSamples / 4) / 2, +(numSamples / 4) / 2] == [ -numSamples/8,+numSamples/8 ]\n"
"    int hw = int( max( 1.0, numSamples / 8.0 ) ); \n"	
"    // ppos is [x pixel, y pixel, depth (0..1) ]\n"
"    float occ = 0.0;\n"
"    int i = -hw;\n"
"    int j = -hw;\n"
"    // vertical edges, j = 0 excluded\n"
"    for( ; j != 0; ++j ) occ += occlusion( vec2( float( i ), float( j ) ) );\n"
"    for( j = 1; j != hw + 1; ++j ) occ += occlusion( vec2( float( i ), float( j ) ) );\n"
"    i = hw;\n"
"    for( ; j != hw + 1; ++j ) occ += occlusion( vec2( float( i ), float( j ) ) );\n"
"    for( j = 1; j != hw + 1; ++j ) occ += occlusion( vec2( float( i ), float( j ) ) );\n"
"    // horizontal edges, i = 0 excluded\n"
"    j = -hw;\n"
"    for( i = -hw + 1; i != 0; ++i ) occ += occlusion( vec2( float( i ), float( j ) ) );\n"
"    for( i = 1; i != hw; ++i ) occ += occlusion( vec2( float( i ), float( j ) ) );\n"
"    j = hw;\n"
"    for( i = -hw + 1; i != 0; ++i ) occ += occlusion( vec2( float( i ), float( j ) ) );\n"
"    for( i = 1; i != hw; ++i ) occ += occlusion( vec2( float( i ), float( j ) ) ); \n"
"    // i = 0 and j = 0\n"
"    occ += hocclusion( -dstep );\n"
"    occ += hocclusion( dstep );\n"
"    occ += vocclusion( -dstep );\n"
"    occ += vocclusion( dstep );\n"
"    // divide occlusion by the number of shot rays\n"
"    occ /= max( 1.0, float( 8 * hw - 2 ) );\n"
"    return occ;\n"
"}\n"
"// computes intersection of ray with sphere\n"
"I ComputeRaySphereIntersection( vec3 C, vec3 D )\n"
"{\n"
"  I i;\n"
"  i.t = -1.0;\n"
"  vec3 P = C - spherepos;\n"
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
"  i.P = C + t * D;\n"
"  i.N = normalize( i.P - spherepos );\n"
"  i.t = t;\n"
"  return i;\n"
"}\n"
"float kd = 1.0;\n"
"float ka = 0.0;\n"
"float ks = 0.0;\n"
"float sh = 90.0;\n"
"vec4 ComputeColor( vec3 n )\n"
"{\n"
"  return vec4( color.rgb * dot( normalize( raydir ), -n ), color.a );\n"
"  vec3 N = faceforward( n, lightDir, n );\n"
"  float d = dot( N, -normalize( lightDir ) );\n"
"  float s = pow( max( 0.0, dot( vec3( 0, 0, 1 ), reflect( lightDir, N ) ) ), sh );\n"
"  if( length( color.rgb ) < 0.75 ) ks = 0.;\n"
"  return vec4(  ks * s * refcolor + kd * d * color.rgb + ka * color.rgb, color.a );\n"
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
"  if( bool( shade ) ) gl_FragColor = ComputeColor( i.N );\n"
"  float z = dot( vec4( i.P, 1 ), gl_ProjectionMatrixTranspose[ 2 ] );\n"
"  float w = dot( vec4( i.P, 1 ), gl_ProjectionMatrixTranspose[ 3 ] );\n"
"  gl_FragDepth = 0.5 * ( z / w + 1.0 );\n"
"  if( bool( ssao ) && bool( shade ) )\n"
"  {\n"
"    normal = i.N;\n"
"    worldPosition = i.P;\n"
"    screenPosition = vec3( gl_FragCoord.xy, gl_FragDepth );\n"
"    ComputeRadiusAndOcclusionAttenuationCoeff();\n"
"    gl_FragColor.rgb *= 1.0 - smoothstep( 0., 1., ComputeOcclusion() * occlusionFactor );\n"
"  }\n"
"}\n";
