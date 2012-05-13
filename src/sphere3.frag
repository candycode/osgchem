varying float c;
varying float maxsqlength;
varying vec4 color;
varying vec3 spherepos;
varying vec3 rayorigin;
varying vec3 raydir;
varying vec3 RS; 

vec3 lightDir = vec3( 0, 0, -1 );
float kd = 1.0;
float ka = 0.01;
float ks = .5;
float sh = 90.0;
vec3 refcolor = vec3( 1, 1, 1 );

vec4 ComputeColor( vec3 n )
{
  vec3 N = faceforward( n, lightDir, n );
  float d = dot( N, -normalize( lightDir ) );
  float s = pow( max( 0.0, dot( vec3( 0, 0, 1 ), reflect( lightDir, N ) ) ), sh );
  return vec4(  ks * s * refcolor + kd * d * color.rgb + ka * color.rgb, color.a );
}

void main(void)
{
  float a = dot( raydir, raydir );
  if( a > maxsqlength ) discard;
  float b = 2.0 * dot( RS, raydir );
  float delta = ( b * b - 4. * a * c );
  if( delta < 0.0 ) discard;
  float d = sqrt( delta );
  a = 1. / a;
  a *= .5;
  float t2 = ( -b + d ) * a;
  float t1 = ( -b - d ) * a;
  float t = min( t1, t2 );
  if( t < 0.0 ) discard;
  vec3 P = rayorigin + t * raydir;
  vec3 N = normalize( P - spherepos );
  gl_FragColor = ComputeColor( N );	  
  float z = dot( vec4( P, 1 ), gl_ProjectionMatrixTranspose[ 2 ] );
  float w = dot( vec4( P, 1 ), gl_ProjectionMatrixTranspose[ 3 ] );
  gl_FragDepth = 0.5 * ( z / w + 1.0 );
}


