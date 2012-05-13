uniform float radius;
uniform vec2 vport;

vec3 raydir;
vec3 spherepos;
vec3 rayorigin;
float sphereradsq;
vec4 color;




struct I
{
  vec3 P;
  vec3 N;
  float t;  
};

// computes intersection of ray with sphere
I ComputeRaySphereIntersection( vec3 R, vec3 D )
{
  I i;
  i.t = -1;
  vec3 P = R - spherepos;
  float a = D.x * D.x + D.y * D.y + D.z * D.z;
  float b = 2.0 * ( D.x * P.x + D.y * P.y + D.z * P.z );
  float c = P.x * P.x + P.y * P.y + P.z * P.z - sphereradsq;
  float delta = ( b * b - 4. * a * c );
  if( delta < 0.0 ) return i;
  float d = sqrt( delta );
  a = .5 / a;
  float t2 = ( -b + d ) * a;
  float t1 = ( -b - d ) * a;
  float t = min( t1, t2 );
  if( t < 0.0 ) return i;
  i.P = R + t * D;
  i.N = normalize( i.P - spherepos );
  i.t = t;
  return i;
}

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
  color = gl_Color;
      	  
  spherepos = gl_ModelViewMatrix * vec4( 0, 0, 0, 1 );
  sphereradsq = radius * radius;
        
  bool perspective = gl_ProjectionMatrix[ 3 ][ 3 ] < 0.001 && gl_ProjectionMatrix[ 2 ][ 3 ] != 0.0;
  
  vec3 fc = vec3( gl_FragCoord );// * gl_FragCoord.w;
  fc.x = gl_FragCoord.x / vport.x;
  fc.y = gl_FragCoord.y / vport.y;
  fc.z = gl_FragCoord.z;
  fc *= 2.0;
  fc -= 1.0;
  vec4 p = gl_ProjectionMatrixInverse * vec4( fc, 1. );
       
  if( perspective )
  {
    raydir = vec3( p ) / p.w;
    rayorigin = vec3( 0, 0, 0 );
  }  
  else
  {
    raydir = vec3( 0, 0, -1 );
    rayorigin = vec3( p.x / p.w, p.y / p.w, 0 );
  }  
  I i = ComputeRaySphereIntersection( rayorigin, raydir );
  if( i.t < 0.0 ) discard;
  gl_FragColor = ComputeColor( i.N );	  
  float z = dot( vec4( i.P, 1 ), gl_ProjectionMatrixTranspose[ 2 ] );
  float w = dot( vec4( i.P, 1 ), gl_ProjectionMatrixTranspose[ 3 ] );
  gl_FragDepth = 0.5 * ( z / w + 1.0 );
}


