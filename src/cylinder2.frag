uniform mat4 objectMatrix;

varying vec3 V;
varying vec3 spherepos;
varying vec3 rayorigin;
//varying float sphereradsq;
varying vec3 oglcolor;
varying vec3 ccenter;

struct I
{
  vec3 P;
  vec3 N;
  float t;  
};

const float sphereradsq = 0.2 * 0.2;

// computes intersection of ray with cylinder centered at the origin and aligned
// along the y axis 
I ComputeRayCylinderIntersection( vec3 R, vec3 D )
{
  I i;
  i.t = -1.0;
  vec3 P = R - ccenter;
  const float a = D.x * D.x + D.z * D.z;
  const float b = 2.0 * ( D.x * P.x + D.z * P.z );
  const float c = P.x * P.x + P.z * P.z - 0.2 * 0.2;
  const float delta = b * b - 4. * a * c;
  if( delta < 0.0 ) return i;
  const float d = sqrt( delta );
  const float t2 = ( -b + d ) / ( 2.0 * a );
  const float t1 = ( -b - d ) / ( 2.0 * a );
  const float t = min( t1, t2 );
  i.P = P + D * t;
  i.N = normalize( vec3( i.P.x, 0, i.P.z ) );
  i.t = t;
  return i;
}

// entry point 
void main(void)
{
  vec3 raydir = normalize( V );
  I i = ComputeRayCylinderIntersection( rayorigin, raydir );
  if( i.t < 0.0 ) discard;//gl_FragColor.g = 1;
  //else gl_FragColor.r = 1;
  vec3 N = i.N;
  float d = dot( N, vec3( 0, 0, 1 ) );
  gl_FragColor = vec4( d, d, d, 1 );
  vec4 P = gl_ModelViewProjectionMatrix * vec4( i.P, 1 );
  P /= P.w;
  gl_FragDepth = 0.5 * ( P.z + 1.0 );
  //vec4 P = gl_ModelViewProjectionMatrix * objectMatrix * ( vec4( i.P, 1 ) );// * objectMatrix ); //objectMatrix is row major
  //vec4 N = gl_ModelViewMatrix * objectMatrix * ( vec4( i.N, 1 ) );//* objectMatrix ); //objectMatrix is row major 
  //P /= P.w;
  //N /= N.w;
  //gl_FragDepth = 0.5 * ( P.z + 1.0 );
  //float d = dot( vec3( i.N ), vec3( 0, 0, 1 ) );
  //if( i.t < 0.0 ) gl_FragColor = vec4( 1, 0, 0, 1 );
  //else gl_FragColor = vec4( d, d, d, 1 );//vec4( d, d, d, 1 );
}


