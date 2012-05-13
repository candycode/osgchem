varying vec3 raydir;
varying vec3 rayorigin;
varying vec3 oglcolor;
varying vec3 ccenter;
varying vec3 tyaxis;
uniform float hlength;
uniform vec4 bcolor;
uniform vec4 tcolor;

varying float a;
varying float b;
varying float c;
varying float d;
varying float e;
varying float f;
varying float g;
varying float h;
varying float i;
varying float j;


struct I
{
  vec3 P;
  vec3 N;
  float t;
  bool top;  
};

vec3 ComputeNormal( vec3 P )
{
  vec3 N;
  N.x = 2*a*P.x + 2*d*P.y + 2*e*P.z + 2*g;
  N.y = 2*b*P.y + 2*d*P.x + 2*f*P.z + 2*h;
  N.z = 2*c*P.z + 2*e*P.x + 2*f*P.y + 2*i;
  return normalize( N ); 		
}

I ComputeRayQuadricIntersection()
{
  I ip;
  ip.t = -1;
  vec3 P = rayorigin - ccenter;
  vec3 D = raydir;
  float A = a*D.x*D.x + b*D.y*D.y + c*D.z*D.z + 2*d*D.x*D.y + 2*e*D.x*D.z + 2*f*D.y*D.z;
  float B = 2*( a*P.x*D.x + b*P.y*D.y + c*P.z*D.z + d*D.x*P.y + d*D.y*P.x + e*P.x*D.z +
            e*P.z*D.x + f*P.y*D.z + f*P.z*D.y + g*D.x + h*D.y + i*D.z );
  float C = a*P.x*P.x + b*P.y*P.y + c*P.z*P.z +
            2*( d*P.x*P.y + e*P.x*P.z + f*P.y*P.z + g*P.x + h*P.y +i*P.z ) + j;
  float delta = B * B - 4.0 * A * C;
  if( delta < 0.0 ) return ip;
  float d = sqrt( delta );
  A = 1. / A;
  float t2 = .5 * A * ( -B + d );
  float t1 = .5 * A * ( -B - d );
  vec3 P1 = rayorigin + D * min( t1, t2 );
  vec3 P2 = rayorigin + D * max( t1, t2 );
  if( abs( dot( P1 - ccenter, tyaxis ) ) <= hlength )
  {
    ip.P = P1;
    ip.N = ComputeNormal( P1 - ccenter );
    ip.t = 0;
    ip.top = dot( P1 - ccenter, tyaxis ) > 0.0;
    return ip; 
  }
  if( abs( dot( P2 - ccenter, tyaxis ) ) <= hlength )
  {
    ip.P = P2;
    ip.N = ComputeNormal( P2 - ccenter );
    ip.t = 0;
    ip.top = dot( P2 - ccenter, tyaxis ) > 0.0;
    return ip; 
  }  
  return ip;         
}

vec3 lightDir = vec3( 0, 0, -1 );
float kd = 1.0;
float ka = 0.01;
float ks = .5;
float sh = 90.0;
vec3 refcolor = vec3( 1, 1, 1 );


vec4 ComputeColor( vec4 color, vec3 n )
{
  vec3 N = faceforward( n, lightDir, n );
  float d = dot( N, -normalize( lightDir ) );
  float s = pow( max( 0.0, dot( vec3( 0, 0, 1 ), reflect( lightDir, N ) ) ), sh );
  return vec4(  ks * s * refcolor + kd * d * color.rgb + ka * color.rgb, color.a );
}

// entry point 
void main(void)
{
  I i = ComputeRayQuadricIntersection();
  if( i.t < 0.0 ) discard;
  if( i.top ) gl_FragColor = ComputeColor( tcolor, i.N );
  else gl_FragColor = ComputeColor( bcolor, i.N );
  float z = dot( vec4( i.P, 1 ), gl_ProjectionMatrixTranspose[ 2 ] );
  float w = dot( vec4( i.P, 1 ), gl_ProjectionMatrixTranspose[ 3 ] );
  gl_FragDepth = 0.5 * ( z / w + 1.0 );
}


