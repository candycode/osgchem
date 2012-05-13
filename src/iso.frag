varying raydir;
varying vec3 rayorigin;
varying vec3 oglcolor;
varying vec3 ccenter;
varying vec3 taxis;
uniform float hlength;
uniform vec4 bcolor;
uniform vec4 tcolor;

struct I
{
  vec3 P;
  vec3 N;
  float t;
  bool top;  
};


float V0 = 0.2 * 0.2;

float IsoFun( vec3 P )
{
   return P.x * P.x + P.y * P.y + P.z * P.z;
}

vec3 ComputeNormal( vec3 P )
{
   return normalize( P - ccenter );	
}

float ComputeRayIsoSurfaceIntersection( vec3 R, vec3 D, float tm, float tM )
{
  if( tm >= tM ) return tm - tM; 
  
  vec3 P = R - ccenter;
    
  float V1 = IsoFun( P + D * tm );
  float V2 = IsoFun( P + D * tM );
  if( V1 * V2 < 0. )
  {
     float t0 = ( tM - tm ) * ( V0 - V1 ) / ( V2 - V1 );
     return t0; 
  }
  else
  {
    float tm1 = tm;
    float tM1 = 0.5 * ( tM + tm );
    float tm2 = tM1;
    float tM2 = tM;
    float t1 = ComputeRayIsoSurfaceIntersection( R, D, tm1, tM1 );
    float t2 = ComputeRayIsoSurfaceIntersection( R, D, tm2, tM2 );
    float t = min( t1, t2 );
    if( t >= 0. ) return t;
    else return max( t1, t2 ); 
  }
  return -1;         
}

// entry point 
void main(void)
{
  vec3 dir = normalize( V );
  I i = ComputeRayQuadricIntersection( rayorigin, raydir );
  if( i.t < 0.0 ) discard;//gl_FragColor.g = 1;
  vec3 N = i.N;
  if( dot( N, vec3( 0, 0, 1 ) ) < 0. ) N = -N;
  
  float d = dot( N, vec3( 0, 0, 1 ) );
  if( i.top ) gl_FragColor = vec4( d * tcolor.rgb, tcolor.a );
  else gl_FragColor = vec4( d * bcolor.rgb, bcolor.a ) ;
  vec4 P = gl_ProjectionMatrix * vec4( i.P, 1 );
  P /= P.w;
  gl_FragDepth = 0.5 * ( P.z + 1.0 );
}


