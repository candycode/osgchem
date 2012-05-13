
uniform mat4 Q;
uniform mat4 osg_ViewMatrixInverse;
uniform mat4 osg_ViewMatrix;
uniform vec3 axis;

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

varying vec3 raydir;
varying vec3 rayorigin;
varying vec3 ccenter;
varying vec3 tyaxis;

void main()
{
  
  mat4 VMI = mat4( mat3( osg_ViewMatrixInverse ) );
  mat4 VMIT = transpose( VMI );
  mat4 M = VMIT * Q * VMI;
  
  a = M[ 0 ][ 0 ];
  b = M[ 1 ][ 1 ];
  c = M[ 2 ][ 2 ];
  d = M[ 1 ][ 0 ];
  e = M[ 2 ][ 0 ];
  f = M[ 2 ][ 1 ];
  g = M[ 3 ][ 0 ];
  h = M[ 3 ][ 1 ];
  i = M[ 3 ][ 2 ];
  j = M[ 3 ][ 3 ];
  	  
  ccenter = gl_ModelViewMatrix * vec4( 0, 0, 0, 1 );
  tyaxis = normalize( mat3( osg_ViewMatrix ) * axis ); // should be transpose( inverse( osg_ViewMatrix ) ) * axis;
  
  bool perspective = gl_ProjectionMatrix[ 3 ][ 3 ] < 0.001 && gl_ProjectionMatrix[ 2 ][ 3 ] != 0.0;
  
  vec4 p = gl_ModelViewMatrix * gl_Vertex;
  gl_ClipVertex = p;
   
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
  gl_Position = gl_ProjectionMatrix * p;
} 