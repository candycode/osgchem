
uniform mat4 invObjectMatrix;

varying vec3 V;
varying vec3 rayorigin;
//varying vec3 oglcolor;
varying vec3 ccenter;

void main()
{
  vec4 o = gl_ModelViewMatrixInverse * ( vec4( 0, 0, 0, 1 ) * invObjectMatrix );
  rayorigin = vec3( o ) / o.w;
  
  ccenter = gl_ModelViewMatrix * vec4( 0, 0, 0, 1 );
  
  vec4 p = gl_Vertex * invObjectMatrix * gl_ModelViewMatrix;
  gl_ClipVertex = p;
  rayorigin = vec3( 0, 0, 0 );  
  V = normalize( vec3( p ) / p.w );// - rayorigin);
  
  gl_Position = ftransform();
} 