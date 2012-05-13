uniform mat4 objectMatrix;
uniform mat4 invObjectMatrix;

varying vec3 RayP;
varying vec3 RayD;

void main()
{
  const mat4 M = invObjectMatrix *  gl_ModelViewMatrixInverse;
  RayP =  M * vec4( 0, 0, 0, 1 );
  RayD =  M * gl_Vertex;
  gl_Position = ftransform();
} 