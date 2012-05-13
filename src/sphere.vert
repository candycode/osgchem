uniform float radius;

varying vec3 raydir;
varying vec4 color;
varying vec3 spherepos;
varying vec3 rayorigin;
varying float sphereradsq;
varying float maxsqlength;
varying vec3 RS;
 
void main()
{
  color = gl_Color;
      	  
  spherepos = gl_ModelViewMatrix * vec4( 0, 0, 0, 1 );
  sphereradsq = radius * radius;
        
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
  
  RS = rayorigin - spherepos; 
  float sqlength = dot( RS, RS );
  maxsqlength = sphereradsq + sqlength;
    
  gl_Position = gl_ProjectionMatrix * p;
} 
