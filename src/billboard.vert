
uniform float radius;

varying vec3 V;
varying vec3 spherepos;
varying vec3 rayorigin;
varying float sphereradsq;
varying vec3 oglcolor;

void main()
{
  vec4 ecpos = gl_ModelViewMatrix * gl_Vertex;
  gl_ClipVertex = ecpos;
  oglcolor = vec3( gl_Color );

  vec4 spos = gl_ModelViewMatrix * vec4( 0, 0, 0, 1.0 );
  spherepos = vec3( spos ) / spos.w;

  vec4 rspos  = gl_ModelViewMatrix * vec4( vec3( 0, 0, 0 ) + vec3( radius, 0, 0 ), 1 ); 	
  sphereradsq = length( spherepos - ( vec3( rspos ) / rspos.w ) );
  sphereradsq *= sphereradsq;
  
  const bool projectionmode = true;
  if ( projectionmode )
  {
     // set view direction vector from eye coordinate of vertex, for 
     // perspective views
     V = normalize( vec3( ecpos ) / ecpos.w );
     rayorigin = vec3( 0, 0, 0 );
   } else {
     // set view direction vector with constant eye coordinate, used for
     // orthographic views
     V = vec3( 0.0, 0.0, -1.0 );
     rayorigin = vec3( ( ecpos.xy / ecpos.w ), 0.0 );  
  }   
  gl_Position = ftransform();
} 