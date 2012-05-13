// computes ray-cylinder intersection and changes fragment z coordinate
// the quad is assumed to be a billboard facing the viewpoint 

uniform float   length; // cylinder height 
uniform float   sqradius; // cylinder squared radius
uniform sampler2D texture; // texture; used in place of actual lighting computation
uniform mat4 objectMatrix;

varying vec3 RayP;
varying vec3 RayD;


// computes intersection of ray with cylinder centered at the origin and aligned
// along the y axis 
vec3 ComputeRayCylinderIntersection( vec3 P, vec3 D )
{
  const float a = D.x * D.x + D.y * D.y;
  const float b = 2.0 * ( D.x * P.x + D.y * P.y );
  const float c = P.x * P.x + P.y * P.y - sqradius;
  const float delta = b * b - 4. * a * c;
  const float d = .5 * sqrt( delta ) / a;
  const float t2 = -b + d;
  const float t1 = -b - d;
  const float t = min( t1, t2 );
  return P + D * t;
}

// entry point 
void main(void)
{
  // which texture unit ?
  int texUnit = 0;	
  vec3 cylIsect = ComputeRayCylinderIntersection( RayP, RayD );
  // discard if above or below cylinder
  if( abs( cylIsect.y ) > length ) discard;
  // compute normal
  //vec3 normal = normalize( vec3( cylIsect.x, 0, cylIsect.z ) );	
  // compute texture coordinate projecting intersection point onto quad plane
  vec2 texCoord = vec2( cylIsect.x + .5, cylIsect.y / length + .5 );
  // set fragment color
  gl_FragColor = texture2D( texture, texCoord );
  
  // compute new fragment z
  
  // transform intersection point to normalized coordinates
  vec4 tip = ( gl_ModelViewProjectionMatrix * objectMatrix ) * vec4( cylIsect, 1.0 );
  // set fragment z
  gl_FragDepth = 0.5 * ( tip.z / tip.w ) + 0.5;    
}


