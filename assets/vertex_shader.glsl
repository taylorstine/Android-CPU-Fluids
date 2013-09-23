uniform mat4 mvp_matrix;

attribute vec3 vPosition;
attribute vec3 vNormal;
attribute vec3 vColor;
attribute vec2 vTexcoord;

varying vec3 icolor;
varying vec3 inormal;
varying vec2 itexcoord;

void main(){
  icolor = vColor;
  inormal = vNormal;
  itexcoord = vTexcoord;
  gl_Position = mvp_matrix * vec4( vPosition, 1.0 );
}
