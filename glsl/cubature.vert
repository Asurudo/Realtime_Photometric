#version 330 core

layout (location = 0) in vec3 inPosition;
layout (location = 1) in vec3 inNormal;
layout (location = 2) in vec2 aTexcoord;

out vec3 wp; // World position
out vec3 n;  // Normal
out vec4 c;  // Color

out vec3 worldPosition;
out vec3 worldNormal;
out vec2 texcoord;

uniform mat4 model;
uniform mat4 view;
uniform mat4 projection;
uniform mat3 normalMatrix;

void main() {
    wp = vec3(model * vec4(inPosition, 1.0));
    n = inNormal;
    // n = mat3(transpose(inverse(model))) * inNormal;
    
    c = vec4(1.0, 1.0, 1.0, 1.0);
    // if(wp.x > 0)
       // wp.x = -(wp.x);
    
    vec4 worldpos = model * vec4(inPosition, 1.0f);
	worldPosition = worldpos.xyz;
	worldNormal = normalMatrix * inNormal;
	texcoord = aTexcoord;

	// gl_Position = projection * view * worldpos;

    gl_Position = projection * view * vec4(wp, 1.0);
}
