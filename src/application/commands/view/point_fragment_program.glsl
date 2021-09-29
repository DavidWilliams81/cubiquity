#version 330 core

in vec4 modelSpacePosition;
in vec4 worldSpacePosition;
in vec4 worldSpaceCentre;
in vec4 voxelNormal;

uniform vec3 cameraPos;

out vec4 color;

bool intersectPlane(vec3 n, vec3 p0, vec3 l0, vec3 l, inout float t) 
{ 
    // assuming vectors are all normalized
    float denom = dot(n, l); 
    if (denom > 1e-6) { 
        vec3 p0l0 = p0 - l0; 
        t = dot(p0l0, n) / denom; 
        return (t >= 0); 
    } 
 
    return false; 
}

void main()
{

	vec3 rayDir = normalize(worldSpacePosition.xyz - cameraPos);
	float distToSurface = length(rayDir);
	
	float t = 0.0f;
	bool hit = intersectPlane(-voxelNormal.xyz, worldSpaceCentre.xyz, cameraPos, rayDir, t);
	
	vec3 intersection = cameraPos + (rayDir * t);

	float dist = length(intersection.xyz - worldSpaceCentre.xyz);
	// For visualization we derive the color of the voxel from its normal.
	vec4 voxelColor = vec4(voxelNormal.xyz * 0.5 + vec3(0.5, 0.5, 0.5), 1.0);
	//vec4 voxelColor = vec4(1.0f, 1.0f, 1.0f, 1.0f);
	
	if(hit)
		color.rgb = voxelColor.xyz;
	else
		color.rgb = vec3(1.0f, 0.0f, 0.0f);
	
	color.a = 1.0 - dist;
}