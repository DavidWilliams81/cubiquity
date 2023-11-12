#version 330 core
out vec4 FragColor;

in vec2 TexCoords;

uniform sampler2D screenTexture;

void main()
{
	int taps = 9; // Must be odd.
	int halfTaps = taps / 2; // Rounds down
    int texWidth = textureSize(screenTexture, 0).x;
	float offset = float(1.0 / texWidth);
	float count = 1.0;
    vec4 col = texture(screenTexture, TexCoords);
	
	vec4 sample;
	for(int tap = 1; tap <= halfTaps; tap++)
	{
		sample = texture(screenTexture, TexCoords + vec2((-offset) * tap, 0));
		if(abs(col.a - sample.a) < 0.5) { col.rgb += sample.rgb; count += 1.0;}

		sample = texture(screenTexture, TexCoords + vec2((+offset) * tap, 0));
		if(abs(col.a - sample.a) < 0.5) { col.rgb += sample.rgb; count += 1.0;}
	}
	
	FragColor = vec4(col.rgb / count, col.a);
}
