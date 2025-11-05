// Tunables
#define BLACK_BLEND_THRESHOLD 0.15
#define COLOR_SPEED 0.1
#define MOVEMENT_SPEED 0.1
#define BRIGHTNESS 1.5

//----------------------------------------------------------
// Smooth union + sphere
//----------------------------------------------------------
float opSmoothUnion(float d1, float d2, float k) {
    float h = clamp(0.5 + 0.5 * (d2 - d1) / k, 0.0, 1.0);
    return mix(d2, d1, h) - k * h * (1.0 - h);
}

float sdSphere(vec3 p, float s) {
    return length(p) - s;
}

//----------------------------------------------------------
// Pseudo-random generator (deterministic per index)
//----------------------------------------------------------
float offsVal(int i) {
    return fract(sin(float(i) * 43758.5453) * 0.723);
}

//----------------------------------------------------------
// Scene map (scalable number of metaballs)
//----------------------------------------------------------
float map(vec3 p) {
    float d = 2.0;

    const int NUM_SPHERES = 12; // freely change this value

    for (int i = 0; i < NUM_SPHERES; i++) {
        float o = offsVal(i);
        float time = iTime * (o - 0.5) * 2.0;

        vec3 move = sin(time * MOVEMENT_SPEED + float(i) * vec3(5.1, 6.3, 9.2))
                  * vec3(2.0, 2.0, 0.8);

        float rad = mix(0.5, 1.0, o);
        d = opSmoothUnion(sdSphere(p + move, rad), d, 0.35);
    }
    return d;
}

//----------------------------------------------------------
// Normal calc (3 taps instead of 4)
//----------------------------------------------------------
vec3 calcNormal(vec3 p) {
    float h = 5e-4;
    vec2 k = vec2(1.0, -1.0);
    return normalize(
        k.xyy * map(p + k.xyy * h) +
        k.yyx * map(p + k.yyx * h) +
        k.yxy * map(p + k.yxy * h)
    );
}

//----------------------------------------------------------
// Main image
//----------------------------------------------------------
void mainImage(out vec4 fragColor, in vec2 fragCoord) {
    vec2 uv = fragCoord / iResolution.xy;

    vec3 rayOri = vec3((uv - 0.5) * vec2(iResolution.x / iResolution.y, 1.0) * 6.0, 3.0);
    vec3 rayDir = vec3(0.0, 0.0, -1.0);

    float depth = 0.0;
    vec3 p = rayOri;

    // reduced ray steps (32 instead of 64)
    for (int i = 0; i < 32; i++) {
        p = rayOri + rayDir * depth;
        float dist = map(p);
        depth += dist;
        if (dist < 5e-5 || depth > 6.0) break;
    }

    depth = min(6.0, depth);

    vec3 n = calcNormal(p);
    float b = max(0.0, dot(n, vec3(0.577)));

    // palette
    float t = iTime * COLOR_SPEED * 3.0;
    vec3 col = (0.5 + 0.5 * cos(b + t + uv.xyx * 2.0 + vec3(0, 2, 4))) * (0.85 + b * 0.35);
    col *= exp(-depth * 0.15);

    // apply brightness
    col *= BRIGHTNESS;
    col = clamp(col, 0.0, 1.0);

    // terminal blending (safe even if iChannel0 is missing)
    vec4 terminalColor = vec4(0.0);
#ifdef CHANNEL0_DEFINED
    terminalColor = texture(iChannel0, uv);
#endif

    // fallback if iChannel0 isn’t defined or black
    float alpha = step(length(terminalColor.rgb), BLACK_BLEND_THRESHOLD);
    vec3 blendedColor = mix(terminalColor.rgb, col * 0.3, alpha);

    fragColor = vec4(blendedColor, 1.0);
}
