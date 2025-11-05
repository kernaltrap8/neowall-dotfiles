// Created by inigo quilez - iq/2013
//   https://www.youtube.com/c/InigoQuilez
//   https://iquilezles.org
// I share this piece (art and code) here in Shadertoy and through its Public API, only for educational purposes.
// You cannot use, sell, share or host this piece or modifications of it as part of your own commercial or non-commercial product, website or project.
// You can share a link to it or an unmodified screenshot of it provided you attribute "by Inigo Quilez, @iquilezles and iquilezles.org".
// If you are a teacher, lecturer, educator or similar and these conditions are too restrictive for your needs, please contact me and we'll work it out.

// Optimized version - reduced iterations for better performance
#define AA 2
#define MAX_ITER 256  // Reduced from 512 for better performance

float mandelbrot(in vec2 c)
{
    float c2 = dot(c, c);
    // Early exit optimizations
    if(256.0*c2*c2 - 96.0*c2 + 32.0*c.x - 3.0 < 0.0) return 0.0;
    if(16.0*(c2 + 2.0*c.x + 1.0) - 1.0 < 0.0) return 0.0;
    
    const float B2 = 65536.0; // B*B pre-computed
    float n = 0.0;
    vec2 z = vec2(0.0);
    
    for(int i = 0; i < MAX_ITER; i++)
    {
        // Optimized complex multiplication
        float zx2 = z.x * z.x;
        float zy2 = z.y * z.y;
        
        if(zx2 + zy2 > B2) break;
        
        z = vec2(zx2 - zy2, 2.0 * z.x * z.y) + c;
        n += 1.0;
    }
    
    if(n >= float(MAX_ITER)) return 0.0;
    
    // Smooth iteration count
    float sn = n - log2(log2(dot(z, z))) + 4.0;
    return sn;
}

void mainImage(out vec4 fragColor, in vec2 fragCoord)
{
    vec3 col = vec3(0.0);
    
#if AA>1
    for(int m = 0; m < AA; m++)
    for(int n = 0; n < AA; n++)
    {
        vec2 p = (-iResolution.xy + 2.0*(fragCoord.xy + vec2(float(m), float(n))/float(AA))) / iResolution.y;
        float w = float(AA*m + n);
        float time = iTime + 0.5 * (1.0/24.0) * w / float(AA*AA);
#else
        vec2 p = (-iResolution.xy + 2.0*fragCoord.xy) / iResolution.y;
        float time = iTime;
#endif
        
        float zoo = 0.62 + 0.38 * cos(0.07 * time);
        float angle = 0.05 * time;
        float coa = cos(angle);
        float sia = sin(angle);
        zoo = pow(zoo, 8.0);
        
        vec2 xy = vec2(p.x*coa - p.y*sia, p.x*sia + p.y*coa);
        vec2 c = vec2(-0.745, 0.186) + xy * zoo;

        float l = mandelbrot(c);
        
        // Base colors
        //vec3 bgPurple = vec3(0.2, 0.0, 0.4);
        //vec3 fgPink  = vec3(1.0, 0.4, 0.8);
        vec3 bgPurple = vec3(0.1, 0.0, 0.2);
        vec3 fgPink  = vec3(0.8, 0.2, 0.6);
        
        // Smooth interpolation for escaped points
        vec3 pinkInterp = bgPurple + (fgPink - bgPurple) * (0.5 + 0.5 * cos(l * 0.15));
        
        // Assign color: purple for Mandelbrot set, interpolated pink for edges
        col += (l < 0.5) ? bgPurple : pinkInterp;
        
#if AA>1
    }
    col /= float(AA*AA);
#endif
    
    fragColor = vec4(col, 1.0);
}
