﻿#version 330 core

struct Vertex {
    vec3 wp;    // World position
    vec3 n;     // Normal
    vec4 c;     // Color
};

uniform vec3 CameraLocation;
uniform vec3 PolygonNormal;
uniform float PolygonArea;
uniform int VertexCount;
uniform vec3 Vertices[8];
uniform float ldtdc;
uniform float ldtdg;
uniform float IntensityMulti;
uniform float maxLDTValue;

uniform sampler2D LDTLUT; 
uniform float LUT_SIZE_X;
uniform float LUT_SIZE_Y;
float LUT_PIXEL_X = 1.0/LUT_SIZE_X;
float LUT_PIXEL_Y = 1.0/LUT_SIZE_Y;
float LUT_SCALE_X = (LUT_SIZE_X - 1.0)/LUT_SIZE_X;
float LUT_SCALE_Y = (LUT_SIZE_Y - 1.0)/LUT_SIZE_Y;
float LUT_BIAS_X  = 0.5/LUT_SIZE_X;
float LUT_BIAS_Y  = 0.5/LUT_SIZE_Y;

in vec3 wp;        // World position
in vec3 n;         // Normal
in vec4 c;         // Color

in vec3 worldPosition;
in vec3 worldNormal;
in vec2 texcoord;

struct Light
{
    vec3 color;
    vec3 points[4];
    bool twoSided;
};
uniform Light areaLight;
uniform vec3 areaLightTranslate;

struct Material
{
    sampler2D diffuse;
    vec4 albedoRoughness; // (x,y,z) = color, w = roughness
};
uniform Material material;

uniform vec3 viewPosition;
uniform sampler2D LTC1; // for inverse M
uniform sampler2D LTC2; // GGX norm, fresnel, 0(unused), sphere

const float LUT_SIZE  = 64.0; // ltc_texture size
const float LUT_SCALE = (LUT_SIZE - 1.0)/LUT_SIZE;
const float LUT_BIAS  = 0.5/LUT_SIZE;

// Vector form without project to the plane (dot with the normal)
// Use for proxy sphere clipping
vec3 IntegrateEdgeVec(vec3 v1, vec3 v2)
{
    // Using built-in acos() function will result flaws
    // Using fitting result for calculating acos()
    float x = dot(v1, v2);
    float y = abs(x);

    float a = 0.8543985 + (0.4965155 + 0.0145206*y)*y;
    float b = 3.4175940 + (4.1616724 + y)*y;
    float v = a / b;

    float theta_sintheta = (x > 0.0) ? v : 0.5*inversesqrt(max(1.0 - x*x, 1e-7)) - v;

    return cross(v1, v2)*theta_sintheta;
}

float IntegrateEdge(vec3 v1, vec3 v2)
{
    return IntegrateEdgeVec(v1, v2).z;
}

// P is fragPos in world space (LTC distribution)
vec3 LTC_Evaluate(vec3 N, vec3 V, vec3 P, mat3 Minv, vec3 points[4], bool twoSided)
{
    // construct orthonormal basis around N
    vec3 T1, T2;
    T1 = normalize(V - N * dot(V, N));
    T2 = cross(N, T1);

    // rotate area light in (T1, T2, N) basis
    Minv = Minv * transpose(mat3(T1, T2, N));

    // polygon (allocate 4 vertices for clipping)
    vec3 L[4];
    // transform polygon from LTC back to origin Do (cosine weighted)
    L[0] = Minv * (points[0] - P);
    L[1] = Minv * (points[1] - P);
    L[2] = Minv * (points[2] - P);
    L[3] = Minv * (points[3] - P);

    // use tabulated horizon-clipped sphere
    // check if the shading point is behind the light
    vec3 dir = points[0] - P; // LTC space
    vec3 lightNormal = cross(points[1] - points[0], points[3] - points[0]);
    bool behind = (dot(dir, lightNormal) < 0.0);

    // cos weighted space
    L[0] = normalize(L[0]);
    L[1] = normalize(L[1]);
    L[2] = normalize(L[2]);
    L[3] = normalize(L[3]);

    // integrate
    vec3 vsum = vec3(0.0);
    vsum += IntegrateEdgeVec(L[0], L[1]);
    vsum += IntegrateEdgeVec(L[1], L[2]);
    vsum += IntegrateEdgeVec(L[2], L[3]);
    vsum += IntegrateEdgeVec(L[3], L[0]);

    // form factor of the polygon in direction vsum
    float len = length(vsum);

    float z = vsum.z/len;
    if (behind)
        z = -z;

    vec2 uv = vec2(z*0.5f + 0.5f, len); // range [0, 1]
    uv = uv*LUT_SCALE + LUT_BIAS;

    // Fetch the form factor for horizon clipping
    float scale = texture(LTC2, uv).w;

    float sum = len*scale;
    if (!behind && !twoSided)
        sum = 0.0;

    // Outgoing radiance (solid angle) for the entire polygon
    vec3 Lo_i = vec3(sum, sum, sum);
    return Lo_i;
}

// PBR-maps for roughness (and metallic) are usually stored in non-linear
// color space (sRGB), so we use these functions to convert into linear RGB.
vec3 PowVec3(vec3 v, float p)
{
    return vec3(pow(v.x, p), pow(v.y, p), pow(v.z, p));
}

const float gamma = 2.2;
vec3 ToLinear(vec3 v) { return PowVec3(v, gamma); }
vec3 ToSRGB(vec3 v)   { return PowVec3(v, 1.0/gamma); }

vec3 getLTCSpec()
{
    // gamma correction
    vec3 mDiffuse = vec3(0.7f, 0.8f, 0.96f);// * texture(material.diffuse, texcoord).xyz;
    vec3 mSpecular = ToLinear(vec3(1.f, 1.f, 1.f)); // mDiffuse

    vec3 result = vec3(0.0f);

    vec3 N = normalize(worldNormal);
    vec3 V = normalize(viewPosition - worldPosition);
    vec3 P = worldPosition;
    float dotNV = clamp(dot(N, V), 0.0f, 1.0f);

    // use roughness and sqrt(1-cos_theta) to sample M_texture
    vec2 uv = vec2(material.albedoRoughness.w, sqrt(1.0f - dotNV));
    uv = uv*LUT_SCALE + LUT_BIAS;

    // get 4 parameters for inverse_M
    vec4 t1 = texture(LTC1, uv);

    // Get 2 parameters for Fresnel calculation
    vec4 t2 = texture(LTC2, uv);

    mat3 Minv = mat3(
        vec3(t1.x, 0, t1.y),
        vec3(  0,  1,    0),
        vec3(t1.z, 0, t1.w)
    );

    // translate light source for testing
    vec3 translatedPoints[4];
    translatedPoints[0] = areaLight.points[0] + areaLightTranslate;
    translatedPoints[1] = areaLight.points[1] + areaLightTranslate;
    translatedPoints[2] = areaLight.points[2] + areaLightTranslate;
    translatedPoints[3] = areaLight.points[3] + areaLightTranslate;

    // Evaluate LTC shading
    vec3 diffuse = LTC_Evaluate(N, V, P, mat3(1), translatedPoints, areaLight.twoSided);
    vec3 specular = LTC_Evaluate(N, V, P, Minv, translatedPoints, areaLight.twoSided);

    // GGX BRDF shadowing and Fresnel
    // t2.x: shadowedF90 (F90 normally it should be 1.0)
    // t2.y: Smith function for Geometric Attenuation Term, it is dot(V or L, H).
    specular *= mSpecular*t2.x + (1.0f - mSpecular) * t2.y;

    // result = areaLight.color * areaLight.intensity * (specular + mDiffuse * diffuse);
    // result = areaLight.color * areaLight.intensity * mDiffuse * diffuse;
    result = areaLight.color * specular;
    //return ToSRGB(result);
    return result;
}

out vec4 fragColor;

const int MAX_VERTEXCOUNT = 8;
const int MAX_VERTEXCOUNT_PLUS_ONE = MAX_VERTEXCOUNT + 1;

vec3 normalize(vec3 v) {
    return v / length(v);
}
float dot(vec3 a, vec3 b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}
vec3 cross(vec3 a, vec3 b) {
    return vec3(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x);
}
vec3 mix(vec3 a, vec3 b, float t) {
    return a * (1.0 - t) + b * t;
}

float getRadiance_World(vec3 dir){
    const float M_PI = 3.14159265359;
    vec3 v = normalize(dir);
    float C = atan(-v.y, -v.z) + M_PI, gamma = M_PI - acos(v.x);
    
    int sz = int(180.0/ldtdg)+1;
    int Cindex = int(C/M_PI*180.0/ldtdc);
    int gammaindex = int(gamma/M_PI*180.0/ldtdg);
    //if(gamma==0.0 || gamma==M_PI)
      //  return intensityDis[Cindex*sz+gammaindex];
    float d = 0.0, e = 0.0;
    while(d+ldtdg<gamma/M_PI*180.0)
        d += ldtdg;
    while(e+ldtdc<C/M_PI*180.0)
        e += ldtdc;
    float a = 1.0-(gamma/M_PI*180.0-d)/ldtdg;
    float b = 1.0-(C/M_PI*180.0-e)/ldtdc;
    float value1 = (a*
    texelFetch(LDTLUT, ivec2(gammaindex, Cindex), 0).r*maxLDTValue
    +(1-a)*
    texelFetch(LDTLUT, ivec2(gammaindex+1, Cindex), 0).r*maxLDTValue
    );

    float value2 = (a*
    texelFetch(LDTLUT, ivec2(gammaindex, Cindex+1), 0).r*maxLDTValue
    +(1-a)*
    texelFetch(LDTLUT, ivec2(gammaindex+1, Cindex+1), 0).r*maxLDTValue
    );
    return 30 * (b*value1 + (1-b)*value2)/683;
}

// Compute solid angle for a planar triangle as seen from the origin
float computeSolidAngle_Norm(vec3 va, vec3 vb, vec3 vc) {
    float numerator = abs(dot(va, cross(vb, vc)));
    float denom = 1.0 + dot(va, vb) + dot(va, vc) + dot(vb, vc);
    float halfSA = atan(numerator, denom);
    return 2.0 * ((halfSA >= 0.0) ? halfSA : (halfSA + 3.14159265359));
}

struct ClosestPoint {
    vec3 point;
    int index;
    int caseType;
};

// Clamp point to polygon
ClosestPoint clampPointToPolygon(vec3 polygonVertices[MAX_VERTEXCOUNT_PLUS_ONE], int polygonVertexCount, bool polygonClockwiseOrder, vec3 p) {
    ClosestPoint result;
    result.caseType = 0;
    result.index = -1;
    result.point = p;
    float smallestDist = 100000.0;

    float flipSng = polygonClockwiseOrder ? -1.0 : 1.0;
    
    vec3 v0 = polygonVertices[polygonVertexCount - 1];
    for (int i = 0; i < polygonVertexCount; ++i) {
        vec3 v1 = polygonVertices[i];

        vec3 edgePlaneN = cross(v0, v1);
        float dotPlane = dot(edgePlaneN, p);

        if (flipSng * dotPlane > -1e-9) {
            vec3 ab = v1 - v0;
            vec3 ap = p - v0;
            float lenSq = dot(ab, ab);
            float t = (lenSq > 1e-5) ? dot(ab, ap) / lenSq : 0.0;
            vec3 projectedPoint = v0 + clamp(t, 0.0, 1.0) * ab;

            float dist = length(projectedPoint - p)*length(projectedPoint - p);
            if (dist < smallestDist) {
                result.point = projectedPoint;
                result.caseType = (t > 0.001 && t < 0.999) ? 1 : 2;
                result.index = (t > 0.999) ? i : ((i > 0) ? (i - 1) : (polygonVertexCount - 1));
                smallestDist = dist;

                if (result.caseType == 1) break;
            }
        }
        v0 = v1;
    }

    return result;
}

// Main shading procedure
void main() {
    vec3 P = wp;
    if(P.x <= 0.0001 && P.x >=-0.0001)
    {
        fragColor = vec4(0, 0, 0, 1.0);
        return ;
    }
        
    vec3 n = normalize(n);

    // Create orthonormal basis around N
    vec3 o = normalize(CameraLocation - P);
    float dotNV = dot(n, o);

    // Shift shading point away from polygon plane if within epsilon to avoid numerical issues
    vec3 l = Vertices[0] - P;
    float height = dot(PolygonNormal, l);
    float planeEps = 1e-5;
    if (abs(height) < planeEps) {
        float shiftDist = planeEps - abs(height);
        vec3 shiftDir = (height < 0.0 ? 1.0 : -1.0) * PolygonNormal;
        vec3 shift = shiftDist * shiftDir;
        P += shift;
    }

    vec3 clippedVa[MAX_VERTEXCOUNT_PLUS_ONE];
    int clippedVc = 0;

    const float eps = 1e-9;
    vec3 vb = Vertices[0] - wp;
    float hb = vb.x;
    if(hb < 0)
        hb = abs(hb);
    bool hbv = hb > eps;
    bool hbn = hb < -eps;

    if (hb >= -eps) 
        clippedVa[clippedVc++] = vb;

    vec3 v0 = vb;
    float h0 = hb;
    bool h0v = hbv;
    bool h0n = hbn;

    for (int vi = 1; vi < VertexCount; vi++) {
        vec3 v1 = Vertices[vi] - wp;
        float h1 = v1.x;
        if(h1 < 0)
            h1 = abs(hb);
        bool h1v = h1 > eps;
        bool h1n = h1 < -eps;

        if ((h0v && h1n) || (h0n && h1v)) {
            vec3 intersection = mix(v0, v1, h0 / (h0 - h1));
            clippedVa[clippedVc++] = intersection;
        }

        if (h1 >= -eps) 
           clippedVa[clippedVc++] = v1;

        v0 = v1;
        h0 = h1;
        h0v = h1v;
        h0n = h1n;
    }

    // Handle the last edge to vertices[0]
    if ((h0v && hbn) || (h0n && hbv)) {
        vec3 intersection = mix(v0, vb, h0 / (h0 - hb));
        clippedVa[clippedVc++] = intersection;
    }

    vec3 color = vec3(0.0);
    if (clippedVc > 2) {
        vec3 lightPlaneN = normalize(PolygonNormal);

        // Find closest point limited to upper hemisphere
        float t = dot(clippedVa[0], lightPlaneN);
        vec3 closestPoint = t * lightPlaneN;

        // Clamp closest point to clipped polygon
        bool ccw = (t > 0.0);
        ClosestPoint closest = clampPointToPolygon(clippedVa, clippedVc, ccw, closestPoint);
        vec3 closestPointDir = normalize(closest.point);

        // Initialize triangle count
        int tc = clippedVc - int(closest.caseType);

        // Initialize fixed vertex data of triangle fan v0
        vec3 v0 = (closest.caseType != 2) ? closestPointDir : normalize(clippedVa[closest.index]);

        float v0out = abs(dot(PolygonNormal, v0));
        float v0Le = getRadiance_World(P-v0) / v0out;

        // Initialize 2nd vertex of first triangle v1
        int i1 = (closest.index + 1) % clippedVc;
        vec3 v1 = normalize(clippedVa[i1]);
        float v1out = abs(dot(PolygonNormal, v1));
        float v1Le = getRadiance_World(P-v1) / v1out;

        float denom = 0.0;
        float Ld = 0.0;

        for (int i = 1; i <= tc; ++i) {
            int i2 = (closest.index + i + 1) % clippedVc;
            vec3 v2 = normalize(clippedVa[i2]);
            float v2out = abs(dot(PolygonNormal, v2));
            float v2Le = getRadiance_World(P-v2) / v2out;

            float sphEx = computeSolidAngle_Norm(v0, v1, v2);

            //Ld += sphEx * (v0Le  + v1Le  + v2Le ) / 3.0;
            // denom += sphEx * (-v0.x + -v1.x + -v2.x);

            float avgLe = (v0Le + v1Le + v2Le) / 3.0;
            float avgG = (v0.x + v1.x + v2.x) / 3.0;
            if(avgG < 0)
                avgG = abs(avgG);
            float G = sphEx * avgG;
            Ld += avgLe * G;
            denom += G;

            v1 = v2;
            v1out = v2out;
            v1Le = v2Le;
        }

        if (Ld > 0.0) {
            vec3 brdf = c.xyz / 3.14159265359;
            color += Ld / PolygonArea * brdf;

            if(denom > 0){
                float Le = Ld / denom;
                color += Le * getLTCSpec();
            }
        }
    }
    color.rgb *= IntensityMulti;
    color.rgb = pow(color.rgb, vec3(1.0 / 2.2));
    fragColor = vec4(color, 1.0);
}
