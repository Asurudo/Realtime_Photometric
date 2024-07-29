#version 330 core

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
uniform float intensityDis[950];
uniform float ldtdc;
uniform float ldtdg;
uniform float maxIntensity;

in vec3 wp;        // World position
in vec3 n;         // Normal
in vec4 c;         // Color
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
    float value1 = (a*intensityDis[Cindex*sz+gammaindex]+(1-a)*intensityDis[Cindex*sz+gammaindex+1]);
    float value2 = (a*intensityDis[(Cindex+1)*sz+gammaindex]+(1-a)*intensityDis[(Cindex+1)*sz+gammaindex+1]);
    return 6*(b*value1 + (1-b)*value2)/683;
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
            //denom += sphEx * (-v0.x + -v1.x + -v2.x);

            float avgLe = (v0Le + v1Le + v2Le) / 3.0;
            float avgG = (v0.x + v1.x + v2.x) / 3.0;
            if(avgG < 0)
                avgG = abs(avgG);
            float G = sphEx * avgG;
            Ld += avgLe * G;

            v1 = v2;
            v1out = v2out;
            v1Le = v2Le;
        }
        
        if (Ld > 0.0) {
            vec3 brdf = c.xyz / 3.14159265359;
            color += Ld / PolygonArea * brdf;
        }
    }
    color.rgb = pow(color.rgb, vec3(1.0 / 2.4));
    fragColor = vec4(color, 1.0);
}
