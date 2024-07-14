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
uniform vec3 Vertices[16];

in vec3 wp;        // World position
in vec3 n;         // Normal
in vec4 c;         // Color
out vec4 fragColor;

const int MAX_VERTEXCOUNT = 16;
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
    return 500.0;
}

// Compute solid angle for a planar triangle as seen from the origin
float computeSolidAngle_Norm(vec3 va, vec3 vb, vec3 vc) {
    float numerator = abs(dot(va, cross(vb, vc)));
    float denom = 1.0 + dot(va, vb) + dot(va, vc) + dot(vb, vc);
    float halfSA = atan(numerator, denom);
    return 2.0 * (halfSA >= 0.0 ? halfSA : halfSA + 3.14159265359);
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

    float smallestDist = 10000000.0;
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

            float dist = length(projectedPoint - p);
            if (dist < smallestDist) {
                result.point = projectedPoint;
                result.caseType = (t > 0.001 && t < 0.999) ? 1 : 2;
                result.index = (t > 0.999) ? i : (i > 0 ? i - 1 : polygonVertexCount - 1);
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

    vec3 T1 = normalize(o - n * dotNV);
    vec3 T2 = normalize(cross(n, T1));
    mat3 w2t = mat3(T1, T2, n);

    // Shift shading point away from polygon plane if within epsilon to avoid numerical issues
    vec3 l = Vertices[0] - P;
    float height = dot(PolygonNormal, l);
    float planeEps = 1e-3;
    if (abs(height) < planeEps) {
        float shiftDist = planeEps - abs(height);
        vec3 shiftDir = (height < 0.0 ? 1.0 : -1.0) * PolygonNormal;
        vec3 shift = shiftDist * shiftDir;
        P += shift;
    }

    // Clip polygon by tangent plane
    vec3 clippedVa[MAX_VERTEXCOUNT_PLUS_ONE];
    int clippedVc = 0;

    // Initialize the clipped polygon vertices (simplified example)
    for (int i = 0; i < VertexCount; ++i) {
        clippedVa[i] = Vertices[i] - P;
    }
    clippedVc = VertexCount;

    vec3 color = vec3(0.0);
    if (clippedVc > 2) {
        vec3 lightPlaneN = normalize(w2t * PolygonNormal);

        // Find closest point limited to upper hemisphere
        float t = dot(clippedVa[0], lightPlaneN);
        vec3 closestPoint = t * lightPlaneN;

        // Clamp closest point to clipped polygon
        bool ccw = t > 0.0;
        ClosestPoint closest = clampPointToPolygon(clippedVa, clippedVc, ccw, closestPoint);
        vec3 closestPointDir = normalize(closest.point);

        // Initialize triangle count
        int tc = clippedVc - int(closest.caseType);

        // Initialize fixed vertex data of triangle fan v0
        vec3 v0 = (closest.caseType != 2) ? closestPointDir : normalize(clippedVa[closest.index]);

        float v0out = abs(dot(PolygonNormal, v0));
        float v0Le = getRadiance_World(v0) / v0out;

        // Initialize 2nd vertex of first triangle v1
        int i1 = (closest.index + 1) % clippedVc;
        vec3 v1 = normalize(clippedVa[i1]);
        float v1out = abs(dot(PolygonNormal, v1));
        float v1Le = getRadiance_World(v1) / v1out;

        float denom = 0.0;
        float Ld = 0.0;
        
        for (int i = 1; i <= tc; ++i) {
            int i2 = (closest.index + i + 1) % clippedVc;
            vec3 v2 = normalize(clippedVa[i2]);
            float v2out = abs(dot(PolygonNormal, v2));
            float v2Le = getRadiance_World(v2) / v2out;

            float sphEx = computeSolidAngle_Norm(v0, v1, v2);
            
            Ld += sphEx * (v0Le * -v0.x + v1Le * -v1.x + v2Le * -v2.x) / 3.0;
            denom += sphEx * (-v0.x + -v1.x + -v2.x);
            
            v1 = v2;
            v1out = v2out;
            v1Le = v2Le;
        }
        
        if (Ld > 0.0) {
            vec3 brdf = c.xyz / 3.14159265359;
            color += Ld / PolygonArea * brdf;
        }
    }

    fragColor = vec4(color, c.w);
}
