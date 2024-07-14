#version 330 core

in vec3 wp;
in vec3 n;
in vec4 c;

out vec4 fragColor;

uniform vec3 CameraLocation;
uniform vec3 PolygonNormal;
// 区域光顶点个数与顶点坐标
uniform vec3 Vertices[16];
uniform int VertexCount;
// 矩形区域面积
uniform float PolygonArea;

vec3 normalize(vec3 v) {
    return v / length(v);
}
float dot(vec3 a, vec3 b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}
vec3 cross(vec3 a, vec3 b) {
    return vec3(a.y * b.z - a.z * b.y,
                a.z * b.x - a.x * b.z,
                a.x * b.y - a.y * b.x);
}
vec3 mix(vec3 a, vec3 b, float t) {
    return a * (1.0 - t) + b * t;
}

// 根据世界坐标的方向得到对应的光强值
float getRadiance_World(vec3 dir){
    return 500.f;
}
// 计算大Ω
float computeSolidAngle_Norm(vec3 va, vec3 vb, vec3 vc) {
    float numerator = abs(dot(va, cross(vb, vc)));
    float denom = 1.0 + dot(va, vb) + dot(va, vc) + dot(vb, vc);
    float halfSA = atan(numerator, denom);
    return 2.0 * (halfSA >= 0.0 ? halfSA : halfSA + 3.14159265359); 
}

//float evalLTCSpec(vec3 P, mat3 w2t, float roughness, float dotNV, vec3 Vertices[], int VertexCount);

void main() {
    const float PI = 3.14159265359;
    
    vec3 P = wp;
    vec3 n = normalize(n);

    // create orthonormal basis around N
    vec3 o = normalize(CameraLocation - P);
    float dotNV = dot(n, o);
    
    vec3 T1 = normalize(o - n * dotNV);
    vec3 T2 = normalize(cross(n, T1));
    // 世界空间转换到切线空间的三维矩阵
    mat3 w2t = mat3(T1, T2, n);



    // 为了防止一些数值问题做一些调整
    // shift shading point away from polygon plane if within epsilon to avoid numerical issues
    vec3 l = Vertices[0] - P;
    float height = dot(PolygonNormal, l);
    float planeEps = 1e-3;
    if (abs(height) < planeEps) {
        float shiftDist = planeEps - abs(height);
        vec3 shiftDir = (height < 0.0 ? 1.0 : -1.0) * PolygonNormal;
        vec3 shift = shiftDist * shiftDir;
        P += shift;
    }

    // 区域光法线转换为切线坐标系
    // polygon normal in tangent space
    vec3 lightPlaneN = normalize(w2t * PolygonNormal);

    // 切线坐标系下shading point到第一个顶点
    // first vertex of polygon in tangent space
    vec3 v0 = w2t * (Vertices[0] - P);

    // 寻找切线坐标系下最近点CP
    // find closest on polygon plane point limited to upper hemisphere
    float t = dot(v0, lightPlaneN);
    vec3 closestPoint = t * lightPlaneN;

    // perform all initialization steps in one loop
    //  - polygon transformation to tangent space  将区域光的每个顶点转换到切平面
    //  - clip by z=0 不可见部分裁剪
    //  - clamp closest point 看看CP点是不是在边界上
    //  - project to hemisphere 投影到半球面
    //  - initialize radiance 计算放射度

    // 最近点在区域内还是边界还是顶点处
    int pointcase = 0; // Inside
    vec3 closestPointClamped = closestPoint;
    int closestIndex = -1; // inside
    float smallestDist = 10000000.0;
    float flipSng = (t > 0.0) ? -1.0 : 1.0;

    float eps = 1e-9;
    int vc = 0;
    vec4 va[16];
    
    vec3 vb = w2t * (Vertices[0] - P);
    
    vec3 v00;
    vec3 vl;


    if (vb.z >= -eps) {
        v00 = vb;
        vl = vb;
        vec3 dir = normalize(vb);
        va[0] = vec4(dir, getRadiance_World(-transpose(w2t) * dir));
        vc = 1;
    }

    v0 = vb;
    bool h0v = vb.z > eps;
    bool h0n = vb.z < -eps;
    
    // 遍历每个顶点
    for (int vi = 1; vi < VertexCount; ++vi) {
        vec3 v1 = w2t * (Vertices[vi] - P);
        bool h1v = v1.z > eps;
        bool h1n = v1.z < -eps;
        if ((h0v && h1n) || (h0n && h1v)) 
        {
            vec3 ve = mix(v0, v1, v0.z / (v0.z - v1.z));

            // clamp closest to new edge
            if (vc > 0) {
                vec3 edgePlaneN = cross(vl, ve);
                float dotPlane = dot(edgePlaneN, closestPoint);

                // check if point is outside the polygon
                if (flipSng * dotPlane > -1e-9) {
                    vec3 ab = ve - vl;
                    vec3 ap = closestPoint - vl;
                    float lenSq = dot(ab, ab);
                    float t = (lenSq > 1e-5) ? dot(ab, ap) / lenSq : 0.0;

                    vec3 projectedPoint = vl + clamp(t, 0.0, 1.0) * ab;

                    // check for projected point distance -> take closest
                    float dist = length(projectedPoint - closestPoint);
                    if (dist < smallestDist) {
                        closestPointClamped = projectedPoint;
                        pointcase = (t > 0.001 && t < 0.999) ? 1 : 2; // Edge or Vertex
                        closestIndex = (t > 0.999) ? vc : vc - 1;
                        smallestDist = dist;
                    }
                } else {
                    v00 = ve;
                }
            }

            // set next vertex
            vl = ve;
            vec3 dir = normalize(ve);
            va[vc] = vec4(dir, getRadiance_World(-transpose(w2t) * dir));
            vc++;
        }

        if (v1.z >= -eps) {
            // clamp closest point to new edge
            if (vc > 0) {
                vec3 edgePlaneN = cross(vl, v1);
                float dotPlane = dot(edgePlaneN, closestPoint);

                // check if point is outside the polygon
                if (flipSng * dotPlane > -1e-9) {
                    vec3 ab = v1 - vl;
                    vec3 ap = closestPoint - vl;
                    float lenSq = dot(ab, ab);
                    float t = (lenSq > 1e-5) ? dot(ab, ap) / lenSq : 0.0;

                    vec3 projectedPoint = vl + clamp(t, 0.0, 1.0) * ab;

                    // check for projected point distance -> take closest
                    float dist = length(projectedPoint - closestPoint);
                    if (dist < smallestDist) {
                        closestPointClamped = projectedPoint;
                        pointcase = (t > 0.001 && t < 0.999) ? 1 : 2; // Edge or Vertex
                        closestIndex = (t > 0.999) ? vc : vc - 1;
                        smallestDist = dist;
                    }
                } else {
                    v00 = v1;
                }
            }

            // set next vertex
            vl = v1;
            vec3 dir = normalize(v1);
            va[vc] = vec4(dir, getRadiance_World(-transpose(w2t) * dir));
            vc++;
        }

        v0 = v1;
        h0v = h1v;
        h0n = h1n;
    }

    // 最后一个顶点和第一个顶点连成的边
    // last edge to vertices[0]
    bool hbv = vb.z > eps;
    bool hbn = vb.z < -eps;
    if ((h0v && hbn) || (h0n && hbv)) {
        vec3 v1 = mix(v0, vb, v0.z / (v0.z - vb.z));
        // clamp closest point to new edge
        if (vc > 0) {
            vec3 edgePlaneN = cross(vl, v1);
            float dotPlane = dot(edgePlaneN, closestPoint);

            // check if point is outside the polygon
            if (flipSng * dotPlane > -1e-9) {
                vec3 ab = v1 - vl;
                vec3 ap = closestPoint - vl;
                float lenSq = dot(ab, ab);
                float t = (lenSq > 1e-5) ? dot(ab, ap) / lenSq : 0.0;

                vec3 projectedPoint = vl + clamp(t, 0.0, 1.0) * ab;

                // check for projected point distance -> take closest
                float dist = length(projectedPoint - closestPoint);
                if (dist < smallestDist) {
                    closestPointClamped = projectedPoint;
                    pointcase = (t > 0.001 && t < 0.999) ? 1 : 2; // Edge or Vertex
                    closestIndex = (t > 0.999) ? 0 : vc - 1;
                    smallestDist = dist;
                }
            } else {
                v00 = v1;
            }
        }

        // set next vertex
        vl = v1;
        vec3 dir = normalize(v1);
        va[vc] = vec4(dir, getRadiance_World(-transpose(w2t) * dir));
        vc++;
    }

    vec3 color = vec3(0.0);
    if (vc > 2) {
        // clamp closest point to last edge
        vec3 edgePlaneN = cross(vl, v00);
        float dotPlane = dot(edgePlaneN, closestPoint);

        // check if point is outside the polygon
        if (flipSng * dotPlane > -1e-9) {
            vec3 ab = v00 - vl;
            vec3 ap = closestPoint - vl;
            float lenSq = dot(ab, ab);
            float t = (lenSq > 1e-5) ? dot(ab, ap) / lenSq : 0.0;

            vec3 projectedPoint = vl + clamp(t, 0.0, 1.0) * ab;

            // check for projected point distance -> take closest
            float dist = length(projectedPoint - closestPoint);
            if (dist < smallestDist) {
                closestPointClamped = projectedPoint;
                pointcase = (t > 0.001 && t < 0.999) ? 1 : 2; // Edge or Vertex
                closestIndex = (t > 0.999) ? 0 : vc - 1;
                smallestDist = dist;
            }
        }

        vec3 closestPointDir = normalize(closestPointClamped);

        // 计算三角形个数
        // init triangle count: VertexCount in case of closest point is inside polygon
        //                      VertexCount-1 triangles in case of edge
        //                      VertexCount-2 triangles in case of corner
        int tc = vc - pointcase;
        
        // init fixed vertex data of triangle fan
        vec4 v0;
        if (pointcase != 2) {
            vec3 dir = closestPointDir;
            vec3 iw = -transpose(w2t) * dir;
            float Le = getRadiance_World(iw); // note: includes 1/dotOut
            v0 = vec4(closestPointDir, Le);
        } else {
            v0 = va[closestIndex];
        }
        
        float denom = 0.0;
        float Ld = 0.0;
        
        for (int i = 1; i <= tc; ++i) {
            int i1 = (closestIndex + i) % vc;
            int i2 = (closestIndex + i + 1) % vc;
                    
            float sphEx = computeSolidAngle_Norm(v0.xyz, va[i1].xyz, va[i2].xyz);
                  
            float avgLe = (v0.w + va[i1].w + va[i2].w) / 3.0;
            float avgG = (v0.z + va[i1].z + va[i2].z) / 3.0;
            float G = sphEx * avgG;
            Ld += avgLe * G;
            denom += G;
    
        }

        if (Ld > 0.0) {
            // diffuse shading
            vec3 brdf = c.xyz * (1.0 / PI);
            color += Ld / PolygonArea * brdf;

            // specular
            // if (ltcSpecular && denom > 0.0) {
            //    float Le = Ld / denom;

            //    vec3 ks = vec3(1.0);
            //    float roughness = 0.1; // "linear" roughness instead of trowbridge-reitz parameter ?? -> assume so, as Point shader matches then
                                                    
            //    float ltcSpec = evalLTCSpec(P, w2t, roughness, dotNV, Vertices, VertexCount);
                                    
            //    color += ks * (ltcSpec * Le);
            //}
        }
    }

    fragColor = vec4(color, c.w);
}
