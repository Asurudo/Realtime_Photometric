#version 330 core

struct Vertex {
    vec3 wp;    // World position
    vec3 n;     // Normal
    vec4 c;     // Color
};

// uniform vec3 CameraLocation;
uniform vec3 PolygonNormal;
uniform float PolygonArea;
uniform int VertexCount;
uniform vec3 Vertices[8];
uniform float ldtdc;
uniform float ldtdg;
uniform float IntensityMulti;
uniform float maxLDTValue;
uniform float randNum1;
uniform float randNum2;

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
    float CosTheta = dot(v1, v2);
    float Theta = acos(CosTheta);
    return (cross(v1, v2) * ((Theta > 0.0001) ? Theta/sin(Theta) : 1.0)).z;
    
    return IntegrateEdgeVec(v1, v2).z;
}

// P is fragPos in world space (LTC distribution)
vec3 LTC_Evaluate(vec3 N, vec3 V, vec3 P, mat3 Minv, vec3 points[4], bool twoSided)
{
//    // construct orthonormal basis around N
    vec3 T1, T2;
    T1 = normalize(V - N * dot(V, N));
    T2 = cross(T1, N);

    // rotate area light in (T1, T2, N) basis
    Minv = Minv * transpose(mat3(T1, T2, N));
    //Minv = Minv * transpose(mat3(N, T2, T1));

    // polygon (allocate 4 vertices for clipping)
    vec3 L[5];
    // transform polygon from LTC back to origin Do (cosine weighted)
    L[0] = Minv * (points[0] - P);
    L[1] = Minv * (points[1] - P);
    L[2] = Minv * (points[2] - P);
    L[3] = Minv * (points[3] - P);
    L[4] = L[3];

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
    L[4] = normalize(L[4]);

    // integrate
    vec3 vsum = vec3(0.0);
    vsum += IntegrateEdgeVec(L[0], L[1]);
    vsum += IntegrateEdgeVec(L[1], L[2]);
    vsum += IntegrateEdgeVec(L[2], L[3]);
    vsum += IntegrateEdgeVec(L[3], L[4]);
    vsum += IntegrateEdgeVec(L[4], L[0]);

    // form factor of the polygon in direction vsum
    float len = length(vsum);

    float z = vsum.z/len;
    if (behind)
        z = -z;

//    vec2 uv = vec2(z*0.5f + 0.5f, len); // range [0, 1]
//    uv = uv*LUT_SCALE + LUT_BIAS;
//
//    // Fetch the form factor for horizon clipping
//    float scale = texture(LTC2, uv).x;
////
//    float sum = len*scale;
    if (!behind && !twoSided)
        len = 0.0;

    // Outgoing radiance (solid angle) for the entire polygon
    vec3 Lo_i = vec3(len, len, len);
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

float g(vec3 lwi, vec3 lwo)
{
	float tan_i = 1.0 / ( lwi.y * lwi.y ) - 1.0; 
	float tan_o = 1.0 / ( lwo.y * lwo.y ) - 1.0; 
	float lambda_i = ( - 1.0 + sqrt( 1.0 + material.albedoRoughness.w * material.albedoRoughness.w * tan_i ) ) / 2.0; 
	float lambda_o = ( - 1.0 + sqrt( 1.0 + material.albedoRoughness.w * material.albedoRoughness.w * tan_o ) ) / 2.0; 
	return 1.0 / ( 1.0 + lambda_i + lambda_o );
}

float d(vec3 lwi, vec3 lwo)
{
  vec3 h = normalize( lwi + lwo );
  float cos2 = h.y * h.y;
  float sin2 = max(0.0, 1.0 - cos2);
  float m_alpha =  material.albedoRoughness.w;
  float denom = 3.141592653 * m_alpha * m_alpha * ( cos2 + sin2 / ( m_alpha * m_alpha ) ) * ( cos2 + sin2 / ( m_alpha * m_alpha ) );
  return 1.0 / denom;

}

float f(vec3 lwi, vec3 lwo, float f0)
{
    vec3 h = normalize( lwi + lwo );
    float cosine = dot( h, lwi );
    float tmp = ( 1.0 - cosine ) * ( 1.0 - cosine ) * ( 1.0 - cosine ) * ( 1.0 - cosine ) * ( 1.0 - cosine );
    return f0 + ( 1.0 - f0 ) * tmp;
}

float eval(vec3 V, vec3 L, float alpha)
{
		if(V.y <= 0)
		{
			return 0;
		}

		// masking
		float a_V = 1.0f / alpha / tan(acos(V.y));
		float LambdaV = (V.y<1.0f) ? 0.5f * (-1.0f + sqrt(1.0f + 1.0f/a_V/a_V)) : 0.0f;
	    float G1 = 1.0f / (1.0f + LambdaV);

		// shadowing
		float G2;
		if(L.y <= 0.0f)
			G2 = 0;
		else
		{
			float a_L = 1.0f / alpha / tan(acos(L.y));
			float LambdaL = (L.y<1.0f) ? 0.5f * (-1.0f + sqrt(1.0f + 1.0f/a_L/a_L)) : 0.0f;
			G2 = 1.0f / (1.0f + LambdaV + LambdaL);
		}

		// D
		vec3 H = normalize(V+L);
		float slopex = H.x/H.y;
		float slopez = H.z/H.y;
		float D = 1.0f / (1.0f + (slopex*slopex+slopez*slopez)/alpha/alpha);
		D = D*D;
		D = D / (3.14159f * alpha * alpha * H.y*H.y*H.y*H.y);

		float res = D * G2 / 4.0f / V.y;
		return res;
}

float ReflectivityAdjust(float dotNV){
    return clamp((1.0-dotNV)+0.3,  0, 1)/2.4;
}

vec3 getLTCSpec()
{
    // gamma correction
    // vec3 mDiffuse = vec3(0.7f, 0.8f, 0.96f);// * texture(material.diffuse, texcoord).xyz;
     

    vec3 result = vec3(0.0f);

    vec3 N = normalize(worldNormal);
    vec3 V = normalize(viewPosition-worldPosition);
    //vec3 L = normalize(vec3(0, 3.5, 0) - worldPosition);
    vec3 P = worldPosition;

    float dotNV = clamp(dot(N, V), 0.0f, 1.0f);
    //float dotNL = clamp(dot(N, L), 0.0f, 1.0f);

    // use roughness and sqrt(1-cos_theta) to sample M_texture
    // vec2 uv = vec2(material.albedoRoughness.w, clamp(acos(dotNV)/3.141592653*2, 0, 1));
    vec2 uv = vec2(material.albedoRoughness.w,sqrt(1-dotNV));
    //return vec3(uv, 1.0);

    //vec2 uv = vec2(material.albedoRoughness.w, sqrt(1.0f-dotNV));
    uv = uv*LUT_SCALE + LUT_BIAS;
    
    // get 4 parameters for inverse_M
    vec4 t1 = texture(LTC1, uv);

    // Get 2 parameters for Fresnel calculation
    vec4 t2 = texture(LTC2, uv);
    
    float a = t1.x, b = t1.y, c = t1.z, d = t1.w;

//    mat3 Minv = mat3(
//        vec3(a*c/(a*a-a*b*d), 0, -b*c/(a-b*d)),
//        vec3(  0,  1,  0),
//        vec3(-d*c/(a-b*d), 0, a*c/(a-b*d))
//    );

    mat3 Minv = mat3(
        vec3( a, 0, c),
        vec3( 0, 1, 0),
        vec3( b, 0, d)
    );

//    mat3 Minv = mat3(
//        vec3( a, 0, d),
//        vec3( 0, c, 0),
//        vec3( b, 0, 1)
//    );
    //Minv = inverse(Minv);

    // translate light source for testing
    vec3 translatedPoints[4];
    translatedPoints[0] = areaLight.points[0] + areaLightTranslate;
    translatedPoints[1] = areaLight.points[1] + areaLightTranslate;
    translatedPoints[2] = areaLight.points[2] + areaLightTranslate;
    translatedPoints[3] = areaLight.points[3] + areaLightTranslate;

//     //construct orthonormal basis around N
//    vec3 T1, T2;
//    T1 = normalize(V - N * dot(V, N));
//    T1 = normalize(cross(N, (abs(N.x) < abs(N.z)) ? vec3(0.0, -N.z, N.x) : vec3(-N.y, N.x, 0.0)));
//    T2 = cross(N, T1);
//
//     //rotate area light in (T1, T2, N) basis
//    Minv = Minv * transpose(mat3(T1, T2, N));

    // Evaluate LTC shading
    vec3 diffuse = LTC_Evaluate(N, V, P, mat3(1), translatedPoints, areaLight.twoSided);
    vec3 specular = LTC_Evaluate(N, V, P, Minv, translatedPoints, areaLight.twoSided);
    specular *= t2.x;
    // GGX BRDF shadowing and Fresnel
    // t2.x: shadowedF90 (F90 normally it should be 1.0)
    // t2.y: Smith function for Geometric Attenuation Term, it is dot(V or L, H).
//    T1 = normalize(cross(N, (abs(N.x) < abs(N.z)) ? vec3(0.0, -N.z, N.x) : vec3(-N.y, N.x, 0.0)));
//    T2 = cross(N, T1);
//    L *= transpose(mat3(T1, N, T2));
//    V *= transpose(mat3(T1, N, T2));
    //specular *= 0.63*t2.x + (1.0f - 0.63) * t2.y;
    // float brdfValue = f(V,L,F0)*g(V,L)*d(V,L)/ (4.0*dotNV*dotNL);
    // vec3 brdf = vec3(brdfValue, brdfValue, brdfValue);
    // specular *= brdf*t2.x + (1.0f - brdf) * t2.y;
    // specular *= f(V,L,F0)*d(V,L)*g(V,L)/ (4.0*dotNV*dotNL);
     //specular *= t2.x + 0.9*t2.y;
    //specular *= f(V,L,F0)*d(V,L)*g(V,L)/ (4.0*dotNV*dotNL);
    //specular *= eval(V, L, material.albedoRoughness.w);
    // result = areaLight.color * areaLight.intensity * (specular + mDiffuse * diffuse);
    // result = areaLight.color * areaLight.intensity * mDiffuse * diffuse;
    //specular *= (1.0-ReflectivityAdjust(dotNV))*10;
    vec3 ks = vec3((1-material.albedoRoughness.w)/8, (1-material.albedoRoughness.w)/8, (1-material.albedoRoughness.w)/8); 
    result = areaLight.color * (ks*specular + diffuse*(1-ks)/2);
    //result = areaLight.color * (ks*specular);

    return result/3.141592653/2;
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
    return 50;
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

    float tex1 = texelFetch(LDTLUT, ivec2(gammaindex, Cindex), 0).r;
    float tex2 = texelFetch(LDTLUT, ivec2(gammaindex+1, Cindex), 0).r;

    float value1 = (a*tex1*maxLDTValue+(1-a)*tex2*maxLDTValue);

    float tex3 = texelFetch(LDTLUT, ivec2(gammaindex, Cindex+1), 0).r;
    float tex4 = texelFetch(LDTLUT, ivec2(gammaindex+1, Cindex+1), 0).r;

    float value2 = (a*tex3*maxLDTValue+(1-a)*tex4*maxLDTValue);
    return 600 * (b*value1 + (1-b)*value2)/683;
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
    float smallestDist = 1000.0;

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
    if(wp.x>0){
        fragColor = vec4(0,0,0,1.0);
        return ;
    }
        
    fragColor = vec4(getRadiance_World(vec3(0))*getLTCSpec().rgb,1.0);
    //fragColor = vec4(pow(getRadiance_World(vec3(0))*getLTCSpec().rgb, vec3(1.0 / 1.5)), 1.0);
    return ;

    vec3 P = wp;
    vec3 n = normalize(n);

    // Create orthonormal basis around N
    // vec3 o = normalize(viewPosition - P);
    //float dotNV = dot(n, o);

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
        float h1 = v1.y;
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

            // Ld += sphEx * (v0Le  + v1Le  + v2Le ) / 3.0;
            // denom += sphEx * (-v0.x + -v1.x + -v2.x);

            float avgLe = (v0Le + v1Le + v2Le) / 3.0;
            float avgG = (v0.y + v1.y + v2.y) / 3.0;
            if(avgG < 0)
                avgG = abs(avgG);
            float G = sphEx * avgG;
            Ld += avgLe * G;
            denom += G;

            v1 = v2;
            v1out = v2out;
            v1Le = v2Le;
        }
        vec3 N = normalize(worldNormal);
        vec3 V = normalize(viewPosition-worldPosition);
        if (Ld > 0.0) {
            vec3 brdf = c.xyz / 3.14159265359;
            // the diffuse
            vec3 ks = vec3((1-material.albedoRoughness.w)/8, (1-material.albedoRoughness.w)/8, (1-material.albedoRoughness.w)/8); 
             color += Ld / PolygonArea * (1-ks) * brdf * ReflectivityAdjust(dot(N ,V));

            if(denom > 0){
                float Le = Ld / denom;
                // the specular
                color += Le * getLTCSpec() * ReflectivityAdjust(dot(N ,V));
            }
        }
    }
    color.rgb *= IntensityMulti;
    //color.rgb /= 2.4;
    color.rgb = pow(color.rgb, vec3(1.0 / 1.5));
    fragColor = vec4(color, 1.0);
}
