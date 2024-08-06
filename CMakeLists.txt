cmake_minimum_required(VERSION 3.8)
project ("Realtime_Photometric")

# 包括头文件目录
include_directories(${PROJECT_SOURCE_DIR}/include ${PROJECT_SOURCE_DIR}/src)

# 设置库目录
link_directories(${PROJECT_SOURCE_DIR}/include/lib)

add_library(imgui STATIC 
    ${PROJECT_SOURCE_DIR}/include/imgui/imgui.cpp 
    ${PROJECT_SOURCE_DIR}/include/imgui/imgui_impl_glfw.cpp 
    ${PROJECT_SOURCE_DIR}/include/imgui/imgui_impl_opengl3.cpp
    ${PROJECT_SOURCE_DIR}/include/imgui/imgui_draw.cpp 
    ${PROJECT_SOURCE_DIR}/include/imgui/imgui_tables.cpp 
    ${PROJECT_SOURCE_DIR}/include/imgui/imgui_widgets.cpp
    ${PROJECT_SOURCE_DIR}/include/imgui/imgui_demo.cpp
 "src/App.h")

link_libraries(imgui)

# 添加源代码到此项目的可执行文件
add_executable(${PROJECT_NAME} "Realtime_Photometric.cpp" "Realtime_Photometric.h" "include/glad/gl.h"   "include/GLFW/glfw3.h" "include/GLFW/glfw3native.h" "include/glm/detail/glm.cpp" "include/glm/simd/common.h" "include/glm/simd/exponential.h" "include/glm/simd/geometric.h" "include/glm/simd/integer.h" "include/glm/simd/matrix.h" "include/glm/simd/neon.h" "include/glm/simd/packing.h" "include/glm/simd/platform.h" "include/glm/simd/trigonometric.h" "include/glm/simd/vector_relational.h" "include/KHR/khrplatform.h" "src/shader.h" "include/assimp/config.h" "include/assimp/revision.h" "include/assimp/aabb.h" "include/assimp/ai_assert.h" "include/assimp/anim.h" "include/assimp/AssertHandler.h" "include/assimp/BaseImporter.h" "include/assimp/Bitmap.h" "include/assimp/BlobIOSystem.h" "include/assimp/ByteSwapper.h" "include/assimp/camera.h" "include/assimp/cexport.h" "include/assimp/cfileio.h" "include/assimp/cimport.h" "include/assimp/ColladaMetaData.h" "include/assimp/color4.h" "include/assimp/commonMetaData.h" "include/assimp/CreateAnimMesh.h" "include/assimp/DefaultIOStream.h" "include/assimp/DefaultIOSystem.h" "include/assimp/defs.h" "include/assimp/Exceptional.h" "include/assimp/fast_atof.h" "include/assimp/GenericProperty.h" "include/assimp/GltfMaterial.h" "include/assimp/Hash.h" "include/assimp/importerdesc.h" "include/assimp/IOStreamBuffer.h" "include/assimp/light.h" "include/assimp/LineSplitter.h" "include/assimp/LogAux.h" "include/assimp/material.h" "include/assimp/MathFunctions.h" "include/assimp/matrix3x3.h" "include/assimp/matrix4x4.h" "include/assimp/MemoryIOWrapper.h" "include/assimp/mesh.h" "include/assimp/metadata.h" "include/assimp/ObjMaterial.h" "include/assimp/ParsingUtils.h" "include/assimp/pbrmaterial.h" "include/assimp/postprocess.h" "include/assimp/Profiler.h" "include/assimp/qnan.h" "include/assimp/quaternion.h" "include/assimp/RemoveComments.h" "include/assimp/scene.h" "include/assimp/SceneCombiner.h" "include/assimp/SGSpatialSort.h" "include/assimp/SkeletonMeshBuilder.h" "include/assimp/SmallVector.h" "include/assimp/SmoothingGroups.h" "include/assimp/SpatialSort.h" "include/assimp/StandardShapes.h" "include/assimp/StreamReader.h" "include/assimp/StreamWriter.h" "include/assimp/StringComparison.h" "include/assimp/StringUtils.h" "include/assimp/Subdivision.h" "include/assimp/texture.h" "include/assimp/TinyFormatter.h" "include/assimp/types.h" "include/assimp/vector2.h" "include/assimp/vector3.h" "include/assimp/version.h" "include/assimp/Vertex.h" "include/assimp/XmlParser.h" "include/assimp/XMLTools.h" "include/assimp/ZipArchiveIOSystem.h" "include/assimp/Compiler/poppack1.h" "include/assimp/Compiler/pstdint.h" "include/assimp/Compiler/pushpack1.h" "include/assimp/port/AndroidJNI/AndroidJNIIOSystem.h" "include/assimp/port/AndroidJNI/BundledAssetIOSystem.h" "include/assimp/config.h" "include/assimp/revision.h" "src/App.h")

# 链接库
target_link_libraries(${PROJECT_NAME} PRIVATE glfw3 ${OPENGL_LIBRARY} assimp-vc143-mt)

# 设置 C++20 标准
if (CMAKE_VERSION VERSION_GREATER 3.12)
  set_property(TARGET ${PROJECT_NAME} PROPERTY CXX_STANDARD 20)
endif()
