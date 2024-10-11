#include<iostream>
#include<cstdio>
#include<string>
#include <random>
#include<algorithm>
#include<cmath>
#include <filesystem>

#define GLAD_GL_IMPLEMENTATION
#include <glad/gl.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include "imgui/imgui.h"
#include "imgui/imgui_impl_glfw.h"
#include "imgui/imgui_impl_opengl3.h"

#include "shader.h"
#include "App.h"
#include "camera.h"
#include "model.h"
#include "mesh.h"
#include "colors.hpp"
#include "ltc_matrix.hpp"
#include "tiny_ldt.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

const float M_PI = 3.141592653;

bool screenShotFlag = true;

// 屏幕大小
const unsigned int SCR_WIDTH = 800;
const unsigned int SCR_HEIGHT = 600;

// 光源
const glm::vec3 LIGHT_COLOR = Color::White;
// 光源乘数
float IntensityMulti = 1.0;
std::string lightType = "LINETIK-S_42184482";
// 粗糙度
static float roughness = 0.6f;

// std::string lightType = "PERLUCE_42182932";
glm::vec3 areaLightTranslate;
Shader* ltcShaderPtr;

// 控制
bool keys[1024]; // activated keys

// 摄像机
// Camera camera(glm::vec3(-28.8, 7.4f, 12.0), glm::vec3(0.0f, 1.0f, 0.0f), -39.3f, -22.6f);
Camera camera(glm::vec3(-25, 2.f, 0), glm::vec3(0.0f, 1.0f, 0.0f), 0.0f, 0.0f);
// Camera camera(glm::vec3(0, 60.f, 0), glm::vec3(0.0f, 1.0f, 0.0f), 0.0f,-90.0f);
float lastX = (float)SCR_WIDTH / 2.0f;
float lastY = (float)SCR_HEIGHT / 2.0f;
bool firstMouse = true;

// 时间
float deltaTime = 0.0f;	// 这一帧与上一帧之间的间隔时间
float lastFrame = 0.0f;

// 光度学文件
tiny_ldt<float>::light ldt;
float maxLDTValue = -1.0;
std::vector<std::vector<float>> intensityDis;

// 面光源与平面定义
// 2---3-5
// |  / /|
// | / / |
// |/ /  |
// 1-4---6
//
GLuint planeVBO, planeVAO;
GLuint areaLightVBO, areaLightVAO;
struct VertexAL {
	glm::vec3 position;
	glm::vec3 normal;
	glm::vec2 texcoord;
};

const GLfloat psize = 40.0f;
VertexAL planeVertices[6] = {
	{ {-psize, 0.0f, -psize}, {0.0f, 1.0f, 0.0f}, {0.0f, 0.0f} },
	{ {-psize, 0.0f,  psize}, {0.0f, 1.0f, 0.0f}, {0.0f, 1.0f} },
	{ { psize, 0.0f,  psize}, {0.0f, 1.0f, 0.0f}, {1.0f, 1.0f} },
	{ {-psize, 0.0f, -psize}, {0.0f, 1.0f, 0.0f}, {0.0f, 0.0f} },
	{ { psize, 0.0f,  psize}, {0.0f, 1.0f, 0.0f}, {1.0f, 1.0f} },
	{ { psize, 0.0f, -psize}, {0.0f, 1.0f, 0.0f}, {1.0f, 0.0f} }
};
VertexAL areaLightVertices[6] = {
	{ {0.0f, 3.4f, -1.5f}, {1.0f, 0.0f, 0.0f}, {0.0f, 0.0f} }, // 0 1 5 4
	{ {0.0f, 3.4f,  1.5f}, {1.0f, 0.0f, 0.0f}, {0.0f, 1.0f} },
	{ {0.0f, 0.4f,  1.5f}, {1.0f, 0.0f, 0.0f}, {1.0f, 1.0f} },
	{ {0.0f, 3.4f, -1.5f}, {1.0f, 0.0f, 0.0f}, {0.0f, 0.0f} },
	{ {0.0f, 0.4f,  1.5f}, {1.0f, 0.0f, 0.0f}, {1.0f, 1.0f} },
	{ {0.0f, 0.4f, -1.5f}, {1.0f, 0.0f, 0.0f}, {1.0f, 0.0f} }
};


void configureMockupData()
{
	// PLANE
	glGenVertexArrays(1, &planeVAO);
	glGenBuffers(1, &planeVBO);

	glBindVertexArray(planeVAO);
	glBindBuffer(GL_ARRAY_BUFFER, planeVBO);
	glBufferData(GL_ARRAY_BUFFER, sizeof(planeVertices), planeVertices, GL_STATIC_DRAW);

	// position
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 8 * sizeof(GLfloat),
		(GLvoid*)0);
	glEnableVertexAttribArray(0);

	// normal
	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 8 * sizeof(GLfloat),
		(GLvoid*)(3 * sizeof(GLfloat)));
	glEnableVertexAttribArray(1);

	// texcoord
	glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, 8 * sizeof(GLfloat),
		(GLvoid*)(6 * sizeof(GLfloat)));
	glEnableVertexAttribArray(2);
	glBindVertexArray(0);

	// AREA LIGHT
	glGenVertexArrays(1, &areaLightVAO);
	glBindVertexArray(areaLightVAO);

	glGenBuffers(1, &areaLightVBO);
	glBindBuffer(GL_ARRAY_BUFFER, areaLightVBO);
	glBufferData(GL_ARRAY_BUFFER, sizeof(areaLightVertices), areaLightVertices, GL_STATIC_DRAW);

	// position
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 8 * sizeof(GLfloat),
		(GLvoid*)0);
	glEnableVertexAttribArray(0);

	// normal
	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 8 * sizeof(GLfloat),
		(GLvoid*)(3 * sizeof(GLfloat)));
	glEnableVertexAttribArray(1);

	// texcoord
	glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, 8 * sizeof(GLfloat),
		(GLvoid*)(6 * sizeof(GLfloat)));
	glEnableVertexAttribArray(2);
	glBindVertexArray(0);

}
void renderPlane()
{
	glBindVertexArray(planeVAO);
	glDrawArrays(GL_TRIANGLES, 0, 6);
	glBindVertexArray(0);
}
void renderAreaLight()
{
	glBindVertexArray(areaLightVAO);
	glDrawArrays(GL_TRIANGLES, 0, 6);
	glBindVertexArray(0);
}

void framebuffer_size_callback(GLFWwindow* window, int width, int height);
void key_callback(GLFWwindow* window, int key, int scancode, int action, int mode);
void mouse_callback(GLFWwindow* window, double xpos, double ypos);
void scroll_callback(GLFWwindow* window, double xoffset, double yoffset);
void do_movement(GLfloat deltaTime);
int SaveScreenshot(const char* filename) {

  char* data = (char*)malloc(
      (size_t)(SCR_WIDTH * SCR_HEIGHT * 4));  // 3 components (R, G, B)

  if (!data) return 0;

  glPixelStorei(GL_PACK_ALIGNMENT, 1);
  glReadPixels(0, 0, SCR_WIDTH, SCR_HEIGHT, GL_RGBA, GL_UNSIGNED_BYTE, data);
  stbi_flip_vertically_on_write(1);

  int saved = stbi_write_png(filename, SCR_WIDTH, SCR_HEIGHT, 4, data, 0);
  free(data);

  return saved;
}


struct LTC_matrices {

	GLuint mat1;
	GLuint mat2;
	GLuint mat;
};

GLuint loadLDTTexture() {
	GLuint texture = 0;
	glGenTextures(1, &texture);
	glBindTexture(GL_TEXTURE_2D, texture);

	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, intensityDis[0].size(),
                     intensityDis.size(), 0, GL_RGBA, GL_FLOAT, LDTLUT);

	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

	glBindTexture(GL_TEXTURE_2D, 0);
	return texture;
}

GLuint loadMTexture()
{
	GLuint texture = 0;
	glGenTextures(1, &texture);
	glBindTexture(GL_TEXTURE_2D, texture);

	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, 64, 64,
		0, GL_RGBA, GL_FLOAT, LTC1);

	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

	glBindTexture(GL_TEXTURE_2D, 0);
	return texture;
}
GLuint loadLUTTexture()
{
	GLuint texture = 0;
	glGenTextures(1, &texture);
	glBindTexture(GL_TEXTURE_2D, texture);

	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, 64, 64,
		0, GL_RGBA, GL_FLOAT, LTC2);

	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

	glBindTexture(GL_TEXTURE_2D, 0);
	return texture;
}

void incrementRoughness(float step)
{
	static glm::vec3 color = Color::White;
	roughness += step;
    roughness = glm::clamp(roughness, 0.0f, 1.0f);
    float roughnessChousei = 1 - pow((1 - roughness), 3.3);
	//std::cout << "roughness: " << roughness << '\n';
	ltcShaderPtr->use();
        ltcShaderPtr->setVec4("material.albedoRoughness",
                              glm::vec4(color, roughnessChousei));
	glUseProgram(0);
}
void incrementLightIntensity(float step)
{
	static float intensity = 1.0f;
	intensity += step;
	intensity = glm::clamp(intensity, 0.0f, 10.0f);
	//std::cout << "intensity: " << intensity << '\n';
	ltcShaderPtr->use();
	ltcShaderPtr->setFloat("areaLight.intensity", intensity);
	glUseProgram(0);
}
void switchTwoSided(bool doSwitch)
{
	static bool twoSided = true;
	if (doSwitch) twoSided = !twoSided;
	//std::cout << "twoSided: " << std::boolalpha << twoSided << '\n';
	ltcShaderPtr->use();
	ltcShaderPtr->setFloat("areaLight.twoSided", twoSided);
	glUseProgram(0);
}

void readLDT() {
	for (auto& innerVec : intensityDis)
		innerVec.clear();
	intensityDis.clear();
	std::string err;
	std::string warn;
	if (!tiny_ldt<float>::load_ldt("../../../photometry/" + lightType + ".LDT", err, warn, ldt)) {
		std::cout << "failed" << std::endl;
	}
	if (!err.empty())
		std::cout << err << std::endl;
	if (!warn.empty())
		std::cout << warn << std::endl;

	int cnt = 0;
	std::cout << ldt.dc << std::endl;//15
	std::cout << ldt.dg << std::endl;//5

	for (float i = 0.0; i <= 360.0; i += ldt.dc) {
		// int j = i;
		// // 105/15=7 259/37-1=6
		// if((i/ldt.dc)>ldt.luminous_intensity_distribution.size()/((int)(180.0/ldt.dg)+1)-1 &&
		// (int)((ldt.luminous_intensity_distribution.size()/((int)(180.0/ldt.dg)+1)-1)*ldt.dc))
		// // 105 %= (259/37-1)=6*15
		//   j %= (int)((ldt.luminous_intensity_distribution.size()/((int)(180.0/ldt.dg)+1)-1)*ldt.dc);
		// else if((i/ldt.dc)>ldt.luminous_intensity_distribution.size()/((int)(180.0/ldt.dg)+1)-1)
		//   j = 0;
		// if(i==270 && (int)((ldt.luminous_intensity_distribution.size()/((int)(180.0/ldt.dg)+1)-1)*ldt.dc>=90))
		//   j = 90;
		int sz = (180 / ldt.dg) + 1;
		// cout << "i" << i << " j" << j << endl;
		// cout << (int)(j/ldt.dc)*((int)(180.0/ldt.dg)+1) << endl;
		// cout << (int)(j/ldt.dc)*((int)(180.0/ldt.dg)+1)+(int)(180.0/ldt.dg)+1 << endl;
		int st = (int)(i / ldt.dc);
		if (st * sz >= ldt.luminous_intensity_distribution.size())
			st %= ldt.luminous_intensity_distribution.size() / sz;
		intensityDis.emplace_back(
			std::vector<float>(ldt.luminous_intensity_distribution.begin() + st * (sz)
				, ldt.luminous_intensity_distribution.begin() + st * (sz)+sz)
		);
	}
	for (auto v : intensityDis) {
		for (auto p : v) {
			maxLDTValue = std::max(p, maxLDTValue);
			std::cout << p << " ";
		}
		std::cout << std::endl;
	}
	std::cout << intensityDis.size() << std::endl;
}

int main()
{	
	readLDT();
	// OpenGL与窗口的初始化 -------------------------------------------------------------------------------------------------start
	glfwInit();
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 0);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

	GLFWwindow* window = glfwCreateWindow(SCR_WIDTH, SCR_HEIGHT, "Real-Time Area Light photometric light", NULL, NULL);
	if (window == NULL)
	{
		std::cout << "Failed to create GLFW window" << std::endl;
		glfwTerminate();
		return -1;
	}
	glfwMakeContextCurrent(window);

	const int version = gladLoadGL(glfwGetProcAddress);
	if (version == 0) {
		fprintf(stderr, "Failed to load OpenGL 3.x/4.x libraries!\n");
		return 1;
	}
	printf("Load OpenGL %d.%d\n", GLAD_VERSION_MAJOR(version), GLAD_VERSION_MINOR(version));

	glfwMakeContextCurrent(window);
    //glEnable(GL_FRAMEBUFFER_SRGB);
	glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);
	glfwSetCursorPosCallback(window, mouse_callback);
	glfwSetScrollCallback(window, scroll_callback);
	glfwSetKeyCallback(window, key_callback);

	// 隐藏鼠标
	glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_HIDDEN);

	// imgui库
	IMGUI_CHECKVERSION();
	ImGui::CreateContext(); //创建上下文
	ImGuiIO& io = ImGui::GetIO();
	io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard; // 允许键盘控制
	io.ConfigFlags |= ImGuiConfigFlags_NavEnableGamepad;  // 允许游戏手柄控制

	// 设置渲染器后端
	ImGui_ImplGlfw_InitForOpenGL(window, true);
	ImGui_ImplOpenGL3_Init();
	// OpenGL与窗口的初始化 -------------------------------------------------------------------------------------------------end

	// 预处理 -------------------------------------------------------------------------------------------------start
	glEnable(GL_DEPTH_TEST);
	

	// SHADERS
	Shader shaderPlane("../../../glsl/cubature.vert", "../../../glsl/cubature.frag");
	
	// Shader shaderPlane("../../../glsl/plane.vert", "../../../glsl/plane.frag");
	ltcShaderPtr = &shaderPlane;
	Shader shaderLight("../../../glsl/area_light.vert", "../../../glsl/area_light.frag");

	// SHADER CONFIGURATION
	shaderPlane.use();
	shaderPlane.setVec3("Vertices[0]", areaLightVertices[0].position);
	shaderPlane.setVec3("Vertices[1]", areaLightVertices[1].position);
	shaderPlane.setVec3("Vertices[2]", areaLightVertices[4].position);
	shaderPlane.setVec3("Vertices[3]", areaLightVertices[5].position);
	std::cout << intensityDis.size() * intensityDis[0].size() << std::endl;
	int iter = 0;
	for(int i = 0; i < intensityDis.size(); i ++)
		for (int j = 0; j < intensityDis[0].size(); j ++)
		{
			// std::string name = "intensityDis[" + std::to_string(i * intensityDis[0].size() + j) + "]";
			// shaderPlane.setFloat(name, intensityDis[i][j]);
			LDTLUT[iter++] = intensityDis[i][j] / maxLDTValue;
			LDTLUT[iter++] = intensityDis[i][j] / maxLDTValue;
			LDTLUT[iter++] = intensityDis[i][j] / maxLDTValue;
			LDTLUT[iter++] = intensityDis[i][j] / maxLDTValue;
		}
	
	// LUT textures
	LTC_matrices mLTC;
	
    mLTC.mat = loadLDTTexture();
	mLTC.mat1 = loadMTexture();
    mLTC.mat2 = loadLUTTexture();
	
	shaderPlane.setFloat("ldtdc", ldt.dc);
	shaderPlane.setFloat("ldtdg", ldt.dg);

	shaderPlane.setVec3("PolygonNormal", areaLightVertices[0].normal);
	shaderPlane.setInt("VertexCount", 4);
	shaderPlane.setFloat("PolygonArea", 9.f);
	shaderPlane.setFloat("IntensityMulti", IntensityMulti);
	shaderPlane.setFloat("maxLDTValue", maxLDTValue);
	shaderPlane.setInt("LDTLUT", 0);
	shaderPlane.setFloat("LUT_SIZE_X", intensityDis.size());
	shaderPlane.setFloat("LUT_SIZE_Y", intensityDis[0].size());

	shaderPlane.setVec3("areaLight.points[0]", areaLightVertices[0].position);
	shaderPlane.setVec3("areaLight.points[1]", areaLightVertices[1].position);
	shaderPlane.setVec3("areaLight.points[2]", areaLightVertices[4].position);
	shaderPlane.setVec3("areaLight.points[3]", areaLightVertices[5].position);
	shaderPlane.setVec3("areaLight.color", LIGHT_COLOR);
	shaderPlane.setInt("LTC1", 1);
	shaderPlane.setInt("LTC2", 2);
	// shaderPlane.setInt("material.diffuse", 2);
	
	incrementRoughness(0.0f);
	// incrementLightIntensity(0.0f);
	switchTwoSided(false);

	glUseProgram(0);

	shaderLight.use();
	glm::mat4 model(1.0f);
	shaderLight.setMat4("model", model);
	shaderLight.setVec3("lightColor", LIGHT_COLOR);
	glUseProgram(0);

	// 3D OBJECTS
	configureMockupData();
	areaLightTranslate = glm::vec3(0.0f, 0.0f, 0.0f);
    std::random_device rd;  
    std::mt19937 gen(rd());  
    std::uniform_real_distribution<> dis(0.0, 1.0);
	// 预处理 -------------------------------------------------------------------------------------------------end
	
	while (!glfwWindowShouldClose(window))
	{
		// 处理时间、输入以及清屏 --------------------------------------------------------------------------------------------------start
		float currentFrame = static_cast<float>(glfwGetTime());
		deltaTime = currentFrame - lastFrame;
		lastFrame = currentFrame;

		glfwPollEvents();
		do_movement(deltaTime);

        
		glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		ImGui_ImplOpenGL3_NewFrame();
		ImGui_ImplGlfw_NewFrame();
		ImGui::NewFrame();
		App::IntensityMultiPtr = &IntensityMulti;
		App::CameraPosition = &camera.Position;
		App::CameraYaw = &camera.Yaw;
		App::CameraPitch = &camera.Pitch;
        App::Roughness = &roughness;
		App::RenderUI();
		if (App::lightChanged)
		{
			lightType = App::lightType;
			iter = maxLDTValue = 0;
			readLDT();
			shaderPlane.use();
			for (int i = 0; i < intensityDis.size(); i++)
				for (int j = 0; j < intensityDis[0].size(); j++)
					maxLDTValue = std::max(maxLDTValue, intensityDis[i][j]);
			for (int i = 0; i < intensityDis.size(); i++)
				for (int j = 0; j < intensityDis[0].size(); j++)
				{
					LDTLUT[iter++] = intensityDis[i][j] / maxLDTValue;
					LDTLUT[iter++] = intensityDis[i][j] / maxLDTValue;
					LDTLUT[iter++] = intensityDis[i][j] / maxLDTValue;
					LDTLUT[iter++] = intensityDis[i][j] / maxLDTValue;
				}
			mLTC.mat = loadLDTTexture();
            mLTC.mat1 = loadMTexture();
            mLTC.mat2 = loadLUTTexture();
			shaderPlane.setFloat("ldtdc", ldt.dc);
			shaderPlane.setFloat("ldtdg", ldt.dg);
			shaderPlane.setFloat("maxLDTValue", maxLDTValue);
			shaderPlane.setInt("LDTLUT", 0);
            shaderPlane.setInt("LTC1", 1);
            shaderPlane.setInt("LTC2", 2);
			shaderPlane.setFloat("LUT_SIZE_X", intensityDis.size());
			shaderPlane.setFloat("LUT_SIZE_Y", intensityDis[0].size());
            float randomValue1 = dis(gen);
            float randomValue2 = dis(gen);
            shaderPlane.setFloat("randNum1", randomValue1);
            shaderPlane.setFloat("randNum2", randomValue2);
			App::lightChanged = false; 
		}
			
		// 处理时间、输入以及清屏 --------------------------------------------------------------------------------------------------end

		shaderPlane.use();
		shaderPlane.setFloat("IntensityMulti", IntensityMulti);
		glm::mat4 model(1.0f);
		glm::mat3 normalMatrix = glm::mat3(model);
		shaderPlane.setMat4("model", model);
		shaderPlane.setMat3("normalMatrix", normalMatrix);
		glm::mat4 view = camera.GetViewMatrix();
		shaderPlane.setMat4("view", view);
		glm::mat4 projection = glm::perspective(
		glm::radians(camera.Zoom), (float)SCR_WIDTH / (float)SCR_HEIGHT, 0.1f, 300.0f);
		shaderPlane.setMat4("projection", projection);
		
		shaderPlane.setVec3("viewPosition", camera.Position);
		// ADD
		shaderPlane.setVec3("areaLightTranslate", areaLightTranslate);

		glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D, mLTC.mat);  
		glActiveTexture(GL_TEXTURE1);
		glBindTexture(GL_TEXTURE_2D, mLTC.mat1);
		glActiveTexture(GL_TEXTURE2);
		glBindTexture(GL_TEXTURE_2D, mLTC.mat2);

		renderPlane();
		glUseProgram(0);

		shaderLight.use();
		model = glm::translate(model, areaLightTranslate);
		shaderLight.setMat4("model", model);
		shaderLight.setMat4("view", view);
		shaderLight.setMat4("projection", projection);
		renderAreaLight();
		glUseProgram(0);

		ImGui::Render();
		ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

		if(screenShotFlag){
            screenShotFlag = false;
            SaveScreenshot("output.png");
		}

		glfwSwapBuffers(window);
	}
	glfwTerminate();
	return 0;
}

void key_callback(GLFWwindow* window, int key, int scancode, int action, int mode)
{
	static unsigned short wireframe = 0;

	if (action == GLFW_PRESS)
	{
		switch (key)
		{
		case GLFW_KEY_ESCAPE:
			glfwSetWindowShouldClose(window, GL_TRUE);
			return;
		case GLFW_KEY_B:
			switchTwoSided(true);
			break;
		default:
			keys[key] = true;
			break;
		}
	}

	if (action == GLFW_RELEASE)
	{
		if (key == GLFW_KEY_SPACE) {
			switch (wireframe)
			{
			case 0:
				glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
				wireframe = 1;
				break;
			default:
				glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
				wireframe = 0;
				break;
			}
		}
		else {
			keys[key] = false;
		}
	}
}
void do_movement(GLfloat deltaTime)
{
	float cameraSpeed = deltaTime * 3.0f;

	if (keys[GLFW_KEY_W]) {
		camera.ProcessKeyboard(FORWARD, cameraSpeed);
	}
	else if (keys[GLFW_KEY_S]) {
		camera.ProcessKeyboard(BACKWARD, cameraSpeed);
	}
	if (keys[GLFW_KEY_A]) {
		camera.ProcessKeyboard(LEFT, cameraSpeed);
	}
	else if (keys[GLFW_KEY_D]) {
		camera.ProcessKeyboard(RIGHT, cameraSpeed);
	}

	if (keys[GLFW_KEY_R]) {
		if (keys[GLFW_KEY_LEFT_SHIFT]) incrementRoughness(0.01f);
		else incrementRoughness(-0.01f);
	}

	if (keys[GLFW_KEY_I]) {
		if (keys[GLFW_KEY_LEFT_SHIFT]) incrementLightIntensity(0.025f);
		else incrementLightIntensity(-0.025f);
	}

	if (keys[GLFW_KEY_LEFT]) {
		areaLightTranslate.z += 0.01f;
	}
	if (keys[GLFW_KEY_RIGHT]) {
		areaLightTranslate.z -= 0.01f;
	}
	if (keys[GLFW_KEY_UP]) {
		areaLightTranslate.y += 0.01f;
	}
	if (keys[GLFW_KEY_DOWN]) {
		areaLightTranslate.y -= 0.01f;
	}
}
// glfw: whenever the window size changed (by OS or user resize) this callback function executes
void framebuffer_size_callback(GLFWwindow* window, int width, int height)
{
	// make sure the viewport matches the new window dimensions; note that width and 
	// height will be significantly larger than specified on retina displays.
	glViewport(0, 0, width, height);
}
// glfw: whenever the mouse moves, this callback is called
void mouse_callback(GLFWwindow* window, double xposIn, double yposIn)
{
	float xpos = static_cast<float>(xposIn);
	float ypos = static_cast<float>(yposIn);

	if (firstMouse)
	{
		lastX = xpos;
		lastY = ypos;
		firstMouse = false;
	}

	float xoffset = xpos - lastX;
	float yoffset = lastY - ypos; // reversed since y-coordinates go from bottom to top

	lastX = xpos;
	lastY = ypos;

	camera.ProcessMouseMovement(xoffset, yoffset);
}
// glfw: whenever the mouse scroll wheel scrolls, this callback is called
void scroll_callback(GLFWwindow* window, double xoffset, double yoffset)
{
	camera.ProcessMouseScroll(static_cast<float>(yoffset));
}

