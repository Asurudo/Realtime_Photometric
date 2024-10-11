#include <iostream>
#include <string>
#include <glm/glm.hpp>
#include "imgui/imgui.h"

namespace App
{
    //选中结果
    bool lightChanged = false;
    //字符串结果
    std::string lightType = "LINETIK-S_42184482";
    //拖拽值
    float *IntensityMultiPtr;
    glm::vec3 *CameraPosition;
    float *CameraYaw, *CameraPitch;
    float *Roughness;

    void RenderUI()
    {
        //显示一些文本（也可以使用字符串格式）
        ImGui::Text("Avg fps: %.3f", ImGui::GetIO().Framerate);
        
        //创建一个设置窗口
        // ImGui::Begin("slider button");
        //按钮在单击时返回true（大多数小部件在编辑/激活时返回true）
        if (ImGui::Button("LINETIK-S"))
            lightChanged = true, lightType = "LINETIK-S_42184482";
         if (ImGui::Button("MIREL"))
             lightChanged = true, lightType = "MIREL_42925637";
         if (ImGui::Button("PERLUCE"))
             lightChanged = true, lightType = "PERLUCE_42182932";
         if (ImGui::Button("SLOTLIGHT"))
             lightChanged = true, lightType = "SLOTLIGHT_42184612";
        
        // 缓冲区用于存储文本输入值
        // char buffer[256] = "";
        // ImGui::InputText("input", buffer, sizeof(buffer));

        // ImGui::Checkbox("drag", &isShowDrag);
        // if (isShowDrag)
        // {
        //     float value = 10.0f;
        //     ImGui::DragFloat("val", &value);
        // }
        //使用从0.0f到1.0f的滑块编辑1个浮动
        ImGui::SliderFloat("Exposure", IntensityMultiPtr, 0.0f, 2.0f);

        ImGui::Text("Camera Poisiton: (%.1f, %.1f, %.1f)", CameraPosition->x, CameraPosition->y, CameraPosition->z);
        ImGui::Text("Camera Yaw: %.3f", *CameraYaw);
        ImGui::Text("Camera Pitch: %.3f", *CameraPitch);
        ImGui::Text("Roughness: %.3f", *Roughness);
        // ImGui::SameLine();
        // ImGui::Text("Value %f", *IntensityMultiPtr);
    }
}