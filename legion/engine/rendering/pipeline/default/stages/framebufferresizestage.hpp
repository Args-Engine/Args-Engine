#pragma once
#include <rendering/pipeline/base/renderstage.hpp>
#include <rendering/pipeline/base/pipeline.hpp>

namespace legion::rendering
{
    class FramebufferResizeStage : public RenderStage<FramebufferResizeStage>
    {
        static std::atomic<float> m_renderScale;

        math::ivec2 m_framebufferSize;
        texture_handle m_colorTexture;
        texture_handle m_depthTexture;
        renderbuffer m_stencilbuffer;

    public:
        static void setRenderScale(float renderScale);
        static float getRenderScale();

        virtual void setup(app::window& context) override;
        virtual void render(app::window& context, camera& cam, const camera::camera_input& camInput, time::span deltaTime) override;
        virtual priority_type priority() override;

    };
}
