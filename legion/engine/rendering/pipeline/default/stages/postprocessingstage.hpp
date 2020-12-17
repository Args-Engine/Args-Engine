#pragma once
#include <rendering/pipeline/base/renderstage.hpp>
#include <rendering/pipeline/base/pipeline.hpp>
#include <rendering/data/vertexarray.hpp>
#include <rendering/data/buffer.hpp>
#include <rendering/data/postprocessingeffect.hpp>

namespace legion::rendering
{
    /**
     * @class PostProcessingStage
     * @brief The stage that runs all of the post processing effects.
     */
    class PostProcessingStage : public RenderStage<PostProcessingStage>
    {
    private:
        /**
         * @brief A multimap with priority as key and postprocessing effect as value.
         */
        static std::multimap<priority_type, std::unique_ptr<PostProcessingEffectBase>,std::greater<>> m_effects;
        vertexarray m_quadVAO;
        buffer m_quadVBO;

        framebuffer m_drawFBO;

        texture_handle m_swapTexture;

        shader_handle m_screenShader;

    public:
        template<typename effect_type, typename ...Args, inherits_from<effect_type, PostProcessingEffect<effect_type>> = 0>
        static void addEffect(priority_type priority = default_priority, Args&&...args)
        {
            m_effects.emplace(priority, std::unique_ptr<PostProcessingEffectBase>(new effect_type(std::forward<Args>(args)...)));
        }

        template<typename effect_type, inherits_from<effect_type, PostProcessingEffect<effect_type>> = 0>
        static void removeEffect()
        {
            for (auto iter = m_effects.begin(); iter != m_effects.end();)
            {
                const auto eraseIter = iter++;
                if (eraseIter->second->getId() == effect_type::m_id)
                {
                    m_effects.erase(eraseIter);
                }
            }
        }

        virtual void setup(app::window& context) override;
        virtual void render(app::window& context, camera& cam, const camera::camera_input& camInput, time::span deltaTime) override;
        virtual priority_type priority() override;
    };

}

