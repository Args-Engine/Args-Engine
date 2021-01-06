#pragma once
#include <rendering/data/particle_system_base.hpp>
#include <rendering/debugrendering.hpp>
#include<core/core.hpp>
#include <rendering/data/Octree.hpp>
using namespace legion;
/**
 * @struct pointCloudParameters
 * @brief A struct that simplifies the parameter input of the particle system constructor.
 */
struct pointCloudParameters
{
    math::vec3 startingSize;
    rendering::material_handle particleMaterial;
    rendering::model_handle particleModel;

    bool looping = false;
    float maxLifeTime = 0.0f;
    float startingLifeTime = 0.0f;
    uint spawnRate = 0;
    uint particleCount = 0;
    uint maxParticles = 0;
    math::vec3 startingVel = math::vec3(0);
    float sizeOverLifeTime = 0.0f;
};

class PointCloudParticleSystem : public rendering::ParticleSystemBase
{
public:
    /**
     * @brief Constructor of the point cloud particle system.
     * @param params A struct with a bunch of default parameters and some parameters needed to be set.
     * @param positions A list of positions that the particle system uses to create its particles at.
     */
    PointCloudParticleSystem(pointCloudParameters params, const std::vector<math::vec3>& positions)
    {
        m_looping = params.looping;
        m_maxLifeTime = params.maxLifeTime;
        m_startingLifeTime = params.startingLifeTime;
        m_spawnRate = params.spawnRate;
        m_particleCount = params.particleCount;
        m_maxParticles = params.maxParticles;
        m_startingVelocity = params.startingVel;
        m_startingSize = params.startingSize;
        m_sizeOverLifetime = params.sizeOverLifeTime;
        m_particleMaterial = params.particleMaterial;
        m_particleModel = params.particleModel;
        m_positions = positions;
    }
    /**
     * @brief Setup function that will be called to populate the emitter with the required particles.
     * @param emitter_handle The emitter that you are populating.
     */
    void setup(ecs::component_handle<rendering::particle_emitter> emitter_handle) const override
    {
        //Reads emitter.
        rendering::particle_emitter emitter = emitter_handle.read();

        float minX = std::numeric_limits<float>().max();
        float maxX = std::numeric_limits<float>().min();

        float minY = std::numeric_limits<float>().max();
        float maxY = std::numeric_limits<float>().min();

        float minZ = std::numeric_limits<float>().max();
        float maxZ = std::numeric_limits<float>().min();

        for (auto position : m_positions)
        {
            //Checks the emitter if it has a recycled particle to use, if not it creates a new one.
            ecs::component_handle<rendering::particle> particleComponent = checkToRecycle(emitter_handle);
            auto ent = particleComponent.entity;
            //Checks if the entity has a transform, if not it adds one.
            if (!ent.has_components<transform>())
                ent.add_components<transform>();

            //Gets position, rotation and scale of entity.
            auto trans = ent.get_component_handles<transform>();
            auto& [pos, _, scale] = trans;
            //auto currentPos = pos.read();
            if (position.x < minX) minX = position.x;
            if (position.x > maxX) maxX = position.x;
            if (position.y < minY) minY = position.y;
            if (position.y > maxY) maxY = position.y;
            if (position.z < minZ) minZ = position.z;
            if (position.z > maxZ) maxZ = position.z;

            //Sets the particle scale to the right scale.
            pos.write(position);
            scale.write(math::vec3(m_startingSize));

            //Populates the particle with the appropriate stuffs.
            createParticle(particleComponent, trans);
        }
        math::vec3 A = math::vec3(minX, minY, minZ);
        math::vec3 B = math::vec3(minX, minY, maxZ);
        //   debug::drawCube(math::vec3(minX, minY, minZ), math::vec3(maxX, maxY, maxZ), math::colors::red, 1.0f);
        auto tree = rendering::Octree<uint8>(8, math::vec3(0, 0, 0), math::vec3(3, 3, 3));

        tree.insertNode(0, math::vec3(0.0f, 2.5f, 0.5f));
        tree.insertNode(0, math::vec3(0.0f, 2.5f, 2.5f));
        tree.insertNode(0, math::vec3(0.0f, 2.5f, 0.5f));
        tree.insertNode(0, math::vec3(0.0f, 2.5f, 0.5f));
        tree.insertNode(0, math::vec3(0.0f, 2.5f, 0.5f));
        tree.insertNode(0, math::vec3(0.0f, 2.5f, 0.5f));
        tree.insertNode(0, math::vec3(0.0f, 2.5f, 0.5f));
        tree.insertNode(0, math::vec3(0.0f, 2.5f, 0.5f));
        tree.insertNode(0, math::vec3(0.0f, 2.5f, 0.5f));
        tree.insertNode(0, math::vec3(0.0f, 2.5f, 0.5f));
        tree.insertNode(0, math::vec3(0.0f, 2.5f, 0.5f));


        tree.DrawTree();
    }

    void update(std::vector<ecs::entity_handle>, ecs::component_handle<rendering::particle_emitter>, time::span) const override
    {

    }

private:
    std::vector<math::vec3> m_positions;
};
