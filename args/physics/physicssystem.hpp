#pragma once

#ifndef ARGS_IMPORT
#define ARGS_IMPORT
#include <core/core.hpp>
#include <core/platform/args_library.hpp>
#else
#include <core/core.hpp>
#endif // !ARGS_IMPORT



namespace args::physics
{
    class PhysicsSystem final : public System<PhysicsSystem>
    {
    public:

        virtual void setup()
        {
            createProcess<&PhysicsSystem::fixedUpdate>("Physics", m_timeStep);
        }

        void fixedUpdate(time::time_span<fast_time> deltaTime)
        {
            runPhysicsPipeline();

            integrateRigidbodies(deltaTime);

        }

    private:

        const float m_timeStep = 0.02f;

        /** @brief Performs the entire physics pipeline (
         * Broadphase Collision Detection, Narrowphase Collision Detection, and the Collision Resolution)
        */
        void runPhysicsPipeline()
        {
            //Broadphase Optimization

            //Narrophase 

            //Collision Resolution

        }

        /** @brief gets all the entities with a rigidbody component and calls the integrate function on them
        */
        void integrateRigidbodies(float deltaTime)
        {

        }
    };
}
