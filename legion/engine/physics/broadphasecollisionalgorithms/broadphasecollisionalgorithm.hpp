#pragma once
#include <core/core.hpp>
#include <physics/components/physics_component.hpp>
#include <physics/data/physics_manifold_precursor.hpp>
#include <physics/data/physics_manifold.hpp>

namespace legion::physics
{
    /**@class BroadPhaseCollisionAlgorithm
     * @brief The base class for a broad phase collision detection
     */
    class BroadPhaseCollisionAlgorithm
    {
    public:
        BroadPhaseCollisionAlgorithm()
        {

        }

        /**@brief Collects collider pairs that have a chance of colliding and should be checked in narrow-phase collision detection
         * @param manifoldPrecursors all the physics components 
         * @param manifoldPrecursorGrouping a list-list of colliders that have a chance of colliding and should be checked
         */
        virtual void collectPairs(std::vector<physics_manifold_precursor>&& manifoldPrecursors,
            std::vector<std::vector<physics_manifold_precursor>>& manifoldPrecursorGrouping) = 0;
    };
}
