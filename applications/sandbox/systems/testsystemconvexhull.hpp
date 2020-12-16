#pragma once
#include <core/core.hpp>
#include <application/application.hpp>
#include <core/filesystem/filesystem.hpp>
#include <rendering/data/material.hpp>
#include <core/logging/logging.hpp>
#include <rendering/debugrendering.hpp>

#include <rendering/components/renderable.hpp>
#include <physics/components/physics_component.hpp>


#include "pointcloud_particlesystem.hpp"
#include <rendering/components/particle_emitter.hpp>

using namespace legion;

struct convex_hull_step : public app::input_action<convex_hull_step> {};
struct convex_hull_draw : public app::input_action<convex_hull_draw> {};
struct convex_hull_info : public app::input_action<convex_hull_info> {};

class TestSystemConvexHull final : public System<TestSystemConvexHull>
{
public:
    std::shared_ptr<legion::physics::ConvexCollider> collider = nullptr;

    TestSystemConvexHull()
    {
        log::filter(log::severity::debug);
        app::WindowSystem::requestWindow(world_entity_id, math::ivec2(1360, 768), "LEGION Engine", "Legion Icon", nullptr, nullptr, 1);
    }

    core::ecs::entity_handle physicsEnt;
    core::mesh_handle meshH;

    virtual void setup()
    {
        app::InputSystem::createBinding<convex_hull_step>(app::inputmap::method::ENTER);
        app::InputSystem::createBinding<convex_hull_draw>(app::inputmap::method::M);
        bindToEvent<convex_hull_step, &TestSystemConvexHull::convexHullStep>();
        bindToEvent<convex_hull_draw, &TestSystemConvexHull::convexHullDraw>();

        createProcess<&TestSystemConvexHull::update>("Update");

        rendering::model_handle cube;
        rendering::model_handle model;

        rendering::material_handle flatGreen;
        rendering::material_handle vertexColor;
        rendering::material_handle directionalLightMH;
        rendering::material_handle wireFrameH;

        app::window window = m_ecs->world.get_component_handle<app::window>().read();
        {
            async::readwrite_guard guard(*window.lock);
            app::ContextHelper::makeContextCurrent(window);

            cube = rendering::ModelCache::create_model("cube", "assets://models/cube.obj"_view);
            model = rendering::ModelCache::create_model("model", "assets://models/suzanne.obj"_view);
            wireFrameH = rendering::MaterialCache::create_material("wireframe", "assets://shaders/wireframe.shs"_view);
            vertexColor = rendering::MaterialCache::create_material("vertexColor", "assets://shaders/vertexcolor.shs"_view);

            // Create physics entity
            {
                physicsEnt = createEntity();
                physicsEnt.add_components<rendering::renderable>({ model, wireFrameH });
                physicsEnt.add_components<transform>(position(0, 0, 0), rotation(), scale(1));
                meshH = model.get_mesh();

                physicsEnt.add_component<physics::physicsComponent>();
               /* auto rbH = ent.add_component<physics::rigidbody>();
                auto rb = rbH.read();
                rb.setMass(1.0f);
                rbH.write(rb);*/
            }
            // Create physics entity
            /*{
                auto ent = createEntity();
                ent.add_components<rendering::renderable>({ cube, wireFrameH });
                ent.add_components<transform>(position(0, -5, 0), rotation(), scale(2));
                auto pcH = ent.add_component<physics::physicsComponent>();
                auto pc = pcH.read();

                pc.AddBox(physics::cube_collider_params(1.0f, 1.0f, 1.0f));
                pcH.write(pc);
            }*/
            // Create entity for reference
            {
                auto ent = createEntity();
                ent.add_components<rendering::renderable>({ cube, vertexColor });
                ent.add_components<transform>(position(5.0f, 0, 0), rotation(), scale(1));
            }
            // Create entity for reference
            {
                auto ent = createEntity();
                ent.add_components<rendering::renderable>({ cube, wireFrameH });
                ent.add_components<transform>(position(0, 0, -5.0f), rotation(), scale(1));
            }
            {
                math::vec3 p = math::vec3(11, 4, 5);
                std::vector<math::vec3> points{ math::vec3(0,0,0), math::vec3(2, 5, 5), math::vec3(12, 5, 5), math::vec3(10, 0, 0) };
                math::mat4 planeMat = math::planeMatrix(points.at(0), points.at(1), points.at(3), math::vec3(0,0,0));
                std::vector<math::vec3> newPoints{
                    math::vec3(math::inverse(planeMat) * math::vec4(points.at(0),1)),
                    math::vec3(math::inverse(planeMat) * math::vec4(points.at(1),1)),
                    math::vec3(math::inverse(planeMat) * math::vec4(points.at(2),1)),
                    math::vec3(math::inverse(planeMat) * math::vec4(points.at(3),1))
                };
                log::debug("mapped: {} {} {} {}", newPoints.at(0), newPoints.at(1), newPoints.at(2), newPoints.at(3));
            }
        }
    }

    void update(time::span deltaTime)
    {
        //physics::PhysicsSystem::IsPaused = false;

        //drawPhysicsColliders();
    }

    void drawPhysicsColliders()
    {
        static auto physicsQuery = createQuery< physics::physicsComponent>();

        for (auto entity : physicsQuery)
        {
            auto rotationHandle = entity.get_component_handle<rotation>();
            auto positionHandle = entity.get_component_handle<position>();
            auto scaleHandle = entity.get_component_handle<scale>();
            auto physicsComponentHandle = entity.get_component_handle<physics::physicsComponent>();

            bool hasTransform = rotationHandle && positionHandle && scaleHandle;
            bool hasNecessaryComponentsForPhysicsManifold = hasTransform && physicsComponentHandle;

            if (hasNecessaryComponentsForPhysicsManifold)
            {
                auto rbColor = math::color(0.0, 0.5, 0, 1);
                auto statibBlockColor = math::color(0, 1, 0, 1);

                rotation rot = rotationHandle.read();
                position pos = positionHandle.read();
                scale scale = scaleHandle.read();

                auto usedColor = statibBlockColor;
                bool useDepth = false;

                if (entity.get_component_handle<physics::rigidbody>())
                {
                    usedColor = rbColor;
                    useDepth = true;
                }


                //assemble the local transform matrix of the entity
                math::mat4 localTransform;
                math::compose(localTransform, scale, rot, pos);

                auto physicsComponent = physicsComponentHandle.read();

                for (auto physCollider : *physicsComponent.colliders)
                {
                    //--------------------------------- Draw Collider Outlines ---------------------------------------------//

                    for (auto face : physCollider->GetHalfEdgeFaces())
                    {
                        //face->forEachEdge(drawFunc);
                        physics::HalfEdgeEdge* initialEdge = face->startEdge;
                        physics::HalfEdgeEdge* currentEdge = face->startEdge;

                        math::vec3 faceStart = localTransform * math::vec4(face->centroid, 1);
                        math::vec3 faceEnd = faceStart + math::vec3((localTransform * math::vec4(face->normal, 0)));

                        debug::drawLine(faceStart, faceEnd, math::colors::green, 5.0f);

                        if (!currentEdge) { return; }

                        do
                        {
                            physics::HalfEdgeEdge* edgeToExecuteOn = currentEdge;
                            currentEdge = currentEdge->nextEdge;

                            math::vec3 worldStart = localTransform * math::vec4(edgeToExecuteOn->edgePosition, 1);
                            math::vec3 worldEnd = localTransform * math::vec4(edgeToExecuteOn->nextEdge->edgePosition, 1);

                            debug::drawLine(worldStart, worldEnd, usedColor, 2.0f, 0.0f, useDepth);

                        } while (initialEdge != currentEdge && currentEdge != nullptr);
                    }
                }

            }

        }
    }

    void convexHullStep(convex_hull_step* action)
    {
        if (action->value)
        {
            auto pc = physicsEnt.read_component<physics::physicsComponent>();
            collider = pc.ConstructConvexHull(meshH);
            physicsEnt.write_component(pc);
        }
    }

    void convexHullDraw(convex_hull_draw* action)
    {
        if (action->value)
        {
                auto debugDrawEdges = [](legion::physics::HalfEdgeEdge* edge)
                {
                    if (!edge || !edge->nextEdge) return;
                    math::vec3 pos0 = edge->edgePosition;
                    math::vec3 pos1 = edge->nextEdge->edgePosition;
                    debug::drawLine(pos0, pos1, math::colors::red);
                };

                auto faces = collider->GetHalfEdgeFaces();
                for (int i = 0; i < faces.size(); ++i)
                {
                    faces.at(i)->forEachEdge(debugDrawEdges);
                    // Draw normals
                    debug::drawLine(faces.at(i)->centroid, faces.at(i)->centroid + faces.at(i)->normal * 0.3f, math::colors::white);
                }
        }
    }
};
