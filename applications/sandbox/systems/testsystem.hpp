#pragma once
#include <core/core.hpp>
#include <physics/halfedgeedge.hpp>
#include <application/application.hpp>
#include <core/math/math.hpp>

#include <core/logging/logging.hpp>
#include <physics/components/physics_component.hpp>
#include <physics/components/rigidbody.hpp>
#include <physics/cube_collider_params.hpp>
#include <physics/data/physics_manifold_precursor.h>
#include <physics/systems/physicssystem.hpp>
#include <physics/halfedgeface.hpp>
#include <physics/data/penetrationquery.h>

#include <core/compute/context.hpp>
#include <core/compute/kernel.hpp>
#include <core/compute/high_level/function.hpp>
#include <rendering/debugrendering.hpp>

#include <physics/physics_statics.hpp>
#include <physics/data/identifier.hpp>
#include <audio/audio.hpp>
#include <Voro++/voro++.hh>
#include <Voro++/common.hh>

using namespace legion;



struct sah
{
    int value;

    sah operator+(const sah& other)
    {
        return { value + other.value };
    }

    sah operator*(const sah& other)
    {
        return { value * other.value };
    }
};

enum class PhysicsUnitTestMode
{
    CollisionDetectionMode,
    CollisionResolution
};

// Move and strive for the wireframe sphere (which holds the only audio source)
// For testing the movement of audio sources (Spatial audio/doppler)
struct audio_move : public app::input_axis<audio_move> {};
struct audio_strive : public app::input_axis<audio_strive> {};
// Gain and Pitch knob
struct gain_change : public app::input_axis<gain_change> {};
struct pitch_change : public app::input_axis<pitch_change> {};
struct play_audio_source : public app::input_action<play_audio_source> {};
struct pause_audio_source : public app::input_action<pause_audio_source> {};
struct stop_audio_source : public app::input_action<stop_audio_source> {};
struct rewind_audio_source : public app::input_action<rewind_audio_source> {};
struct audio_test_input : public app::input_action<audio_test_input> {};

//some test stuff just so i can change things

struct physics_test_move : public app::input_axis<physics_test_move> {};

struct light_switch : public app::input_action<light_switch> {};
struct tonemap_switch : public app::input_action<tonemap_switch> {};

struct activate_CRtest2 : public app::input_action<activate_CRtest2> {};
struct activate_CRtest3 : public app::input_action<activate_CRtest3> {};

struct activateFrictionTest : public app::input_action<activateFrictionTest> {};

struct extendedPhysicsContinue : public app::input_action<extendedPhysicsContinue> {};
struct nextPhysicsTimeStepContinue : public app::input_action<nextPhysicsTimeStepContinue> {};


//scene loading binds
struct savescene1 : public app::input_action<savescene1> {};
struct savescene2 : public app::input_action<savescene2> {};

struct loadscene1 : public app::input_action<loadscene1> {};
struct loadscene2 : public app::input_action<loadscene2> {};

class TestSystem final : public System<TestSystem>
{
public:
    TestSystem()
    {
        log::filter(log::severity::debug);
        app::WindowSystem::requestWindow(world_entity_id, math::ivec2(1920, 1080), "LEGION Engine", "Legion Icon", nullptr, nullptr, 1); // Create the request for the main window.
    }

    ecs::entity_handle audioSphereLeft;
    ecs::entity_handle audioSphereRight;

    std::vector< ecs::entity_handle > physicsUnitTestCD;
    std::vector< ecs::entity_handle > physicsUnitTestCR;

    //Collision Resolution Tests
    ecs::entity_handle staticToAABBEntLinear;
    ecs::entity_handle staticToAABBEntRotation;
    ecs::entity_handle staticToOBBEnt;
    ecs::entity_handle staticToEdgeEnt;
    ecs::entity_handle staticTo2StackEnt;

    ecs::entity_handle sun;
    rendering::material_handle pbrH;
    rendering::material_handle copperH;
    rendering::material_handle aluminumH;
    rendering::material_handle ironH;
    rendering::material_handle slateH;
    rendering::material_handle rockH;
    rendering::material_handle rock2H;
    rendering::material_handle fabricH;
    rendering::material_handle bogH;
    rendering::material_handle paintH;
    rendering::material_handle skyboxH;
    rendering::material_handle gnomeMH;

    //Friction Test
    std::vector<ecs::entity_handle> physicsFrictionTestRotators;
    bool rotate = false;

    ecs::entity_handle NoFrictionBody;
    ecs::entity_handle Point3FrictionBody;
    ecs::entity_handle Point6FrictionBody;
    ecs::entity_handle FullFrictionBody;



    virtual void setup()
    {



#pragma region Input binding
        app::InputSystem::createBinding<physics_test_move>(app::inputmap::method::LEFT, -1.f);
        app::InputSystem::createBinding<physics_test_move>(app::inputmap::method::RIGHT, 1.f);

        app::InputSystem::createBinding<audio_move>(app::inputmap::method::I, 1.f);
        app::InputSystem::createBinding<audio_move>(app::inputmap::method::K, -1.f);
        app::InputSystem::createBinding<audio_strive>(app::inputmap::method::J, 1.f);
        app::InputSystem::createBinding<audio_strive>(app::inputmap::method::U, 10.f); // Power mode of J
        app::InputSystem::createBinding<audio_strive>(app::inputmap::method::L, -1.f);
        app::InputSystem::createBinding<audio_strive>(app::inputmap::method::O, -10.f); // Power mode of O

        app::InputSystem::createBinding<gain_change>(app::inputmap::method::F10, 0.05f);
        app::InputSystem::createBinding<gain_change>(app::inputmap::method::F9, -0.05f);
        app::InputSystem::createBinding<pitch_change>(app::inputmap::method::APOSTROPHE, 1.0f);
        app::InputSystem::createBinding<pitch_change>(app::inputmap::method::SEMICOLON, -1.0f);

        app::InputSystem::createBinding<play_audio_source>(app::inputmap::method::P);
        app::InputSystem::createBinding<pause_audio_source>(app::inputmap::method::LEFT_BRACKET);
        app::InputSystem::createBinding<stop_audio_source>(app::inputmap::method::RIGHT_BRACKET);
        app::InputSystem::createBinding<rewind_audio_source>(app::inputmap::method::BACKSPACE);
        app::InputSystem::createBinding<audio_test_input>(app::inputmap::method::SLASH);

        app::InputSystem::createBinding< activate_CRtest2>(app::inputmap::method::KP_2);
        app::InputSystem::createBinding< activate_CRtest3>(app::inputmap::method::KP_3);

        app::InputSystem::createBinding< activateFrictionTest >(app::inputmap::method::KP_4);

        app::InputSystem::createBinding<light_switch>(app::inputmap::method::F);
        app::InputSystem::createBinding<tonemap_switch>(app::inputmap::method::G);

        app::InputSystem::createBinding< extendedPhysicsContinue>(app::inputmap::method::M);
        app::InputSystem::createBinding<nextPhysicsTimeStepContinue>(app::inputmap::method::N);

        app::InputSystem::createBinding<savescene1>(app::inputmap::method::F1);
        app::InputSystem::createBinding<savescene2>(app::inputmap::method::F2);
        app::InputSystem::createBinding<loadscene1>(app::inputmap::method::F3);
        app::InputSystem::createBinding<loadscene2>(app::inputmap::method::F4);


        bindToEvent<physics_test_move, &TestSystem::onUnitPhysicsUnitTestMove>();

        bindToEvent<light_switch, &TestSystem::onLightSwitch>();
        bindToEvent<tonemap_switch, &TestSystem::onTonemapSwitch>();

        bindToEvent<audio_move, &TestSystem::onSphereAAMove>();
        bindToEvent<audio_strive, &TestSystem::onSphereAAStrive>();
        bindToEvent<gain_change, &TestSystem::onGainChange>();
        bindToEvent<pitch_change, &TestSystem::onPitchChange>();

        bindToEvent<play_audio_source, &TestSystem::playAudioSource>();
        bindToEvent<pause_audio_source, &TestSystem::pauseAudioSource>();
        bindToEvent<stop_audio_source, &TestSystem::stopAudioSource>();
        bindToEvent<rewind_audio_source, &TestSystem::rewindAudioSource>();

        //bindToEvent<physics::trigger_event, &TestSystem::testPhysicsEvent>();

        bindToEvent<audio_test_input, &TestSystem::audioTestInput>();

        //collision resolution test
        bindToEvent< activate_CRtest2, &TestSystem::onActivateUnitTest2>();
        bindToEvent< activate_CRtest3, &TestSystem::onActivateUnitTest3>();

        //friction test
        bindToEvent< activateFrictionTest, &TestSystem::FrictionTestActivate>();

        bindToEvent< extendedPhysicsContinue, &TestSystem::onExtendedPhysicsContinueRequest>();
        bindToEvent<nextPhysicsTimeStepContinue, &TestSystem::onNextPhysicsTimeStepRequest>();

        //load scene

        bindToEvent<savescene1, &TestSystem::saveScene1>();
        bindToEvent<savescene2, &TestSystem::saveScene2>();
        bindToEvent<loadscene1, &TestSystem::loadScene1>();
        bindToEvent<loadscene2, &TestSystem::loadScene2>();


#pragma endregion              

#pragma region Model and material loading
        rendering::model_handle directionalLightH;
        rendering::model_handle spotLightH;
        rendering::model_handle pointLightH;
        rendering::model_handle audioSourceH;
        rendering::model_handle cubeH;
        rendering::model_handle sphereH;
        rendering::model_handle suzanneH;
        rendering::model_handle gnomeH;
        rendering::model_handle uvsphereH;
        rendering::model_handle axesH;
        rendering::model_handle submeshtestH;
        rendering::model_handle planeH;
        //rendering::model_handle cylinderH;

        rendering::material_handle wireframeH;
        rendering::material_handle vertexColorH;

        rendering::material_handle uvH;
        rendering::material_handle textureH;
        rendering::material_handle directionalLightMH;
        rendering::material_handle spotLightMH;
        rendering::material_handle pointLightMH;
        rendering::material_handle gizmoMH;
        rendering::material_handle normalH;

        app::window window = m_ecs->world.get_component_handle<app::window>().read();

        {
            async::readwrite_guard guard(*window.lock);
            app::ContextHelper::makeContextCurrent(window);

            directionalLightH = rendering::ModelCache::create_model("directional light", "assets://models/directional-light.obj"_view);
            spotLightH = rendering::ModelCache::create_model("spot light", "assets://models/spot-light.obj"_view);
            pointLightH = rendering::ModelCache::create_model("point light", "assets://models/point-light.obj"_view);
            audioSourceH = rendering::ModelCache::create_model("audio source", "assets://models/audio-source.obj"_view);
            cubeH = rendering::ModelCache::create_model("cube", "assets://models/cube.obj"_view);
            sphereH = rendering::ModelCache::create_model("sphere", "assets://models/sphere.obj"_view);
            suzanneH = rendering::ModelCache::create_model("suzanne", "assets://models/suzanne.obj"_view);
            gnomeH = rendering::ModelCache::create_model("gnome", "assets://models/wizardgnomeretop.obj"_view);
            uvsphereH = rendering::ModelCache::create_model("uvsphere", "assets://models/uvsphere.obj"_view);
            axesH = rendering::ModelCache::create_model("axes", "assets://models/xyz.obj"_view);
            submeshtestH = rendering::ModelCache::create_model("submeshtest", "assets://models/submeshtest.obj"_view);
            planeH = rendering::ModelCache::create_model("plane", "assets://models/plane.obj"_view);
            //cylinderH = rendering::ModelCache::create_model("cylinder","assets://models/cylinder.obj"_view);

            wireframeH = rendering::MaterialCache::create_material("wireframe", "assets://shaders/wireframe.shs"_view);
            vertexColorH = rendering::MaterialCache::create_material("vertex color", "assets://shaders/vertexcolor.shs"_view);
            uvH = rendering::MaterialCache::create_material("uv", "assets://shaders/uv.shs"_view);

            auto colorshader = rendering::ShaderCache::create_shader("color", "assets://shaders/color.shs"_view);
            directionalLightMH = rendering::MaterialCache::create_material("directional light", colorshader);
            directionalLightMH.set_param("color", math::color(1, 1, 0.8f));

            spotLightMH = rendering::MaterialCache::create_material("spot light", colorshader);
            spotLightMH.set_param("color", math::colors::green);

            pointLightMH = rendering::MaterialCache::create_material("point light", colorshader);
            pointLightMH.set_param("color", math::colors::red);

            gizmoMH = rendering::MaterialCache::create_material("gizmo", colorshader);
            gizmoMH.set_param("color", math::colors::lightgrey);

            textureH = rendering::MaterialCache::create_material("texture", "assets://shaders/texture.shs"_view);
            textureH.set_param("_texture", rendering::TextureCache::create_texture("engine://resources/default/albedo"_view));

            auto pbrShader = rendering::ShaderCache::create_shader("pbr", "assets://shaders/pbr.shs"_view);
            pbrH = rendering::MaterialCache::create_material("pbr", pbrShader);
            pbrH.set_param(SV_ALBEDO, rendering::TextureCache::create_texture("engine://resources/default/albedo"_view));
            pbrH.set_param(SV_NORMALHEIGHT, rendering::TextureCache::create_texture("engine://resources/default/normalHeight"_view));
            pbrH.set_param(SV_MRDAO, rendering::TextureCache::create_texture("engine://resources/default/MRDAo"_view));
            pbrH.set_param(SV_EMISSIVE, rendering::TextureCache::create_texture("engine://resources/default/emissive"_view));
            pbrH.set_param(SV_HEIGHTSCALE, 1.f);
            pbrH.set_param("discardExcess", false);
            pbrH.set_param("skycolor", math::color(0.2f, 0.4f, 1.0f));
            pbrH.set_param("tonemap", false);

            copperH = rendering::MaterialCache::create_material("copper", pbrShader);
            copperH.set_param("material_input.albedo", rendering::TextureCache::create_texture("assets://textures/copper/copper-albedo.png"_view));
            copperH.set_param("material_input.normalHeight", rendering::TextureCache::create_texture("assets://textures/copper/copper-normalHeight.png"_view));
            copperH.set_param("material_input.MRDAo", rendering::TextureCache::create_texture("assets://textures/copper/copper-MRDAo.png"_view));
            copperH.set_param("material_input.emissive", rendering::TextureCache::create_texture("assets://textures/copper/copper-emissive.png"_view));
            copperH.set_param("material_input.heightScale", 0.1f);
            copperH.set_param("discardExcess", false);
            copperH.set_param("skycolor", math::color(0.2f, 0.4f, 1.0f));
            copperH.set_param("tonemap", false);

            aluminumH = rendering::MaterialCache::create_material("aluminum", pbrShader);
            aluminumH.set_param("material_input.albedo", rendering::TextureCache::create_texture("assets://textures/aluminum/aluminum-albedo.png"_view));
            aluminumH.set_param("material_input.normalHeight", rendering::TextureCache::create_texture("assets://textures/aluminum/aluminum-normalHeight.png"_view));
            aluminumH.set_param("material_input.MRDAo", rendering::TextureCache::create_texture("assets://textures/aluminum/aluminum-MRDAo.png"_view));
            aluminumH.set_param("material_input.emissive", rendering::TextureCache::create_texture("assets://textures/aluminum/aluminum-emissive.png"_view));
            aluminumH.set_param("material_input.heightScale", 0.f);
            aluminumH.set_param("discardExcess", false);
            aluminumH.set_param("skycolor", math::color(0.2f, 0.4f, 1.0f));
            aluminumH.set_param("tonemap", false);

            ironH = rendering::MaterialCache::create_material("iron", pbrShader);
            ironH.set_param("material_input.albedo", rendering::TextureCache::create_texture("assets://textures/iron/rustediron-albedo.png"_view));
            ironH.set_param("material_input.normalHeight", rendering::TextureCache::create_texture("assets://textures/iron/rustediron-normalHeight.png"_view));
            ironH.set_param("material_input.MRDAo", rendering::TextureCache::create_texture("assets://textures/iron/rustediron-MRDAo.png"_view));
            ironH.set_param("material_input.emissive", rendering::TextureCache::create_texture("assets://textures/iron/rustediron-emissive.png"_view));
            ironH.set_param("material_input.heightScale", 0.1f);
            ironH.set_param("discardExcess", false);
            ironH.set_param("skycolor", math::color(0.2f, 0.4f, 1.0f));
            ironH.set_param("tonemap", false);

            slateH = rendering::MaterialCache::create_material("slate", pbrShader);
            slateH.set_param("material_input.albedo", rendering::TextureCache::create_texture("assets://textures/slate/slate-albedo.png"_view));
            slateH.set_param("material_input.normalHeight", rendering::TextureCache::create_texture("assets://textures/slate/slate-normalHeight.png"_view));
            slateH.set_param("material_input.MRDAo", rendering::TextureCache::create_texture("assets://textures/slate/slate-MRDAo.png"_view));
            slateH.set_param("material_input.emissive", rendering::TextureCache::create_texture("assets://textures/slate/slate-emissive.png"_view));
            slateH.set_param("material_input.heightScale", 1.f);
            slateH.set_param("discardExcess", true);
            slateH.set_param("skycolor", math::color(0.2f, 0.4f, 1.0f));
            slateH.set_param("tonemap", false);

            rockH = rendering::MaterialCache::create_material("rock", pbrShader);
            rockH.set_param("material_input.albedo", rendering::TextureCache::create_texture("assets://textures/rock/rock-albedo.png"_view));
            rockH.set_param("material_input.normalHeight", rendering::TextureCache::create_texture("assets://textures/rock/rock-normalHeight.png"_view));
            rockH.set_param("material_input.MRDAo", rendering::TextureCache::create_texture("assets://textures/rock/rock-MRDAo.png"_view));
            rockH.set_param("material_input.emissive", rendering::TextureCache::create_texture("assets://textures/rock/rock-emissive.png"_view));
            rockH.set_param("material_input.heightScale", 1.f);
            rockH.set_param("discardExcess", true);
            rockH.set_param("skycolor", math::color(0.2f, 0.4f, 1.0f));
            rockH.set_param("tonemap", false);

            rock2H = rendering::MaterialCache::create_material("rock 2", pbrShader);
            rock2H.set_param("material_input.albedo", rendering::TextureCache::get_handle("rock-albedo.png"));
            rock2H.set_param("material_input.normalHeight", rendering::TextureCache::get_handle("rock-normalHeight.png"));
            rock2H.set_param("material_input.MRDAo", rendering::TextureCache::get_handle("rock-MRDAo.png"));
            rock2H.set_param("material_input.emissive", rendering::TextureCache::get_handle("rock-emissive.png"));
            rock2H.set_param("material_input.heightScale", 0.5f);
            rock2H.set_param("discardExcess", false);
            rock2H.set_param("skycolor", math::color(0.2f, 0.4f, 1.0f));
            rock2H.set_param("tonemap", false);

            fabricH = rendering::MaterialCache::create_material("fabric", pbrShader);
            fabricH.set_param("material_input.albedo", rendering::TextureCache::create_texture("assets://textures/fabric/fabric-lowres-albedo.png"_view));
            fabricH.set_param("material_input.normalHeight", rendering::TextureCache::create_texture("assets://textures/fabric/fabric-lowres-normalHeight.png"_view));
            fabricH.set_param("material_input.MRDAo", rendering::TextureCache::create_texture("assets://textures/fabric/fabric-lowres-MRDAo.png"_view));
            fabricH.set_param("material_input.emissive", rendering::TextureCache::create_texture("assets://textures/fabric/fabric-lowres-emissive.png"_view));
            fabricH.set_param("material_input.heightScale", 0.1f);
            fabricH.set_param("discardExcess", false);
            fabricH.set_param("skycolor", math::color(0.2f, 0.4f, 1.0f));
            fabricH.set_param("tonemap", false);

            bogH = rendering::MaterialCache::create_material("bog", pbrShader);
            bogH.set_param("material_input.albedo", rendering::TextureCache::create_texture("assets://textures/bog/bog-albedo.png"_view));
            bogH.set_param("material_input.normalHeight", rendering::TextureCache::create_texture("assets://textures/bog/bog-normalHeight.png"_view));
            bogH.set_param("material_input.MRDAo", rendering::TextureCache::create_texture("assets://textures/bog/bog-MRDAo.png"_view));
            bogH.set_param("material_input.emissive", rendering::TextureCache::create_texture("assets://textures/bog/bog-emissive.png"_view));
            bogH.set_param("material_input.heightScale", 0.5f);
            bogH.set_param("discardExcess", true);
            bogH.set_param("skycolor", math::color(0.2f, 0.4f, 1.0f));
            bogH.set_param("tonemap", false);

            paintH = rendering::MaterialCache::create_material("paint", pbrShader);
            paintH.set_param("material_input.albedo", rendering::TextureCache::create_texture("assets://textures/paint/paint-peeling-albedo.png"_view));
            paintH.set_param("material_input.normalHeight", rendering::TextureCache::create_texture("assets://textures/paint/paint-peeling-normalHeight.png"_view));
            paintH.set_param("material_input.MRDAo", rendering::TextureCache::create_texture("assets://textures/paint/paint-peeling-MRDAo.png"_view));
            paintH.set_param("material_input.emissive", rendering::TextureCache::create_texture("assets://textures/paint/paint-peeling-emissive.png"_view));
            paintH.set_param("material_input.heightScale", 0.1f);
            paintH.set_param("discardExcess", false);
            paintH.set_param("skycolor", math::color(0.2f, 0.4f, 1.0f));
            paintH.set_param("tonemap", false);

            gnomeMH = rendering::MaterialCache::create_material("gnome", pbrShader);
            gnomeMH.set_param("material_input.albedo", rendering::TextureCache::create_texture("assets://textures/warlock/warlock-albedo.png"_view));
            gnomeMH.set_param("material_input.normalHeight", rendering::TextureCache::create_texture("assets://textures/warlock/warlock-normalHeight.png"_view));
            gnomeMH.set_param("material_input.MRDAo", rendering::TextureCache::create_texture("assets://textures/warlock/warlock-MRDAo.png"_view));
            gnomeMH.set_param("material_input.emissive", rendering::TextureCache::create_texture("assets://textures/warlock/warlock-emissive.png"_view));
            gnomeMH.set_param("material_input.heightScale", 0.f);
            gnomeMH.set_param("discardExcess", false);
            gnomeMH.set_param("skycolor", math::color(0.2f, 0.4f, 1.0f));
            gnomeMH.set_param("tonemap", false);

            normalH = rendering::MaterialCache::create_material("normal", "assets://shaders/normal.shs"_view);
            normalH.set_param("material_input.normalHeight", rendering::TextureCache::create_texture("engine://resources/default/normalHeight"_view));

            skyboxH = rendering::MaterialCache::create_material("skybox", "assets://shaders/skybox.shs"_view);
            skyboxH.set_param("skycolor", math::color(0.2f, 0.4f, 1.0f));
            skyboxH.set_param("tonemap", false);

            app::ContextHelper::makeContextCurrent(nullptr);
        }
#pragma endregion

#pragma region Entities
        {
            auto ent = createEntity();
            ent.add_components<rendering::renderable>(planeH.get_mesh(), rendering::mesh_renderer(slateH));
            ent.add_components<transform>(position(0, 0.01f, 0), rotation(), scale(10));
                        
        }

        {
            auto ent = createEntity();
            ent.add_components<rendering::renderable>(planeH.get_mesh(), rendering::mesh_renderer(copperH));
            ent.add_components<transform>(position(10, 0.01f, 0), rotation(), scale(10));
        }

        {
            auto ent = createEntity();
            ent.add_components<rendering::renderable>(planeH.get_mesh(), rendering::mesh_renderer(aluminumH));
            ent.add_components<transform>(position(10, 0.01f, 10), rotation(), scale(10));
        }

        {
            auto ent = createEntity();
            ent.add_components<rendering::renderable>(planeH.get_mesh(), rendering::mesh_renderer(ironH));
            ent.add_components<transform>(position(10, 0.01f, -10), rotation(), scale(10));
        }

        {
            auto ent = createEntity();
            ent.add_components<rendering::renderable>(planeH.get_mesh(), rendering::mesh_renderer(fabricH));
            ent.add_components<transform>(position(-10, 0.01f, 0), rotation(), scale(10));
        }

        {
            auto ent = createEntity();
            ent.add_components<rendering::renderable>(planeH.get_mesh(), rendering::mesh_renderer(bogH));
            ent.add_components<transform>(position(-10, 0.01f, 10), rotation(), scale(10));
        }

        {
            auto ent = createEntity();
            ent.add_components<rendering::renderable>(planeH.get_mesh(), rendering::mesh_renderer(pbrH));
            ent.add_components<transform>(position(0, 0.01f, 10), rotation(), scale(10));
        }

        {
            auto ent = createEntity();
            ent.add_components<rendering::renderable>(planeH.get_mesh(), rendering::mesh_renderer(rockH));
            ent.add_components<transform>(position(0, 0.01f, -10), rotation(), scale(10));
        }

        {
            auto ent = createEntity();
            ent.add_components<rendering::renderable>(planeH.get_mesh(), rendering::mesh_renderer(paintH));
            ent.add_components<transform>(position(-10, 0.01f, -10), rotation(), scale(10));
        }

        {
            auto ent = createEntity();
            ent.add_components<rendering::renderable>(directionalLightH.get_mesh(), rendering::mesh_renderer(directionalLightMH));
            ent.add_component<rendering::light>(rendering::light::directional(math::color(1, 1, 0.8f), 10.f));
            ent.add_components<transform>(position(10, 10, 10), rotation::lookat(math::vec3(1, 1, 1), math::vec3::zero), scale());
        }

        {
            auto ent = createEntity();
            ent.add_components<rendering::renderable>(spotLightH.get_mesh(), rendering::mesh_renderer(spotLightMH));
            ent.add_component<rendering::light>(rendering::light::spot(math::colors::green, math::deg2rad(45.f), 1.f, 100.f));
            ent.add_components<transform>(position(-10, 0.5, -10), rotation::lookat(math::vec3(0, 0, -1), math::vec3::zero), scale());
        }

        {
            auto ent = createEntity();
            ent.add_components<rendering::renderable>(spotLightH.get_mesh(), rendering::mesh_renderer(spotLightMH));
            ent.add_component<rendering::light>(rendering::light::spot(math::colors::green, math::deg2rad(45.f), 1.f, 100.f));
            ent.add_components<transform>(position(0, 0.5, -10), rotation::lookat(math::vec3(0, 0, -1), math::vec3::zero), scale());
        }

        {
            auto ent = createEntity();
            ent.add_components<rendering::renderable>(spotLightH.get_mesh(), rendering::mesh_renderer(spotLightMH));
            ent.add_component<rendering::light>(rendering::light::spot(math::colors::green, math::deg2rad(45.f), 1.f, 100.f));
            ent.add_components<transform>(position(10, 0.5, -10), rotation::lookat(math::vec3(0, 0, -1), math::vec3::zero), scale());
        }

        {
            auto ent = createEntity();
            ent.add_components<rendering::renderable>(pointLightH.get_mesh(), rendering::mesh_renderer(pointLightMH));
            ent.add_component<rendering::light>(rendering::light::point(math::colors::red, 1.f));
            ent.add_components<transform>(position(0, 1, 0), rotation(), scale());
        }


        {
            auto ent = createEntity();
            ent.add_components<rendering::renderable>(pointLightH.get_mesh(), rendering::mesh_renderer(pointLightMH));
            ent.add_component<rendering::light>(rendering::light::point(math::colors::red, 1.f));
            ent.add_components<transform>(position(-10, 1, 0), rotation(), scale());
        }


        {
            auto ent = createEntity();
            ent.add_components<rendering::renderable>(pointLightH.get_mesh(), rendering::mesh_renderer(pointLightMH));
            ent.add_component<rendering::light>(rendering::light::point(math::colors::red, 1.f));
            ent.add_components<transform>(position(10, 1, 0), rotation(), scale());
        }

        {
            auto ent = createEntity();
            ent.add_components<rendering::renderable>(suzanneH.get_mesh(), rendering::mesh_renderer(normalH));
            //ent.add_component<sah>({});
            ent.add_components<transform>(position(0, 3, 5.1f), rotation(), scale());
        }

        {
            auto ent = createEntity();
            ent.add_components<rendering::renderable>(suzanneH.get_mesh(), rendering::mesh_renderer(wireframeH));
            //ent.add_component<sah>({});
            ent.add_components<transform>(position(0, 3, 8.1f), rotation(), scale());
        }

        {
            auto ent = m_ecs->createEntity();
            ent.add_components<rendering::renderable>(suzanneH.get_mesh(), rendering::mesh_renderer(copperH));
            //ent.add_component<sah>({});

            ent.add_components<transform>(position(0, 3, 11.1f), rotation(), scale());
        }

        {
            auto ent = m_ecs->createEntity();
            ent.add_components<rendering::renderable>(gnomeH.get_mesh(), rendering::mesh_renderer(gnomeMH));
            //ent.add_component<sah>({});

            ent.add_components<transform>(position(0, 3, 2.1f), rotation(), scale());
        }

        /*   {
               auto ent = m_ecs->createEntity();
               ent.add_component<sah>();
               m_ecs->createComponent<rendering::renderable>(ent, { gnomeH, normalH });

               auto [positionH, rotationH, scaleH] = m_ecs->createComponent<transform>(ent);
               positionH.write(math::vec3(0, 2, 15.1f));
               scaleH.write(math::vec3(1.f));
           }*/

        {
            auto ent = createEntity();
            ent.add_components<rendering::renderable>(submeshtestH.get_mesh(), rendering::mesh_renderer(pbrH));
            //ent.add_component<sah>({});
            ent.add_components<transform>(position(0, 10, 0), rotation(), scale());
        }

        {
            auto ent = createEntity();
            ent.add_components<rendering::renderable>(axesH.get_mesh(), rendering::mesh_renderer(vertexColorH));
            ent.add_components<transform>();
        }

        {
            auto ent = createEntity();
            ent.add_components<rendering::renderable>(cubeH.get_mesh(), rendering::mesh_renderer(pbrH));
            //ent.add_component<sah>({});
            ent.add_components<transform>(position(5.1f, 9, 0), rotation(), scale());
        }

        {
            auto ent = createEntity();
            ent.add_components<rendering::renderable>(sphereH.get_mesh(), rendering::mesh_renderer(copperH));
            //ent.add_component<sah>({});
            ent.add_components<transform>(position(0, 3, -5.1f), rotation(), scale());
        }

        {
            auto ent = createEntity();
            ent.add_components<rendering::renderable>(sphereH.get_mesh(), rendering::mesh_renderer(aluminumH));
            //ent.add_component<sah>({});
            ent.add_components<transform>(position(0, 3, -8.f), rotation(), scale());
        }

        {
            auto ent = createEntity();
            ent.add_components<rendering::renderable>(sphereH.get_mesh(), rendering::mesh_renderer(ironH));
            //ent.add_component<sah>({});
            ent.add_components<transform>(position(0, 3, -2.2f), rotation(), scale());
        }

        {
            auto ent = createEntity();
            ent.add_components<rendering::renderable>(sphereH.get_mesh(), rendering::mesh_renderer(rock2H));
            //ent.add_component<sah>({});
            ent.add_components<transform>(position(4, 3, -8.f), rotation(), scale());
        }

        {
            auto ent = createEntity();
            ent.add_components<rendering::renderable>(sphereH.get_mesh(), rendering::mesh_renderer(fabricH));
            //ent.add_component<sah>({});
            ent.add_components<transform>(position(4, 3, -5.1f), rotation(), scale());
        }

        {
            auto ent = createEntity();
            ent.add_components<rendering::renderable>(sphereH.get_mesh(), rendering::mesh_renderer(paintH));
           //ent.add_component<sah>({});
            ent.add_components<transform>(position(4, 3, -2.2f), rotation(), scale());
        }

        {
            auto ent = createEntity();
            ent.add_components<rendering::renderable>(uvsphereH.get_mesh(), rendering::mesh_renderer(copperH));
            //ent.add_component<sah>({});
            ent.add_components<transform>(position(0, 3, -3.6f), rotation(), scale());
        }

        {
            auto ent = createEntity();
            ent.add_components<rendering::renderable>(uvsphereH.get_mesh(), rendering::mesh_renderer(aluminumH));
           //ent.add_component<sah>({});
            ent.add_components<transform>(position(0, 3, -6.5f), rotation(), scale());
        }

        {
            auto ent = createEntity();
            ent.add_components<rendering::renderable>(uvsphereH.get_mesh(), rendering::mesh_renderer(ironH));
            //ent.add_component<sah>({});
            ent.add_components<transform>(position(0, 3, -0.7f), rotation(), scale());
        }

        {
            auto ent = createEntity();
            ent.add_components<rendering::renderable>(uvsphereH.get_mesh(), rendering::mesh_renderer(rock2H));
            //ent.add_component<sah>({});
            ent.add_components<transform>(position(4, 3, -6.5f), rotation(), scale());
        }

        {
            auto ent = createEntity();
            ent.add_components<rendering::renderable>(uvsphereH.get_mesh(), rendering::mesh_renderer(fabricH));
            //ent.add_component<sah>({});
            ent.add_components<transform>(position(4, 3, -3.6f), rotation(), scale());
        }

        {
            auto ent = createEntity();
            ent.add_components<rendering::renderable>(uvsphereH.get_mesh(), rendering::mesh_renderer(paintH));
            //ent.add_component<sah>({});
            ent.add_components<transform>(position(4, 3, -0.7f), rotation(), scale());
            auto ent2 = ent.clone();
            auto pos = ent2.get_component_handle<position>().read();
            pos.y = 6;
            ent2.get_component_handle<position>().write(pos);
        }

        //audioSphereLeft setup
        {
            audioSphereLeft = createEntity();
            audioSphereLeft.add_components<rendering::renderable>(audioSourceH.get_mesh(), rendering::mesh_renderer(gizmoMH));
            audioSphereLeft.add_components<transform>(position(-5, 1, 10), rotation(), scale(0.5));

            auto segment = audio::AudioSegmentCache::createAudioSegment("kilogram", "assets://audio/kilogram-of-scotland_stereo32.wav"_view, { audio::audio_import_settings::channel_processing_setting::split_channels });

            audio::audio_source source;
            source.setAudioHandle(segment);
            audioSphereLeft.add_component<audio::audio_source>(source);
        }
        //audioSphereRight setup
        {
            audioSphereRight = createEntity();
            audioSphereRight.add_components<rendering::renderable>(audioSourceH.get_mesh(), rendering::mesh_renderer(gizmoMH));
            audioSphereRight.add_components<transform>(position(5, 1, 10), rotation(), scale(0.5));

            auto segment = audio::AudioSegmentCache::getAudioSegment("kilogram_channel1");

            audio::audio_source source;
            source.setAudioHandle(segment);
            audioSphereRight.add_component<audio::audio_source>(source);
        }
#pragma endregion



        //---------------------------------------------------------- Physics Collision Unit Test -------------------------------------------------------------------//

        //setupPhysicsCDUnitTest(cubeH, wireframeH);

        //----------- Rigidbody-Collider AABB Test------------//


        //setupPhysicsFrictionUnitTest(cubeH, uvH);

        //setupPhysicsStackingUnitTest(cubeH, uvH);

        physics::cube_collider_params cubeParams;
        cubeParams.breadth = 1.0f;
        cubeParams.width = 1.0f;
        cubeParams.height = 1.0f;
        //setupPhysicsCRUnitTest(cubeH, uvH);

        //CreateCubeStack(3, 2, 2, math::vec3(0, -3.0f, 8.0f), math::vec3(1, 1, 1)
        //    ,cubeParams, 0.1f, cubeH, wireframeH);

        createProcess<&TestSystem::update>("Update");
    }


  /*  void testPhysicsEvent(physics::TriggerEvent* evnt)
    {
        log::debug("received trigger event {}", evnt->manifold.isColliding);
    }*/

    void setupPhysicsCDUnitTest(rendering::model_handle cubeH, rendering::material_handle wireframeH)
    {
        physics::cube_collider_params cubeParams;
        cubeParams.breadth = 2.0f;
        cubeParams.width = 2.0f;
        cubeParams.height = 2.0f;

        physics::cube_collider_params cubeParams2;
        cubeParams2.breadth = 3.0f;
        cubeParams2.width = 3.0f;
        cubeParams2.height = 3.0f;

        //----------- AABB to AABB Test(Single Axis)  ------------//
        //**
        {
            auto ent = m_ecs->createEntity();

            auto entPhyHande = ent.add_component<physics::physicsComponent>();

            physics::physicsComponent physicsComponent2;
            physics::physicsComponent::init(physicsComponent2);


            physicsComponent2.AddBox(cubeParams);
            physicsComponent2.isTrigger = false;
            entPhyHande.write(physicsComponent2);

            ent.add_components<rendering::renderable>(cubeH.get_mesh(), rendering::mesh_renderer(wireframeH));

            auto [positionH, rotationH, scaleH] = m_ecs->createComponents<transform>(ent);
            positionH.write(math::vec3(0, -3.0f, 8.0f));
            scaleH.write(math::vec3(1.0f));
        }

        {
            auto ent = m_ecs->createEntity();
            physicsUnitTestCD.push_back(ent);
            auto entPhyHande = ent.add_component<physics::physicsComponent>();

            physics::physicsComponent physicsComponent2;
            physics::physicsComponent::init(physicsComponent2);


            physicsComponent2.AddBox(cubeParams);
            physicsComponent2.isTrigger = true;
            entPhyHande.write(physicsComponent2);


            ent.add_components<rendering::renderable>(cubeH.get_mesh(), rendering::mesh_renderer(wireframeH));

            auto [positionH, rotationH, scaleH] = m_ecs->createComponents<transform>(ent);
            positionH.write(math::vec3(3.0, -3.0f, 8.0f));
            scaleH.write(math::vec3(1.0f));
        }



        //----------- AABB to AABB Test  ------------//
        //**
        {
            auto ent = m_ecs->createEntity();

            auto entPhyHande = ent.add_component<physics::physicsComponent>();

            physics::physicsComponent physicsComponent2;
            physics::physicsComponent::init(physicsComponent2);


            physicsComponent2.AddBox(cubeParams);
            physicsComponent2.isTrigger = false;
            entPhyHande.write(physicsComponent2);


            ent.add_components<rendering::renderable>(cubeH.get_mesh(), rendering::mesh_renderer(wireframeH));

            auto [positionH, rotationH, scaleH] = m_ecs->createComponents<transform>(ent);
            positionH.write(math::vec3(0, -3.0f, 5.0f));
            scaleH.write(math::vec3(1.0f));
        }

        {
            auto ent = m_ecs->createEntity();
            physicsUnitTestCD.push_back(ent);
            auto entPhyHande = ent.add_component<physics::physicsComponent>();

            physics::physicsComponent physicsComponent2;
            physics::physicsComponent::init(physicsComponent2);


            physicsComponent2.AddBox(cubeParams);
            physicsComponent2.isTrigger = true;
            entPhyHande.write(physicsComponent2);


            ent.add_components<rendering::renderable>(cubeH.get_mesh(), rendering::mesh_renderer(wireframeH));

            auto [positionH, rotationH, scaleH] = m_ecs->createComponents<transform>(ent);
            positionH.write(math::vec3(3.0, -2.0f, 4.0f));
            scaleH.write(math::vec3(1.0f));
        }
        //*/

        //----------- AABB to OBB Test  ------------//
        //**
        {
            auto ent = m_ecs->createEntity();

            auto entPhyHande = ent.add_component<physics::physicsComponent>();

            physics::physicsComponent physicsComponent2;
            physics::physicsComponent::init(physicsComponent2);


            physicsComponent2.AddBox(cubeParams);
            physicsComponent2.isTrigger = false;
            entPhyHande.write(physicsComponent2);


            ent.add_components<rendering::renderable>(cubeH.get_mesh(), rendering::mesh_renderer(wireframeH));

            auto [positionH, rotationH, scaleH] = m_ecs->createComponents<transform>(ent);
            positionH.write(math::vec3(0, -3.0f, -2.0f));
            scaleH.write(math::vec3(1.0f));
        }

        {
            auto ent = m_ecs->createEntity();
            physicsUnitTestCD.push_back(ent);
            auto entPhyHande = ent.add_component<physics::physicsComponent>();

            physics::physicsComponent physicsComponent2;
            physics::physicsComponent::init(physicsComponent2);


            physicsComponent2.AddBox(cubeParams);
            physicsComponent2.isTrigger = true;
            entPhyHande.write(physicsComponent2);


            ent.add_components<rendering::renderable>(cubeH.get_mesh(), rendering::mesh_renderer(wireframeH));

            auto [positionH, rotationH, scaleH] = m_ecs->createComponents<transform>(ent);
            positionH.write(math::vec3(3.0, -3.0f, -2.0f));

            auto rot = rotationH.read();
            rot *= math::angleAxis(45.f, math::vec3(0, 1, 0));
            rot *= math::angleAxis(45.f, math::vec3(0, 0, 1));

            rotationH.write(rot);

            scaleH.write(math::vec3(1.0f));
        }
        //*/

        //----------- OBB to OBB Test  ------------//
        //**
        {
            auto ent = m_ecs->createEntity();

            auto entPhyHande = ent.add_component<physics::physicsComponent>();

            physics::physicsComponent physicsComponent2;
            physics::physicsComponent::init(physicsComponent2);


            physicsComponent2.AddBox(cubeParams);
            physicsComponent2.isTrigger = false;
            entPhyHande.write(physicsComponent2);


            ent.add_components<rendering::renderable>(cubeH.get_mesh(), rendering::mesh_renderer(wireframeH));

            auto [positionH, rotationH, scaleH] = m_ecs->createComponents<transform>(ent);
            positionH.write(math::vec3(0, -3.0f, -7.0f));

            auto rot = rotationH.read();
            rot *= math::angleAxis(20.f, math::vec3(0, 1, 0));
            rotationH.write(rot);

            scaleH.write(math::vec3(1.0f));
        }

        {
            auto ent = m_ecs->createEntity();
            physicsUnitTestCD.push_back(ent);
            auto entPhyHande = ent.add_component<physics::physicsComponent>();

            physics::physicsComponent physicsComponent2;
            physics::physicsComponent::init(physicsComponent2);


            physicsComponent2.AddBox(cubeParams);
            physicsComponent2.isTrigger = true;
            entPhyHande.write(physicsComponent2);

            ent.add_components<rendering::renderable>(cubeH.get_mesh(), rendering::mesh_renderer(wireframeH));

            auto [positionH, rotationH, scaleH] = m_ecs->createComponents<transform>(ent);
            positionH.write(math::vec3(3.0, -3.0f, -7.0f));

            auto rot = rotationH.read();
            rot *= math::angleAxis(45.f, math::vec3(0, 1, 0));
            rot *= math::angleAxis(45.f, math::vec3(0, 0, 1));

            rotationH.write(rot);

            scaleH.write(math::vec3(1.0f));
        }
        //*/

        //----------- OBB Edge-Edge Test  ------------//
        //**
        {
            auto ent = m_ecs->createEntity();
            //physicsUnitTestObjects.push_back(ent);
            auto entPhyHande = ent.add_component<physics::physicsComponent>();

            physics::physicsComponent physicsComponent2;
            physics::physicsComponent::init(physicsComponent2);


            physicsComponent2.AddBox(cubeParams);
            physicsComponent2.isTrigger = false;
            entPhyHande.write(physicsComponent2);


            ent.add_components<rendering::renderable>(cubeH.get_mesh(), rendering::mesh_renderer(wireframeH));

            auto [positionH, rotationH, scaleH] = m_ecs->createComponents<transform>(ent);
            positionH.write(math::vec3(0, -3.0f, -12.0f));

            auto rot = rotationH.read();
            rot *= math::angleAxis(45.f, math::vec3(0, 1, 0));
            rotationH.write(rot);

            scaleH.write(math::vec3(1.0f));
        }
        {
            auto ent = m_ecs->createEntity();
            physicsUnitTestCD.push_back(ent);
            auto entPhyHande = ent.add_component<physics::physicsComponent>();

            physics::physicsComponent physicsComponent2;
            physics::physicsComponent::init(physicsComponent2);


            physicsComponent2.AddBox(cubeParams);
            physicsComponent2.isTrigger = true;
            entPhyHande.write(physicsComponent2);


            ent.add_components<rendering::renderable>(cubeH.get_mesh(), rendering::mesh_renderer(wireframeH));

            auto [positionH, rotationH, scaleH] = m_ecs->createComponents<transform>(ent);
            positionH.write(math::vec3(3.0, -3.0f, -13.0f));

            auto rot = rotationH.read();
            rot *= math::angleAxis(45.f, math::vec3(1, 0, 0));
            rot *= math::angleAxis(45.f, math::vec3(0, 1, 0));
            rotationH.write(rot);

            scaleH.write(math::vec3(1.0f));
        }
        //*/
    }
    //10,0,15
    void setupPhysicsCRUnitTest(rendering::model_handle cubeH, rendering::material_handle wireframeH)
    {
        physics::cube_collider_params cubeParams;
        cubeParams.breadth = 1.0f;
        cubeParams.width = 1.0f;
        cubeParams.height = 1.0f;

        physics::cube_collider_params staticBlockParams;
        staticBlockParams.breadth = 5.0f;
        staticBlockParams.width = 5.0f;
        staticBlockParams.height = 2.0f;

        float testPos = 10.0f;
        //----------- Static Block To AABB Body Stability Test ------------//
        {
            auto ent = m_ecs->createEntity();

            auto [positionH, rotationH, scaleH] = m_ecs->createComponents<transform>(ent);
            positionH.write(math::vec3(testPos, -3.0f, 15.0f));
            scaleH.write(math::vec3(2.5f, 1.0f, 2.5f));

            ent.add_components<rendering::renderable>(cubeH.get_mesh(), rendering::mesh_renderer(wireframeH));
        }

        {
            auto ent = m_ecs->createEntity();

            auto entPhyHande = ent.add_component<physics::physicsComponent>();

            physics::physicsComponent physicsComponent2;
            physics::physicsComponent::init(physicsComponent2);


            physicsComponent2.AddBox(staticBlockParams);
            physicsComponent2.isTrigger = false;
            entPhyHande.write(physicsComponent2);

            //auto renderableHandle = m_ecs->createComponent<rendering::renderable>(ent);
            //renderableHandle.write({ cubeH, wireframeH });

            auto [positionH, rotationH, scaleH] = m_ecs->createComponents<transform>(ent);
            positionH.write(math::vec3(testPos, -3.0f, 15.0f));
            scaleH.write(math::vec3(1.0f));

            auto idHandle = m_ecs->createComponent<physics::identifier>(ent);
            auto idComp = idHandle.read();
            idComp.id = "AABBStaticStable";
            idHandle.write(idComp);


        }

        {
            staticToAABBEntLinear = m_ecs->createEntity();
            physicsUnitTestCD.push_back(staticToAABBEntLinear);
            auto entPhyHande = staticToAABBEntLinear.add_component<physics::physicsComponent>();

            physics::physicsComponent physicsComponent2;
            physics::physicsComponent::init(physicsComponent2);


            physicsComponent2.AddBox(cubeParams);
            entPhyHande.write(physicsComponent2);

            //auto crb = m_ecs->createComponent<physics::rigidbody>(staticToAABBEnt);
            //auto rbHandle = staticToAABBEnt.add_component<physics::rigidbody>();

            //auto renderableHandle = m_ecs->createComponent<rendering::renderable>(staticToAABBEnt);
            //renderableHandle.write({ cubeH, wireframeH });

            auto [positionH, rotationH, scaleH] = m_ecs->createComponents<transform>(staticToAABBEntLinear);
            positionH.write(math::vec3(testPos, -0.0f, 15.0f));
            scaleH.write(math::vec3(1.0f));

            auto idHandle = m_ecs->createComponent<physics::identifier>(staticToAABBEntLinear);
            auto id = idHandle.read();
            id.id = "AABBRbStable";
            idHandle.write(id);
        }

        //{
        //    auto ent = m_ecs->createEntity();
        //    //physicsUnitTestCD.push_back(staticToAABBEntLinear);
        //    auto entPhyHande = ent.add_component<physics::physicsComponent>();

        //    physics::physicsComponent physicsComponent2;
        //    physics::physicsComponent::init(physicsComponent2);


        //    physicsComponent2.AddBox(cubeParams);
        //    entPhyHande.write(physicsComponent2);

        //    auto crb = m_ecs->createComponent<physics::rigidbody>(ent);
        //    auto rbHandle = ent.add_component<physics::rigidbody>();

        //    //auto renderableHandle = m_ecs->createComponent<rendering::renderable>(staticToAABBEnt);
        //    //renderableHandle.write({ cubeH, wireframeH });

        //    auto [positionH, rotationH, scaleH] = m_ecs->createComponents<transform>(ent);
        //    positionH.write(math::vec3(testPos+0.5f, -1.0f, 15.0f));
        //    scaleH.write(math::vec3(1.0f));

        //    auto idHandle = m_ecs->createComponent<physics::identifier>(ent);
        //    auto id = idHandle.read();
        //    id.id = "AABBRbStable";
        //    idHandle.write(id);
        //}

        //----------- Static Block To AABB Body Rotation Test------------//

        {
            auto ent = m_ecs->createEntity();

            auto [positionH, rotationH, scaleH] = m_ecs->createComponents<transform>(ent);
            positionH.write(math::vec3(testPos, -3.0f, 8.0f));
            scaleH.write(math::vec3(2.5f, 1.0f, 2.5f));

            ent.add_components<rendering::renderable>(cubeH.get_mesh(), rendering::mesh_renderer(wireframeH));
        }

        {
            auto ent = m_ecs->createEntity();

            auto entPhyHande = ent.add_component<physics::physicsComponent>();

            physics::physicsComponent physicsComponent2;
            physics::physicsComponent::init(physicsComponent2);


            physicsComponent2.AddBox(staticBlockParams);
            physicsComponent2.isTrigger = false;
            entPhyHande.write(physicsComponent2);

            //auto renderableHandle = m_ecs->createComponent<rendering::renderable>(ent);
            //renderableHandle.write({ cubeH, wireframeH });

            auto [positionH, rotationH, scaleH] = m_ecs->createComponents<transform>(ent);
            positionH.write(math::vec3(testPos, -3.0f, 8.0f));
            scaleH.write(math::vec3(1.0f));

            auto idHandle = m_ecs->createComponent<physics::identifier>(ent);
            auto idComp = idHandle.read();
            idComp.id = "AABBStatic";
            idHandle.write(idComp);


        }

        {
            staticToAABBEntRotation = m_ecs->createEntity();
            physicsUnitTestCD.push_back(staticToAABBEntRotation);
            auto entPhyHande = staticToAABBEntRotation.add_component<physics::physicsComponent>();

            physics::physicsComponent physicsComponent2;
            physics::physicsComponent::init(physicsComponent2);

            physicsComponent2.AddBox(cubeParams);
            entPhyHande.write(physicsComponent2);

            //auto crb = m_ecs->createComponent<physics::rigidbody>(staticToAABBEnt);
            //auto rbHandle = staticToAABBEnt.add_component<physics::rigidbody>();

            //auto renderableHandle = m_ecs->createComponent<rendering::renderable>(staticToAABBEnt);
            //renderableHandle.write({ cubeH, wireframeH });

            auto [positionH, rotationH, scaleH] = m_ecs->createComponents<transform>(staticToAABBEntRotation);
            positionH.write(math::vec3(testPos - 2.7f, -0.0f, 8.0f));
            scaleH.write(math::vec3(1.0f));

            auto idHandle = m_ecs->createComponent<physics::identifier>(staticToAABBEntRotation);
            auto id = idHandle.read();
            id.id = "AABBRb";
            idHandle.write(id);
        }

        //----------- Static Block To Rotated Body ------------//

        {
            auto ent = m_ecs->createEntity();

            auto [positionH, rotationH, scaleH] = m_ecs->createComponents<transform>(ent);
            positionH.write(math::vec3(testPos, -3.0f, 2.0f));
            scaleH.write(math::vec3(2.5f, 1.0f, 2.5f));

            ent.add_components<rendering::renderable>(cubeH.get_mesh(), rendering::mesh_renderer(wireframeH));
        }

        {
            auto ent = m_ecs->createEntity();

            auto entPhyHande = ent.add_component<physics::physicsComponent>();

            physics::physicsComponent physicsComponent2;
            physics::physicsComponent::init(physicsComponent2);



            physicsComponent2.AddBox(staticBlockParams);
            physicsComponent2.isTrigger = false;
            entPhyHande.write(physicsComponent2);

            //auto renderableHandle = m_ecs->createComponent<rendering::renderable>(ent);
            //renderableHandle.write({ cubeH, wireframeH });

            auto [positionH, rotationH, scaleH] = m_ecs->createComponents<transform>(ent);
            positionH.write(math::vec3(testPos, -3.0f, 2.0f));
            scaleH.write(math::vec3(1.0f));

            auto idHandle = m_ecs->createComponent<physics::identifier>(ent);
            auto id = idHandle.read();
            id.id = "OBBStatic";
            idHandle.write(id);
        }

        {
            staticToOBBEnt = m_ecs->createEntity();

            auto entPhyHande = staticToOBBEnt.add_component<physics::physicsComponent>();

            physics::physicsComponent physicsComponent2;
            physics::physicsComponent::init(physicsComponent2);


            physicsComponent2.AddBox(cubeParams);
            physicsComponent2.isTrigger = false;
            entPhyHande.write(physicsComponent2);


            //auto renderableHandle = m_ecs->createComponent<rendering::renderable>(staticToOBBEnt);
            //renderableHandle.write({ cubeH, wireframeH });

            auto [positionH, rotationH, scaleH] = m_ecs->createComponents<transform>(staticToOBBEnt);
            positionH.write(math::vec3(testPos, -0.0f, 2.0f));
            scaleH.write(math::vec3(1.0f));

            auto rot = rotationH.read();
            rot *= math::angleAxis(math::radians(90.f), math::vec3(0, 0, 1));
            rot *= math::angleAxis(math::radians(40.f), math::vec3(1, 0, 0));
            rot *= math::angleAxis(math::radians(42.f), math::vec3(0, 1, 0));
            //rot *= math::angleAxis(45.f, math::vec3(0, 1, 0));
            rotationH.write(rot);

            auto idHandle = m_ecs->createComponent<physics::identifier>(staticToOBBEnt);
            auto id = idHandle.read();
            id.id = "OBBRb";
            idHandle.write(id);


            //log::debug("rb.angularVelocity {}", rb.angularVelocity);

        }

        //----------- Static Block To Rotated Body Edge ------------//

        {
            auto ent = m_ecs->createEntity();

            auto [positionH, rotationH, scaleH] = m_ecs->createComponents<transform>(ent);
            positionH.write(math::vec3(testPos, -3.0f, -4.0f));
            scaleH.write(math::vec3(2.5f, 1.0f, 2.5f));

            ent.add_components<rendering::renderable>(cubeH.get_mesh(), rendering::mesh_renderer(wireframeH));
        }

        {
            auto ent = m_ecs->createEntity();

            auto entPhyHande = ent.add_component<physics::physicsComponent>();

            physics::physicsComponent physicsComponent2;
            physics::physicsComponent::init(physicsComponent2);



            physicsComponent2.AddBox(staticBlockParams);
            physicsComponent2.isTrigger = false;
            entPhyHande.write(physicsComponent2);

            //auto renderableHandle = m_ecs->createComponent<rendering::renderable>(ent);
            //renderableHandle.write({ cubeH, wireframeH });

            auto [positionH, rotationH, scaleH] = m_ecs->createComponents<transform>(ent);
            positionH.write(math::vec3(testPos, -3.0f, -4.0f));
            scaleH.write(math::vec3(1.0f));

            auto idHandle = m_ecs->createComponent<physics::identifier>(ent);
            auto id = idHandle.read();
            id.id = "OBBFoundationStaticEdge";
            idHandle.write(id);
        }

        {
            staticToEdgeEnt = m_ecs->createEntity();

            auto entPhyHande = staticToEdgeEnt.add_component<physics::physicsComponent>();

            physics::physicsComponent physicsComponent2;
            physics::physicsComponent::init(physicsComponent2);


            physicsComponent2.AddBox(cubeParams);
            entPhyHande.write(physicsComponent2);


            //auto renderableHandle = m_ecs->createComponent<rendering::renderable>(staticToEdgeEnt);
            //renderableHandle.write({ cubeH, wireframeH });

            auto [positionH, rotationH, scaleH] = m_ecs->createComponents<transform>(staticToEdgeEnt);
            positionH.write(math::vec3(testPos + 3.0f, 2.0f, -4.0f));
            scaleH.write(math::vec3(1.0f));

            auto rot = rotationH.read();
            rot *= math::angleAxis(45.f, math::vec3(1, 0, 0));
            rot *= math::angleAxis(45.f, math::vec3(0, 1, 0));
            rotationH.write(rot);

            auto idHandle = m_ecs->createComponent<physics::identifier>(staticToEdgeEnt);
            auto id = idHandle.read();
            id.id = "OBBRb";
            idHandle.write(id);


            //log::debug("rb.angularVelocity {}", rb.angularVelocity);

        }


    }
    //15,0,15
    void setupPhysicsFrictionUnitTest(rendering::model_handle cubeH, rendering::material_handle wireframeH)
    {
        float testPos = 20.f;
        physics::cube_collider_params cubeParams;
        cubeParams.breadth = 1.0f;
        cubeParams.width = 1.0f;
        cubeParams.height = 1.0f;

        physics::cube_collider_params staticBlockParams;
        staticBlockParams.breadth = 5.0f;
        staticBlockParams.width = 5.0f;
        staticBlockParams.height = 2.0f;

        //NO Friction
        {
            auto ent = m_ecs->createEntity();
            physicsFrictionTestRotators.push_back(ent);
            auto [positionH, rotationH, scaleH] = m_ecs->createComponents<transform>(ent);
            positionH.write(math::vec3(testPos, -3.0f, 15.0f));
            scaleH.write(math::vec3(2.5f, 1.0f, 2.5f));

            ent.add_components<rendering::renderable>(cubeH.get_mesh(), rendering::mesh_renderer(wireframeH));
        }

        {
            auto ent = m_ecs->createEntity();
            physicsFrictionTestRotators.push_back(ent);
            auto entPhyHande = ent.add_component<physics::physicsComponent>();

            physics::physicsComponent physicsComponent2;
            physics::physicsComponent::init(physicsComponent2);


            physicsComponent2.AddBox(staticBlockParams);
            physicsComponent2.isTrigger = false;
            entPhyHande.write(physicsComponent2);

            //auto renderableHandle = m_ecs->createComponent<rendering::renderable>(ent);
            //renderableHandle.write({ cubeH, wireframeH });

            auto [positionH, rotationH, scaleH] = m_ecs->createComponents<transform>(ent);
            positionH.write(math::vec3(testPos, -3.0f, 15.0f));
            scaleH.write(math::vec3(1.0f));

            auto idHandle = m_ecs->createComponent<physics::identifier>(ent);
            auto idComp = idHandle.read();
            idComp.id = "StaticNoFriction";
            idHandle.write(idComp);


        }

        {
            NoFrictionBody = m_ecs->createEntity();

            auto entPhyHande = NoFrictionBody.add_component<physics::physicsComponent>();

            physics::physicsComponent physicsComponent2;
            physics::physicsComponent::init(physicsComponent2);

            physicsComponent2.AddBox(cubeParams);
            physicsComponent2.isTrigger = false;
            entPhyHande.write(physicsComponent2);

            //auto renderableHandle = m_ecs->createComponent<rendering::renderable>(ent);
            //renderableHandle.write({ cubeH, wireframeH });

            auto [positionH, rotationH, scaleH] = m_ecs->createComponents<transform>(NoFrictionBody);
            positionH.write(math::vec3(testPos - 1.0f, -1.5f, 15.0f));
            scaleH.write(math::vec3(1.0f));

            auto idHandle = m_ecs->createComponent<physics::identifier>(NoFrictionBody);
            auto idComp = idHandle.read();
            idComp.id = "NoFriction";
            idHandle.write(idComp);

        }

        //0.3f Friction

        {
            auto ent = m_ecs->createEntity();
            physicsFrictionTestRotators.push_back(ent);
            auto [positionH, rotationH, scaleH] = m_ecs->createComponents<transform>(ent);
            positionH.write(math::vec3(testPos, -3.0f, 8.0f));
            scaleH.write(math::vec3(2.5f, 1.0f, 2.5f));

            ent.add_components<rendering::renderable>(cubeH.get_mesh(), rendering::mesh_renderer(wireframeH));
        }

        {
            auto ent = m_ecs->createEntity();
            physicsFrictionTestRotators.push_back(ent);
            auto entPhyHande = ent.add_component<physics::physicsComponent>();

            physics::physicsComponent physicsComponent2;
            physics::physicsComponent::init(physicsComponent2);


            physicsComponent2.AddBox(staticBlockParams);
            physicsComponent2.isTrigger = false;
            entPhyHande.write(physicsComponent2);

            //auto renderableHandle = m_ecs->createComponent<rendering::renderable>(ent);
            //renderableHandle.write({ cubeH, wireframeH });

            auto [positionH, rotationH, scaleH] = m_ecs->createComponents<transform>(ent);
            positionH.write(math::vec3(testPos, -3.0f, 8.0f));
            scaleH.write(math::vec3(1.0f));

            auto idHandle = m_ecs->createComponent<physics::identifier>(ent);
            auto idComp = idHandle.read();
            idComp.id = "StaticNoFriction";
            idHandle.write(idComp);


        }

        {
            Point3FrictionBody = m_ecs->createEntity();

            auto entPhyHande = Point3FrictionBody.add_component<physics::physicsComponent>();

            physics::physicsComponent physicsComponent2;
            physics::physicsComponent::init(physicsComponent2);


            physicsComponent2.AddBox(cubeParams);
            physicsComponent2.isTrigger = false;
            entPhyHande.write(physicsComponent2);

            //auto renderableHandle = m_ecs->createComponent<rendering::renderable>(ent);
            //renderableHandle.write({ cubeH, wireframeH });

            auto [positionH, rotationH, scaleH] = m_ecs->createComponents<transform>(Point3FrictionBody);
            positionH.write(math::vec3(testPos - 1.0f, -1.5f, 8.0f));
            scaleH.write(math::vec3(1.0f));

            auto idHandle = m_ecs->createComponent<physics::identifier>(Point3FrictionBody);
            auto idComp = idHandle.read();
            idComp.id = "Point3Friction";
            idHandle.write(idComp);

        }


        //0.6f Friction

        {
            auto ent = m_ecs->createEntity();
            physicsFrictionTestRotators.push_back(ent);
            auto [positionH, rotationH, scaleH] = m_ecs->createComponents<transform>(ent);
            positionH.write(math::vec3(testPos, -3.0f, 1.0f));
            scaleH.write(math::vec3(2.5f, 1.0f, 2.5f));

            ent.add_components<rendering::renderable>(cubeH.get_mesh(), rendering::mesh_renderer(wireframeH));
        }

        {
            auto ent = m_ecs->createEntity();
            physicsFrictionTestRotators.push_back(ent);
            auto entPhyHande = ent.add_component<physics::physicsComponent>();

            physics::physicsComponent physicsComponent2;
            physics::physicsComponent::init(physicsComponent2);


            physicsComponent2.AddBox(staticBlockParams);
            physicsComponent2.isTrigger = false;
            entPhyHande.write(physicsComponent2);

            //auto renderableHandle = m_ecs->createComponent<rendering::renderable>(ent);
            //renderableHandle.write({ cubeH, wireframeH });

            auto [positionH, rotationH, scaleH] = m_ecs->createComponents<transform>(ent);
            positionH.write(math::vec3(testPos, -3.0f, 1.0f));
            scaleH.write(math::vec3(1.0f));

            auto idHandle = m_ecs->createComponent<physics::identifier>(ent);
            auto idComp = idHandle.read();
            idComp.id = "StaticNoFriction";
            idHandle.write(idComp);


        }

        {
            Point6FrictionBody = m_ecs->createEntity();

            auto entPhyHande = Point6FrictionBody.add_component<physics::physicsComponent>();

            physics::physicsComponent physicsComponent2;
            physics::physicsComponent::init(physicsComponent2);


            physicsComponent2.AddBox(cubeParams);
            physicsComponent2.isTrigger = false;
            entPhyHande.write(physicsComponent2);

            //auto renderableHandle = m_ecs->createComponent<rendering::renderable>(ent);
            //renderableHandle.write({ cubeH, wireframeH });

            auto [positionH, rotationH, scaleH] = m_ecs->createComponents<transform>(Point6FrictionBody);
            positionH.write(math::vec3(testPos - 1.0f, -1.5f, 1.0f));
            scaleH.write(math::vec3(1.0f));

            auto idHandle = m_ecs->createComponent<physics::identifier>(Point6FrictionBody);
            auto idComp = idHandle.read();
            idComp.id = "NoFriction";
            idHandle.write(idComp);

        }

        //1.0f Friction

        {
            auto ent = m_ecs->createEntity();
            physicsFrictionTestRotators.push_back(ent);
            auto [positionH, rotationH, scaleH] = m_ecs->createComponents<transform>(ent);
            positionH.write(math::vec3(testPos, -3.0f, -6.0f));
            scaleH.write(math::vec3(2.5f, 1.0f, 2.5f));

            ent.add_components<rendering::renderable>(cubeH.get_mesh(), rendering::mesh_renderer(wireframeH));
        }

        {
            auto ent = m_ecs->createEntity();
            physicsFrictionTestRotators.push_back(ent);
            auto entPhyHande = ent.add_component<physics::physicsComponent>();

            physics::physicsComponent physicsComponent2;
            physics::physicsComponent::init(physicsComponent2);


            physicsComponent2.AddBox(staticBlockParams);
            physicsComponent2.isTrigger = false;
            entPhyHande.write(physicsComponent2);

            //auto renderableHandle = m_ecs->createComponent<rendering::renderable>(ent);
            //renderableHandle.write({ cubeH, wireframeH });

            auto [positionH, rotationH, scaleH] = m_ecs->createComponents<transform>(ent);
            positionH.write(math::vec3(testPos, -3.0f, -6.0f));
            scaleH.write(math::vec3(1.0f));

            auto idHandle = m_ecs->createComponent<physics::identifier>(ent);
            auto idComp = idHandle.read();
            idComp.id = "StaticNoFriction";
            idHandle.write(idComp);


        }

        {
            FullFrictionBody = m_ecs->createEntity();

            auto entPhyHande = FullFrictionBody.add_component<physics::physicsComponent>();

            physics::physicsComponent physicsComponent2;
            physics::physicsComponent::init(physicsComponent2);


            physicsComponent2.AddBox(cubeParams);
            physicsComponent2.isTrigger = false;
            entPhyHande.write(physicsComponent2);

            //auto renderableHandle = m_ecs->createComponent<rendering::renderable>(ent);
            //renderableHandle.write({ cubeH, wireframeH });

            auto [positionH, rotationH, scaleH] = m_ecs->createComponents<transform>(FullFrictionBody);
            positionH.write(math::vec3(testPos - 1.0f, -1.5f, -6.0f));
            scaleH.write(math::vec3(1.0f));

            auto idHandle = m_ecs->createComponent<physics::identifier>(FullFrictionBody);
            auto idComp = idHandle.read();
            idComp.id = "FullFrictionBody";
            idHandle.write(idComp);

        }

        //


    }
    //20,0,15
    void setupPhysicsStackingUnitTest(rendering::model_handle cubeH, rendering::material_handle wireframeH)
    {
        float testPos = 20.f;
        physics::cube_collider_params cubeParams;
        cubeParams.breadth = 1.0f;
        cubeParams.width = 1.0f;
        cubeParams.height = 1.0f;

        physics::cube_collider_params staticBlockParams;
        staticBlockParams.breadth = 5.0f;
        staticBlockParams.width = 5.0f;
        staticBlockParams.height = 2.0f;

        //2 stack
        {
            auto ent = m_ecs->createEntity();

            auto [positionH, rotationH, scaleH] = m_ecs->createComponents<transform>(ent);
            positionH.write(math::vec3(35, -3.0f, 15.0f));
            scaleH.write(math::vec3(2.5f, 1.0f, 2.5f));

            ent.add_components<rendering::renderable>(cubeH.get_mesh(), rendering::mesh_renderer(wireframeH));
        }

        {
            auto ent = m_ecs->createEntity();

            auto entPhyHande = ent.add_component<physics::physicsComponent>();

            physics::physicsComponent physicsComponent2;
            physics::physicsComponent::init(physicsComponent2);


            physicsComponent2.AddBox(staticBlockParams);
            physicsComponent2.isTrigger = false;
            entPhyHande.write(physicsComponent2);

            //auto renderableHandle = m_ecs->createComponent<rendering::renderable>(ent);
            //renderableHandle.write({ cubeH, wireframeH });

            auto [positionH, rotationH, scaleH] = m_ecs->createComponents<transform>(ent);
            positionH.write(math::vec3(35, -3.0f, 15.0f));
            scaleH.write(math::vec3(1.0f));

            auto idHandle = m_ecs->createComponent<physics::identifier>(ent);
            auto idComp = idHandle.read();
            idComp.id = "Stack";
            idHandle.write(idComp);


        }

        CreateCubeStack(2, 1, 1, math::vec3(35, -1.5f, 15.0f), math::vec3(1, 1, 1)
            , cubeParams, 0.1f, cubeH, wireframeH);

        //10 stactk

        //Pyramid Stack
    }

#pragma region input stuff

    void saveScene1(savescene1* action)
    {
        if (action->pressed())
        {
            log::debug("Saving scene: Main");
            scenemanagement::SceneManager::createScene("Main");
            log::debug("Finishded saving scene: Main");
        }
    }
    void saveScene2(savescene2* action)
    {
        if (action->pressed())
        {
            log::debug("Saving scene: Main2");
            scenemanagement::SceneManager::createScene("Main2");
            log::debug("Finishded saving scene: Main2");
        }
    }
    void loadScene1(loadscene1* action)
    {
        if (action->pressed())
        {
            log::debug("Started loading scene: Main");
            scenemanagement::SceneManager::loadScene("Main");
            log::debug("Finished loading a scene");
        }
    }
    void loadScene2(loadscene2* action)
    {
        if (action->pressed())
        {
            log::debug("Started loading scene: Main");
            scenemanagement::SceneManager::loadScene("Main2");
            log::debug("Finished loading a scene");
        }
    }

    void onLightSwitch(light_switch* action)
    {
        static bool on = true;

        if (!action->value)
        {
            //auto light = sun.read_component<rendering::light>();
            if (on)
            {
                /*light.set_intensity(0.f);
                sun.write_component(light);*/

                if (sun)
                    sun.destroy();

                pbrH.set_param("skycolor", math::color(0.002f, 0.01f, 0.05f));
                copperH.set_param("skycolor", math::color(0.002f, 0.01f, 0.05f));
                aluminumH.set_param("skycolor", math::color(0.002f, 0.01f, 0.05f));
                ironH.set_param("skycolor", math::color(0.002f, 0.01f, 0.05f));
                slateH.set_param("skycolor", math::color(0.002f, 0.01f, 0.05f));
                rockH.set_param("skycolor", math::color(0.002f, 0.01f, 0.05f));
                rock2H.set_param("skycolor", math::color(0.002f, 0.01f, 0.05f));
                fabricH.set_param("skycolor", math::color(0.002f, 0.01f, 0.05f));
                bogH.set_param("skycolor", math::color(0.002f, 0.01f, 0.05f));
                paintH.set_param("skycolor", math::color(0.002f, 0.01f, 0.05f));
                skyboxH.set_param("skycolor", math::color(0.002f, 0.01f, 0.05f));
            }
            else
            {
                if (!sun)
                {
                    sun = createEntity();
                    sun.add_components<rendering::renderable>(MeshCache::get_handle("directional light"), rendering::mesh_renderer(rendering::MaterialCache::get_material("directional light")));
                    sun.add_component<rendering::light>(rendering::light::directional(math::color(1, 1, 0.8f), 10.f));
                    sun.add_components<transform>(position(10, 10, 10), rotation::lookat(math::vec3(1, 1, 1), math::vec3::zero), scale());
                }

                pbrH.set_param("skycolor", math::color(0.2f, 0.4f, 1.0f));
                copperH.set_param("skycolor", math::color(0.2f, 0.4f, 1.0f));
                aluminumH.set_param("skycolor", math::color(0.2f, 0.4f, 1.0f));
                ironH.set_param("skycolor", math::color(0.2f, 0.4f, 1.0f));
                slateH.set_param("skycolor", math::color(0.2f, 0.4f, 1.0f));
                rockH.set_param("skycolor", math::color(0.2f, 0.4f, 1.0f));
                rock2H.set_param("skycolor", math::color(0.2f, 0.4f, 1.0f));
                fabricH.set_param("skycolor", math::color(0.2f, 0.4f, 1.0f));
                bogH.set_param("skycolor", math::color(0.2f, 0.4f, 1.0f));
                paintH.set_param("skycolor", math::color(0.2f, 0.4f, 1.0f));
                skyboxH.set_param("skycolor", math::color(0.2f, 0.4f, 1.0f));
            }
            on = !on;
        }
    }

    void onTonemapSwitch(tonemap_switch* action)
    {
        static bool on = false;

        if (!action->value)
        {
            if (on)
            {
                pbrH.set_param("tonemap", false);
                copperH.set_param("tonemap", false);
                aluminumH.set_param("tonemap", false);
                ironH.set_param("tonemap", false);
                slateH.set_param("tonemap", false);
                rockH.set_param("tonemap", false);
                rock2H.set_param("tonemap", false);
                fabricH.set_param("tonemap", false);
                bogH.set_param("tonemap", false);
                paintH.set_param("tonemap", false);
                skyboxH.set_param("tonemap", false);
            }
            else
            {
                pbrH.set_param("tonemap", true);
                copperH.set_param("tonemap", true);
                aluminumH.set_param("tonemap", true);
                ironH.set_param("tonemap", true);
                slateH.set_param("tonemap", true);
                rockH.set_param("tonemap", true);
                rock2H.set_param("tonemap", true);
                fabricH.set_param("tonemap", true);
                bogH.set_param("tonemap", true);
                paintH.set_param("tonemap", true);
                skyboxH.set_param("tonemap", true);
            }
            on = !on;
        }
    }

    void onSphereAAMove(audio_move* action)
    {
        auto posH = audioSphereLeft.get_component_handle<position>();
        math::vec3 move = math::vec3(0.f, 0.f, 1.f);
        move = math::normalize(move * math::vec3(1, 0, 1)) * action->value * action->input_delta * 10.f;
        posH.fetch_add(move);
    }

    void onSphereAAStrive(audio_strive* action)
    {
        auto posH = audioSphereLeft.get_component_handle<position>();
        math::vec3 move = math::vec3(1.f, 0.f, 0.f);
        move = math::normalize(move * math::vec3(1, 0, 1)) * action->value * action->input_delta * 10.f;
        posH.fetch_add(move);
    }

    void onGainChange(gain_change* action)
    {
        if (action->value == 0) return;
        using namespace audio;
        auto sourceH = audioSphereLeft.get_component_handle<audio_source>();
        audio_source source = sourceH.read();
        float g = source.getGain() + action->value * action->input_delta * 10.0f;
        source.setGain(g);
        sourceH.write(source);

        sourceH = audioSphereRight.get_component_handle<audio_source>();
        source = sourceH.read();
        g = source.getGain() + action->value * action->input_delta * 10.0f;
        source.setGain(g);
        sourceH.write(source);
    }

    void onPitchChange(pitch_change* action)
    {
        if (action->value == 0) return;
        using namespace audio;
        auto sourceH = audioSphereLeft.get_component_handle<audio_source>();
        audio_source source = sourceH.read();
        const float p = source.getPitch() + action->value * action->input_delta * 10.0f;
        source.setPitch(p);
        sourceH.write(source);
    }

    void audioTestInput(audio_test_input* action)
    {
        using namespace audio;
        auto sourceH = audioSphereLeft.get_component_handle<audio_source>();
        audio_source source = sourceH.read();
        source.disableSpatialAudio();
        sourceH.write(source);
    }

    void onUnitPhysicsUnitTestMove(physics_test_move* action)
    {
        for (auto entity : physicsUnitTestCD)
        {
            auto posHandle = entity.get_component_handle<position>();

            auto pos = posHandle.read();
            pos.x += 0.005 * action->value;


            posHandle.write(pos);

        }
    }

    void playAudioSource(play_audio_source* action)
    {
        {
            auto sourceH = audioSphereLeft.get_component_handle<audio::audio_source>();
            audio::audio_source source = sourceH.read();
            source.play();
            sourceH.write(source);
        }
        {
            auto sourceH = audioSphereRight.get_component_handle<audio::audio_source>();
            audio::audio_source source = sourceH.read();
            source.play();
            sourceH.write(source);
        }
    }

    void pauseAudioSource(pause_audio_source* action)
    {
        {
            auto sourceH = audioSphereLeft.get_component_handle<audio::audio_source>();
            audio::audio_source source = sourceH.read();
            source.pause();
            sourceH.write(source);
        }
        {
            auto sourceH = audioSphereRight.get_component_handle<audio::audio_source>();
            audio::audio_source source = sourceH.read();
            source.pause();
            sourceH.write(source);
        }
    }

    void stopAudioSource(stop_audio_source* action)
    {
        {
            auto sourceH = audioSphereLeft.get_component_handle<audio::audio_source>();
            audio::audio_source source = sourceH.read();
            source.stop();
            sourceH.write(source);
        }
        {
            auto sourceH = audioSphereRight.get_component_handle<audio::audio_source>();
            audio::audio_source source = sourceH.read();
            source.stop();
            sourceH.write(source);
        }
    }

    void rewindAudioSource(rewind_audio_source* action)
    {
        auto sourceH = audioSphereLeft.get_component_handle<audio::audio_source>();
        audio::audio_source source = sourceH.read();
        source.rewind();
        sourceH.write(source);
    }
#pragma endregion

    void update(time::span deltaTime)
    {

        static auto sahQuery = createQuery<sah, rotation, position>();

        sahQuery.queryEntities();
        for (auto entity : sahQuery)
        {
            auto rot = entity.read_component<rotation>();

            rot *= math::angleAxis(math::deg2rad(45.f * deltaTime), math::vec3(0, 1, 0));

            entity.write_component(rot);

            auto pos = entity.read_component<position>();
            debug::drawLine(pos, pos + rot.forward(), math::colors::magenta);
        }

        if (rotate && !physics::PhysicsSystem::IsPaused)
        {
            for (auto entity : physicsFrictionTestRotators)
            {
                auto rot = entity.read_component<rotation>();

                rot *= math::angleAxis(math::deg2rad(-20.f * deltaTime), math::vec3(0, 0, 1));

                entity.write_component(rot);
            }
        }
    }

    void drawInterval(time::span deltaTime)
    {
        static auto physicsQuery = createQuery< physics::physicsComponent>();
        uint i = 0;

        float duration = 0.02f;

        for (auto penetration : physics::PhysicsSystem::penetrationQueries)
        {
            debug::drawLine(penetration->faceCentroid
                , penetration->faceCentroid + penetration->normal, math::vec4(1, 0, 1, 1), 15.0f, duration);
            auto x = 1;
        }


        //--------------------------------------- Draw contact points ---------------------------------------//



        for (int i = 0; i < physics::PhysicsSystem::contactPoints.size(); i++)
        {
            //ref is red
            //inc is blue

            auto& contact = physics::PhysicsSystem::contactPoints.at(i);

            debug::drawLine(contact.refRBCentroid
                , contact.RefWorldContact, math::vec4(1, 0, 0, 1), 5.0f, duration, true);

            debug::drawLine(contact.incRBCentroid
                , contact.IncWorldContact, math::vec4(0, 0, 1, 1), 5.0f, duration, true);

            debug::drawLine(contact.IncWorldContact
                , contact.IncWorldContact + math::vec3(0, 0.1f, 0), math::vec4(0.5, 0.5, 0.5, 1), 5.0f, duration, true);

            debug::drawLine(contact.refRBCentroid
                , contact.refRBCentroid + math::vec3(0, 0.1f, 0), math::vec4(0, 0, 0, 1), 5.0f, duration, true);

        }

        //--------------------------------------- Draw extreme points ---------------------------------------//

        i = 0;
        for (auto penetration : physics::PhysicsSystem::aPoint)
        {
            debug::drawLine(penetration
                , penetration + math::vec3(0, 0.2, 0), math::vec4(1, 0, 0, 1), 15.0f);

        }
        i = 0;
        for (auto penetration : physics::PhysicsSystem::bPoint)
        {
            debug::drawLine(penetration
                , penetration + math::vec3(0, 0.2, 0), math::vec4(0, 0, 1, 1), 15.0f);

        }

        physics::PhysicsSystem::contactPoints.clear();
        physics::PhysicsSystem::penetrationQueries.clear();
        physics::PhysicsSystem::aPoint.clear();
        physics::PhysicsSystem::bPoint.clear();

        physicsQuery.queryEntities();
        //this is called so that i can draw stuff
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

                i = 0;
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

                            math::vec3 worldStart = localTransform * math::vec4(*(edgeToExecuteOn->edgePositionPtr), 1);
                            math::vec3 worldEnd = localTransform * math::vec4(*(edgeToExecuteOn->nextEdge->edgePositionPtr), 1);

                            debug::drawLine(worldStart, worldEnd, usedColor, 2.0f, 0.0f, useDepth);

                        } while (initialEdge != currentEdge && currentEdge != nullptr);
                    }
                }

            }

        }

        //FindClosestPointsToLineSegment unit test


        math::vec3 p1(5, -0.5, 0);
        math::vec3 p2(5, 0.5, 0);

        math::vec3 p3(6, 0, -0.5);
        math::vec3 p4(6, 0, 0.5);

        math::vec3 p1p2;
        math::vec3 p3p4;

        debug::drawLine(p1, p2, math::colors::red, 5.0f);
        debug::drawLine(p3, p4, math::colors::red, 5.0f);

        physics::PhysicsStatics::FindClosestPointsToLineSegment(p1, p2, p3, p4, p1p2, p3p4);

        debug::drawLine(p1p2, p3p4, math::colors::green, 5.0f);

        p1 = math::vec3(8, 0, 0);
        p2 = p1 + math::vec3(0, 1.0f, 0);

        p3 = math::vec3(10, 0, 0) + math::vec3(1.0f);
        p4 = p3 - math::vec3(1.0f);

        debug::drawLine(p1, p2, math::colors::red, 5.0f);
        debug::drawLine(p3, p4, math::colors::red, 5.0f);

        physics::PhysicsStatics::FindClosestPointsToLineSegment(p1, p2, p3, p4, p1p2, p3p4);

        debug::drawLine(p1p2, p3p4, math::colors::green, 5.0f);

    }

    void onActivateUnitTest2(activate_CRtest2* action)
    {
        static bool activated = false;

        if (action->value)
        {
            activated = true;
            auto crb = m_ecs->createComponent<physics::rigidbody>(staticToOBBEnt);
        }
    }

    void onActivateUnitTest3(activate_CRtest3* action)
    {
        if (action->value)
        {
            auto crb = m_ecs->createComponent<physics::rigidbody>(staticToEdgeEnt);
        }
    }

    void onExtendedPhysicsContinueRequest(extendedPhysicsContinue* action)
    {
        if (action->value)
        {
            physics::PhysicsSystem::IsPaused = false;
        }
        else
        {
            physics::PhysicsSystem::IsPaused = true;
        }

    }

    void onNextPhysicsTimeStepRequest(nextPhysicsTimeStepContinue* action)
    {
        if (!(action->value))
        {
            physics::PhysicsSystem::IsPaused = true;
            physics::PhysicsSystem::oneTimeRunActive = true;
            log::debug(" onNextPhysicsTimeStepRequest");
        }

    }

    void FrictionTestActivate(activateFrictionTest* action)
    {
        if (action->value)
        {
            rotate = true;
            {
                auto rbH = m_ecs->createComponent<physics::rigidbody>(NoFrictionBody);

                auto rb = rbH.read();

                rb.friction = 0.0f;

                rbH.write(rb);
            }

            {
                auto rbH = m_ecs->createComponent<physics::rigidbody>(Point3FrictionBody);

                auto rb = rbH.read();

                rb.friction = 0.1f;

                rbH.write(rb);
            }

            {
                auto rbH = m_ecs->createComponent<physics::rigidbody>(Point6FrictionBody);

                auto rb = rbH.read();

                rb.friction = 0.4f;

                rbH.write(rb);
            }

            {
                auto rbH = m_ecs->createComponent<physics::rigidbody>(FullFrictionBody);

                auto rb = rbH.read();

                rb.friction = 1.0f;

                rbH.write(rb);
            }


        }
    }

    void CreateCubeStack(int height, int width, int breadth, math::vec3 startPosition, math::vec3 offset,
        physics::cube_collider_params cubeParams, float cubeFriction, rendering::model_handle cubeH, rendering::material_handle wireframeH, bool addRigidbody = true)
    {
        int l = 0;
        for (int i = 0; i < height; i++)
        {
            for (int j = 0; j < width; j++)
            {
                for (int k = 0; k < breadth; k++)
                {
                    auto ent = m_ecs->createEntity();

                    auto entPhyHande = ent.add_component<physics::physicsComponent>();

                    physics::physicsComponent physicsComponent;
                    physics::physicsComponent::init(physicsComponent);
                    physicsComponent.AddBox(cubeParams);
                    physicsComponent.isTrigger = false;
                    entPhyHande.write(physicsComponent);

                    auto [positionH, rotationH, scaleH] = m_ecs->createComponents<transform>(ent);
                    positionH.write(startPosition + math::vec3(j * offset.x, i * offset.y, k * offset.z));
                    scaleH.write(math::vec3(1.0f));

                    if (addRigidbody)
                    {
                        auto rbH = m_ecs->createComponent<physics::rigidbody>(ent);

                        auto rb = rbH.read();

                        rb.friction = cubeFriction;

                        rbH.write(rb);

                    }

                    auto idHandle = m_ecs->createComponent<physics::identifier>(ent);
                    auto idComp = idHandle.read();
                    idComp.id = "Stack" + std::to_string(l);
                    l++;
                    idHandle.write(idComp);

                }
            }
        }
    }

};
