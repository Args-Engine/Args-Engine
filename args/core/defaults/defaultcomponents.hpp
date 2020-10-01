#pragma once
#include <core/math/math.hpp>
#include <core/ecs/archetype.hpp>

#include <core/logging/logging.hpp>

namespace args::core
{
    struct position : public math::vec3
    {
        position() : math::vec3(0, 0, 0) {}
        position(const position&) = default;
        position(position&&) = default;
        position(const math::vec3& src) : math::vec3(src) {}
        position(float x, float y, float z) : math::vec3(x, y, z) {}
        position(float v) : math::vec3(v) {}
        position& operator=(const position&) = default;
        position& operator=(position&&) = default;
        position& operator=(const math::vec3& src)
        {
            data = src.data;
            return *this;
        }
        position& operator=(math::vec3&& src)
        {
            data = src.data;
            return *this;
        }
    };

    struct rotation : public math::quat
    {
        rotation() : math::quat(1, 0, 0, 0) {}
        rotation(float w, float x, float y, float z) : math::quat(w, x, y, z) {}
        rotation(const rotation&) = default;
        rotation(rotation&&) = default;
        rotation(const math::quat& src) : math::quat(src) {}
        rotation& operator=(const rotation&) = default;
        rotation& operator=(rotation&&) = default;
        rotation& operator=(const math::quat& src)
        {
            data = src.data;
            return *this;
        }
        rotation& operator=(math::quat&& src)
        {
            data = src.data;
            return *this;
        }

        math::vec3 right()
        {
            return math::toMat3(*this) * math::vec3(1.f, 0.f, 0.f);
        }

        math::vec3 up()
        {
            return math::toMat3(*this) * math::vec3(0.f, 1.f, 0.f);
        }

        math::vec3 forward()
        {
            return math::toMat3(*this) * math::vec3(0.f, 0.f, 1.f);
        }

        math::mat3 matrix()
        {
            return math::toMat3(*this);
        }
    };

    struct scale : public math::vec3
    {
        scale() : math::vec3(1, 1, 1) {}
        scale(float x, float y, float z) : math::vec3(x, y, z) {}
        scale(float v) : math::vec3(v) {}
        scale(const scale&) = default;
        scale(scale&&) = default;
        scale(const math::vec3& src) : math::vec3(src) {}
        scale& operator=(const scale&) = default;
        scale& operator=(scale&&) = default;
        scale& operator=(const math::vec3& src)
        {
            data = src.data;
            return *this;
        }
        scale& operator=(math::vec3&& src)
        {
            data = src.data;
            return *this;
        }
    };

    struct transform : public ecs::archetype<position, rotation, scale>
    {
        using base = ecs::archetype<position, rotation, scale>;

        transform(const base::handleGroup& handles) : base(handles) {}

        std::tuple<position, rotation, scale> get_world_components()
        {
            math::mat4 worldMatrix = get_world_to_local_matrix();
            math::vec3 p;
            math::quat r;
            math::vec3 s;

            math::decompose(worldMatrix, s, r, p);

            return std::make_tuple<position, rotation, scale>(p, r, s);
        }

        math::mat4 get_world_to_local_matrix()
        {
            return math::inverse(get_local_to_world_matrix());
        }

        math::mat4 get_local_to_world_matrix()
        {
            auto& [positionH, rotationH, scaleH] = handles;
            auto parent = positionH.entity.get_parent();

            if (parent && parent.has_components<transform>())
            {
                transform transf = parent.get_component_handles<transform>();
                return transf.get_local_to_world_matrix() * math::compose(positionH.read(), rotationH.read(), scaleH.read());
            }

            return math::compose(positionH.read(), rotationH.read(), scaleH.read());
        }

        math::mat4 get_local_to_parent_matrix()
        {
            auto& [positionH, rotationH, scaleH] = handles;
            return math::compose(positionH.read(), rotationH.read(), scaleH.read());
        }
    };
}

namespace fmt
{
    using namespace args::core;

    template <>
    struct fmt::formatter<position> {
        char presentation = 'f';

        constexpr auto parse(format_parse_context& ctx) {

            auto it = ctx.begin(), end = ctx.end();
            if (it != end && (*it == 'f' || *it == 'e')) presentation = *it++;

            if (it != end && *it != '}')
                throw format_error("invalid format");

            return it;
        }

        template <typename FormatContext>
        auto format(const position& p, FormatContext& ctx) {
            return format_to(
                ctx.out(),

                presentation == 'f' ? "{:f}" : "{:e}",
                static_cast<math::vec3>(p));
        }
    };

    template <>
    struct fmt::formatter<rotation> {
        char presentation = 'f';

        constexpr auto parse(format_parse_context& ctx) {

            auto it = ctx.begin(), end = ctx.end();
            if (it != end && (*it == 'f' || *it == 'e')) presentation = *it++;

            if (it != end && *it != '}')
                throw format_error("invalid format");

            return it;
        }

        template <typename FormatContext>
        auto format(const rotation& r, FormatContext& ctx) {
            return format_to(
                ctx.out(),

                presentation == 'f' ? "{:f}" : "{:e}",
                static_cast<math::quat>(r));
        }
    };

    template <>
    struct fmt::formatter<scale> {
        char presentation = 'f';

        constexpr auto parse(format_parse_context& ctx) {

            auto it = ctx.begin(), end = ctx.end();
            if (it != end && (*it == 'f' || *it == 'e')) presentation = *it++;

            if (it != end && *it != '}')
                throw format_error("invalid format");

            return it;
        }

        template <typename FormatContext>
        auto format(const scale& s, FormatContext& ctx) {
            return format_to(
                ctx.out(),

                presentation == 'f' ? "{:f}" : "{:e}",
                static_cast<math::vec3>(s));
        }
    };
}
