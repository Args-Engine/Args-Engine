#pragma once

#include <core/events/event.hpp>
#include <core/math/math.hpp>
#include <application/input/inputmap.hpp>


namespace args::application
{
    template <class T>
    struct input_action : public events::event<T>
    {

        void set(bool v,inputmap::modifier_keys m,inputmap::method id)
        {
            value = v;
            mods = m;
            identifier = id;
        }

        bool value{};
        inputmap::modifier_keys mods{};
        inputmap::method identifier{};

        A_NODISCARD bool pressed()  const { return value == true;  }
        A_NODISCARD bool released() const { return value == false; }
    };

    template <class T>
    struct input_axis : public events::event<T>
    {
        void set(float v,inputmap::modifier_keys m,inputmap::method id)
        {
            value = v;
            mods = m;
            identifier = id;
        }

        float value{};
        inputmap::modifier_keys mods{};
        inputmap::method identifier{};
    };

}
