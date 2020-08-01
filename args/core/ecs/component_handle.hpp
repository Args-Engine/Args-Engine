#pragma once
#include <atomic>
#include <core/common/exception.hpp>
#include <core/ecs/ecsregistry.hpp>
#include <core/ecs/component_container.hpp>
#include <core/ecs/entity_handle.hpp>
#include <core/platform/platform.hpp>

/**
 * @file component_handle.hpp
 */

namespace args::core::ecs
{
	/**@class component_handle_base
	 * @brief Base class of args::core::ecs::component_handle.
	 * @ref args::core::ecs::component_handle.
	 */
	class component_handle_base
	{
	public:
		// Entity that owns this component.
		const entity_handle entity_handle;

	protected:
		EcsRegistry& m_registry;
		id_type m_ownerId;

	public:
		component_handle_base(id_type entityId, EcsRegistry& registry) : entity_handle(registry.getEntity(entityId)), m_registry(registry), m_ownerId(entityId) {}

		/**@brief Checks if handle still points to a valid component.
		 */
		virtual bool valid() ARGS_IMPURE_RETURN(m_ownerId);

		/**@brief Checks if handle still points to a valid component.
		 */
		operator bool() { return valid(); }
	};

	/**@class component_handle
	 * @brief Handle to components that allow safe component loading and storing.
	 * @tparam component_type Type of targeted component.
	 */
	template<typename component_type>
	class component_handle : public component_handle_base
	{
	public:
		/**@brief Creates component handle for the given entity.
		 */
		component_handle(id_type entityId, EcsRegistry& registry) : component_handle_base(entityId, registry) {}

		/**@brief Atomic read of component.
		 * @param order Memory order at which to load the component.
		 * @returns component_type Current value of component.
		 */
		component_type read(std::memory_order order = std::memory_order_acquire)
		{
			async::transferable_atomic<component_type>* comp = m_registry.getFamily<component_type>()->get_component(m_ownerId);
			if (!comp)
				return component_type();

			async::readonly_guard guard(comp->get_lock());
			if (!valid())
				return component_type();

			return comp->get().load(order);
		}

		/**@brief Atomic write of component.
		 * @param value Value you wish to write.
		 * @param order Memory order at which to write the component.
		 * @returns component_type Current value of component.
		 */
		component_type write(component_type&& value, std::memory_order order = std::memory_order_release)
		{
			async::transferable_atomic<component_type>* comp = m_registry.getFamily<component_type>()->get_component(m_ownerId);
			if (!comp)
				return component_type();

			async::readonly_guard guard(comp->get_lock());
			if (!valid())
				return component_type();

			comp->get().store(value, order);

			return value;
		}

		/**@brief Atomic read modify write with add modification on component.
		 * @param value Value you wish to add.
		 * @param loadOrder Memory order at which to load the component.
		 * @param successOrder Memory order upon success of CAS-loop.
		 * @param failureOrder Memory order upon failure of CAS-loop.
		 * @returns component_type Current value of component.
		 */
		component_type fetch_add(component_type&& value,
			std::memory_order loadOrder = std::memory_order_acquire,
			std::memory_order successOrder = std::memory_order_release,
			std::memory_order failureOrder = std::memory_order_relaxed)
		{
			async::transferable_atomic<component_type>* comp = m_registry.getFamily<component_type>()->get_component(m_ownerId);
			if (!comp)
				return component_type();

			async::readonly_guard guard(comp->get_lock());
			if (!valid())
				return component_type();

			component_type oldVal = comp->get().load(loadOrder);
			component_type newVal = oldVal + value;

			// CAS loop to assure our modification will happen correctly without overwriting some other change.
			while (!comp->get().compare_exchange_strong(oldVal, newVal, successOrder, failureOrder))
				newVal = oldVal + value;

			return newVal;
		}

		/**@brief Atomic read modify write with multiply modification on component.
		 * @param value Value you wish to multiply by.
		 * @param loadOrder Memory order at which to load the component.
		 * @param successOrder Memory order upon success of CAS-loop.
		 * @param failureOrder Memory order upon failure of CAS-loop.
		 * @returns component_type Current value of component.
		 */
		component_type fetch_multiply(component_type&& value,
			std::memory_order loadOrder = std::memory_order_acquire,
			std::memory_order successOrder = std::memory_order_release,
			std::memory_order failureOrder = std::memory_order_relaxed)
		{
			async::transferable_atomic<component_type>* comp = m_registry.getFamily<component_type>()->get_component(m_ownerId);
			if (!comp)
				return component_type();

			async::readonly_guard guard(comp->get_lock());
			if (!valid())
				return component_type();

			component_type oldVal = comp->get().load(loadOrder);
			component_type newVal = oldVal * value;

			// CAS loop to assure our modification will happen correctly without overwriting some other change.
			while (!comp->get().compare_exchange_strong(oldVal, newVal, successOrder, failureOrder))
				newVal = oldVal * value;

			return newVal;
		}

		/**@brief Locks component family and destroys component.
		 * @ref args::core::ecs::component_container::destroy_component
		 */
		void destroy()
		{
			async::transferable_atomic<component_type>* comp = m_registry.getFamily<component_type>()->get_component(m_ownerId);
			if (!comp)
				return;

			async::readwrite_guard guard(comp->get_lock());
			if (valid())
				m_registry.destroyComponent<component_type>(m_ownerId);
		}

		/**@brief Checks if handle still points to a valid component.
		 */
		virtual bool valid() override
		{
			return m_ownerId && m_registry.getFamily<component_type>()->has_component(m_ownerId);
		}
	};
}