#pragma once
#include <core/containers/delegate.hpp>
#include <core/containers/sparse_map.hpp>
#include <core/containers/hashed_sparse_set.hpp>
#include <core/types/types.hpp>
#include <core/events/event.hpp>

#include <memory>

/**@file eventbus.hpp
 */

namespace args::core::events
{
	/**@class EventBus
	 * @brief Central communication channel for events and messages.
	 */
	class EventBus
	{
		sparse_map<id_type, hashed_sparse_set<event_base*>> m_events;
		sparse_map<id_type, multicast_delegate<void(EventBus*)>> m_eventCallbacks;

	public:
		~EventBus()
		{
			for (auto& events : m_events)
			{
				for (auto* event : events)
					delete event;

				events.clear();
			}
		}

		/**@brief Insert event into bus and notify all subscribers.
		 * @tparam event_type Event type to raise.
		 * @param arguments Arguments to pass to the constructor of the event.
		 */
		template<typename event_type, typename... Args, typename = inherits_from<event_type, event<event_type>>>
		void raiseEvent(Args&&... arguments)
		{
			event_type* event = new event_type(arguments...); // Create new event.

			if (event->unique() && m_events[event_type::id].size()) // Check if it's unique and if another one already exists.
			{
				m_eventCallbacks[event_type::id].invoke(this); // Notify and delete.
				delete event;
				return;
			}

			m_events[event_type::id].insert(event); // Keep the event stored. (Or at least keep it somewhere fetch able.)

			m_eventCallbacks[event_type::id].invoke(this); // Notify.

			if (!event->persistent()) // If the event isn't persistent then delete it again.
			{
				m_events[event_type::id].erase(event);
				delete event;
			}
		}

		/**@brief Check if an event is active.
		 * @tparam event_type Event type to check for.
		 */
		template<typename event_type, inherits_from<event_type, event<event_type>> = 0>
		bool checkEvent() const
		{
			return m_events.contains(event_type::id) && m_events[event_type::id].size();
		}

		/**@brief Get the amount of events/messages that are currently in the bus.
		 * @tparam event_type Event type to get the amount of.
		 */
		template<typename event_type, inherits_from<event_type, event<event_type>> = 0>
		size_type getEventCount() const
		{
			if (m_events.contains(event_type::id))
				return m_events[event_type::id].size();
			return 0;
		}

		/**@brief Get a reference to an event/message of a certain index.
		 * @tparam event_type Event type to fetch.
		 */
		template<typename event_type, inherits_from<event_type, event<event_type>> = 0>
		const event_type& getEvent(index_type index = 0) const
		{
			if (checkEvent<event_type>())
				return *static_cast<event_type*>(m_events[event_type::id][index]); // Static cast because we already know that the types are the same.
			return nullptr;
		}

		/**@brief Get a reference to the most recently raised event of this type.
		 * @tparam event_type Event type to fetch.
		 */
		template<typename event_type, inherits_from<event_type, event<event_type>> = 0>
		const event_type& getLastEvent() const
		{
			if (checkEvent<event_type>())
			{
				size_type size = m_events[event_type::id].size();
				return *static_cast<event_type*>(m_events[event_type::id][size - 1]); // Static cast because we already know that the types are the same.
			}
			return nullptr;
		}

		/**@brief Removes a certain event from the bus.
		 * @tparam event_type Event type to clear.
		 */
		template<typename event_type, inherits_from<event_type, event<event_type>> = 0>
		void clearEvent(index_type index = 0)
		{
			if (checkEvent<event_type>())
			{
				auto* event = m_events[event_type::id].dense()[index];
				m_events[event_type::id].erase(event);
				delete event;
			}
		}

		/**@brief Removes a the most recent event from the bus.
		 * @tparam event_type Event type to clear.
		 */
		template<typename event_type, inherits_from<event_type, event<event_type>> = 0>
		void clearLastEvent()
		{
			if (checkEvent<event_type>())
			{
				auto* event = m_events[event_type::id].dense()[m_events[event_type::id].size() - 1];
				m_events[event_type::id].erase(event);
				delete event;
			}
		}

		/**@brief Link a callback to an event type in order to get notified whenever one gets raised.
		 * @tparam event_type Event type to subscribe to.
		 */
		template<typename event_type, inherits_from<event_type, event<event_type>> = 0>
		void bindToEvent(delegate<void(EventBus*)> callback)
		{
			m_eventCallbacks[event_type::id] += callback;
		}
	};
}
