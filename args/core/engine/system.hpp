#pragma once
#include <core/platform/platform.hpp>
#include <core/types/primitives.hpp>
#include <core/types/type_util.hpp>
#include <core/ecs/ecsregistry.hpp>
#include <core/scheduling/scheduler.hpp>
#include <core/events/eventbus.hpp>
#include <memory>

namespace args::core
{
	class SystemBase
	{
		friend class Module;
	protected:
		ecs::EcsRegistry* m_ecs;
		scheduling::Scheduler* m_scheduler;
		events::EventBus* m_eventBus;

		sparse_map<id_type, std::unique_ptr<scheduling::Process>> m_processes;

	public:
		const id_type id;
		const std::string name;

		SystemBase(id_type id, const std::string& name) : id(id), name(name) {}

		virtual void setup() ARGS_PURE;

		virtual ~SystemBase() = default;
	};

	template<typename SelfType>
	class System : public SystemBase
	{
	protected:
		template <void(SelfType::* func_type)(time::time_span<fast_time>), size_type charc>
		void createProcess(const char(&chainName)[charc], time::time_span<fast_time> interval = 0)
		{
			std::string name = std::string(chainName) + undecoratedTypeName<SelfType>() + std::to_string(interval);
			id_type id = nameHash(name);

			std::unique_ptr<scheduling::Process> process = std::make_unique<scheduling::Process>(name, id, interval);
			process->setOperation(delegate<void(time::time_span<fast_time>)>::create<SelfType, func_type>((SelfType*)this));
			m_processes.insert(id, std::move(process));

			m_scheduler->hookProcess<charc>(chainName, m_processes[id].get());
		}

		void createProcess(cstring chainName, delegate<void(time::time_span<fast_time>)>&& operation, time::time_span<fast_time> interval = 0)
		{
			std::string name = std::string(chainName) + undecoratedTypeName<SelfType>() + std::to_string(interval);
			id_type id = nameHash(name);

			std::unique_ptr<scheduling::Process> process = std::make_unique<scheduling::Process>(name, id, interval);
			process->setOperation(operation);
			m_processes.insert(id, std::move(process));

			m_scheduler->hookProcess(chainName, m_processes[id].get());
		}

		template<typename... component_types>
		A_NODISCARD ecs::EntityQuery createQuery()
		{
			return m_ecs->createQuery<component_types...>();
		}

		template<typename event_type, typename... Args, inherits_from<event_type, events::event<event_type>> = 0>
		void raiseEvent(Args... arguments)
		{
			m_eventBus->raiseEvent<event_type>(arguments...);
		}

		template<typename event_type, inherits_from<event_type, events::event<event_type>> = 0>
		A_NODISCARD bool checkEvent() const
		{
			return m_eventBus->checkEvent<event_type>();
		}

		template<typename event_type, inherits_from<event_type, events::event<event_type>> = 0>
		A_NODISCARD size_type getEventCount() const
		{
			return m_eventBus->getEventCount<event_type>();
		}

		template<typename event_type, inherits_from<event_type, events::event<event_type>> = 0>
		A_NODISCARD const event_type& getEvent(index_type index = 0) const
		{
			return m_eventBus->getEvent<event_type>(index);
		}

		template<typename event_type, inherits_from<event_type, events::event<event_type>> = 0>
		A_NODISCARD const event_type& getLastEvent() const
		{
			return m_eventBus->getLastEvent<event_type>();
		}

		template<typename event_type, inherits_from<event_type, events::event<event_type>> = 0>
		void clearEvent(index_type index = 0)
		{
			m_eventBus->clearEvent<event_type>(index);
		}

		template<typename event_type, inherits_from<event_type, events::event<event_type>> = 0>
		void clearLastEvent()
		{
			m_eventBus->clearLastEvent<event_type>();
		}

		template <void(SelfType::* func_type)(events::EventBus*), typename event_type, inherits_from<event_type, events::event<event_type>> = 0>
		void createProcess()
		{
			m_eventBus->bindToEvent<event_type>(delegate<void(events::EventBus*)>::create<SelfType, func_type>((SelfType*)this));
		}

		template<typename event_type, inherits_from<event_type, events::event<event_type>> = 0>
		void bindToEvent(delegate<void(events::EventBus*)> callback)
		{
			m_eventBus->bindToEvent<event_type>(callback);
		}

	public:
		System() : SystemBase(typeHash<SelfType>(), undecoratedTypeName<SelfType>()) {}
		virtual ~System() = default;
	};
}