#include <core/ecs/entity_handle.hpp>
#include <core/ecs/ecsregistry.hpp>
#include <core/ecs/component_handle.hpp>

namespace args::core::ecs
{
	inline entity_handle& entity_handle::operator=(const entity_handle& other)
	{
		m_id = other.m_id;
		m_registry = other.m_registry;
		return *this;
	}

	A_NODISCARD inline const hashed_sparse_set<id_type>& entity_handle::component_composition() const
	{
		if (!m_registry)
			throw args_invalid_entity_error;
		return m_registry->getEntityData(m_id).components;
	}

	A_NODISCARD id_type entity_handle::get_id() const
	{
		if (valid())
			return m_id;
		return invalid_id;
	}

	A_NODISCARD sparse_map<id_type, entity_handle>::const_iterator entity_handle::begin() const
	{
		if (!m_registry)
			throw args_invalid_entity_error;
		return m_registry->getEntityData(m_id).children.begin();
	}

	A_NODISCARD sparse_map<id_type, entity_handle>::const_iterator entity_handle::end() const
	{
		if (!m_registry)
			throw args_invalid_entity_error;
		return m_registry->getEntityData(m_id).children.end();
	}

	A_NODISCARD entity_handle entity_handle::get_parent() const
	{
		if (!m_registry)
			throw args_invalid_entity_error;
		return m_registry->getEntity(m_registry->getEntityData(m_id).parent);
	}

	inline void entity_handle::set_parent(id_type newParent) const
	{
		if (!m_registry)
			throw args_invalid_entity_error;
		entity_data& data = m_registry->getEntityData(m_id);

		if (m_registry->validateEntity(data.parent))
			m_registry->getEntityData(data.parent).children.erase(m_id);

		if (m_registry->validateEntity(newParent))
		{
			data.parent = newParent;

			m_registry->getEntityData(data.parent).children.insert(m_id, *this);
		}
		else
			data.parent = invalid_id;
	}

	A_NODISCARD inline entity_handle entity_handle::operator[](index_type index) const
	{
		return get_child(index);
	}

	A_NODISCARD inline entity_handle entity_handle::get_child(index_type index) const
	{
		if (!m_registry)
			throw args_invalid_entity_error;

		sparse_map<id_type, entity_handle>& children = m_registry->getEntityData(m_id).children;
		if (index >= children.size())
			throw std::out_of_range("Child index out of range.");

		return children.dense()[index];
	}

	A_NODISCARD size_type entity_handle::child_count() const
	{
		if (!m_registry)
			throw args_invalid_entity_error;
		return m_registry->getEntityData(m_id).children.size();
	}

	inline void entity_handle::add_child(id_type childId) const
	{
		if (!m_registry)
			throw args_invalid_entity_error;
		entity_data& data = m_registry->getEntityData(m_id);
		if (!data.children.contains(childId))
			data.children[childId].set_parent(m_id);
	}

	inline void entity_handle::remove_child(id_type childId) const
	{
		if (!m_registry)
			throw args_invalid_entity_error;
		entity_data& data = m_registry->getEntityData(m_id);
		if (!data.children.contains(childId))
			data.children[childId].set_parent(invalid_id);
	}

	A_NODISCARD inline bool entity_handle::has_component(id_type componentTypeId) const
	{
		if (!m_registry)
			throw args_invalid_entity_error;
		return m_registry->getEntityData(m_id).components.contains(componentTypeId);
	}

	A_NODISCARD inline component_handle_base entity_handle::get_component(id_type componentTypeId) const
	{
		if (!m_registry)
			throw args_invalid_entity_error;
		return m_registry->getComponent(m_id, componentTypeId);
	}

	inline component_handle_base entity_handle::add_component(id_type componentTypeId) const
	{
		if (!m_registry)
			throw args_invalid_entity_error;
		return m_registry->createComponent(m_id, componentTypeId);
	}

	void entity_handle::remove_component(id_type componentTypeId) const
	{
		if (!m_registry)
			throw args_invalid_entity_error;
		m_registry->destroyComponent(m_id, componentTypeId);
	}

	inline void entity_handle::destroy(bool recurse) const
	{
		if (!m_registry)
			throw args_invalid_entity_error;
		m_registry->destroyEntity(m_id);
	}

	inline bool entity_handle::valid() const
	{
		if (m_registry && m_id)
			if (m_registry->validateEntity(m_id))
				return true;

		return false;
	}
}