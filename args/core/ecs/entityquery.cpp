#include <core/ecs/entityquery.hpp>
#include <core/ecs/queryregistry.hpp>
#include <core/ecs/ecsregistry.hpp>
#include <core/ecs/entity_handle.hpp>
#include <core/containers/hashed_sparse_set.hpp>

namespace args::core::ecs
{
    EntityQuery::EntityQuery(id_type id, QueryRegistry* registry, EcsRegistry* ecsRegistry) : m_registry(registry), m_ecsRegistry(ecsRegistry), m_id(id)
    {
        m_registry->addReference(m_id);
    }

    EntityQuery::EntityQuery(EntityQuery&& other)
    {
        m_id = other.m_id; 
        m_registry = other.m_registry;
        m_ecsRegistry = other.m_ecsRegistry;
        other.m_id = invalid_id;
    }

    EntityQuery::EntityQuery(const EntityQuery& other)
    {
        m_id = other.m_id; 
        m_registry = other.m_registry;
        m_ecsRegistry = other.m_ecsRegistry; 
        m_registry->addReference(m_id);
    }

    EntityQuery EntityQuery::operator=(EntityQuery&& other)
    {
        m_id = other.m_id;
        m_registry = other.m_registry;
        m_ecsRegistry = other.m_ecsRegistry;
        other.m_id = invalid_id;
        return *this;
    }

    EntityQuery EntityQuery::operator=(const EntityQuery& other)
    {
        m_id = other.m_id;
        m_registry = other.m_registry;
        m_ecsRegistry = other.m_ecsRegistry;
        m_registry->addReference(m_id);
        return *this;
    }

    EntityQuery::~EntityQuery()
    {
        if (QueryRegistry::isValid(m_registry))
        {
            m_registry->removeReference(m_id);
        }
    }

    sparse_map<id_type, entity_handle>::const_iterator EntityQuery::begin() const
    {
        return m_registry->getEntities(m_id).begin();
    }

    sparse_map<id_type, entity_handle>::const_iterator EntityQuery::end() const
    {
        return m_registry->getEntities(m_id).end();
    }

    inline void EntityQuery::addComponentType(id_type componentTypeId)
    {
        hashed_sparse_set<id_type> componentTypes;
        if (m_id)
            componentTypes = m_registry->getComponentTypes(m_id); // If this query is a valid query fetch a copy of the component types we're already querying.

        if (componentTypes.contains(componentTypeId)) // If we're already querying this component type then we don't need to do anything.
            return;

        componentTypes.insert(componentTypeId); // Add new type to query for.

        id_type newId = m_registry->getQueryId(componentTypes); // Check registry for an existing query with the desired component combination.
        if (newId) // If an existing query was found move this handle to that query instead.
        {
            m_registry->removeReference(m_id);
            m_id = newId;
            m_registry->addReference(m_id);
        }
        else if (m_id != invalid_id && m_registry->getReferenceCount(m_id) == 1) // If the query is valid and there are no other handles then we can safely add another component type to this query.
        {
            m_registry->addComponentType(m_id, componentTypeId);
        }
        else // If we cannot add a component type, and there is no existing query for us, then we should just create a new query with the desired component types.
        {
            m_registry->removeReference(m_id);
            m_id = m_registry->addQuery(componentTypes);
            m_registry->addReference(m_id);
        }

    }

    inline void EntityQuery::removeComponentType(id_type componentTypeId)
    {
        if (!m_id) // We're not pointing to a valid query, there's nothing to remove from.
            return;

        hashed_sparse_set<id_type> componentMap = m_registry->getComponentTypes(m_id); // Fetch the list of components we're querying.

        if (!componentMap.contains(componentTypeId)) // If we're not querying the component we wish to stop querying then we don't need to do anything.
            return;

        componentMap.erase(componentTypeId); // Remove component type from the list.

        id_type newId = m_registry->getQueryId(componentMap); // Check registry for an existing query with the desired component combination.
        if (newId) // If an existing query was found move this handle to that query instead.
        {
            m_registry->removeReference(m_id);
            m_id = newId;
            m_registry->addReference(m_id);
        }
        else if (m_id != invalid_id && m_registry->getReferenceCount(m_id) == 1) // If the query is valid and there are no other handles then we can safely remove the component type from this query.
        {
            m_registry->removeComponentType(m_id, componentTypeId);
        }
        else // If we cannot remove a component type, and there is no existing query for us, then we should just create a new query with the desired component types.
        {
            m_registry->removeReference(m_id);
            m_id = m_registry->addQuery(componentMap);
            m_registry->addReference(m_id);
        }

    }

    inline entity_handle EntityQuery::operator[](size_type index)
    {
        return m_ecsRegistry->getEntity(m_registry->getEntities(m_id).dense()[index]);
    }

    inline size_type EntityQuery::size()
    {
        return m_registry->getEntities(m_id).size();
    }
}
