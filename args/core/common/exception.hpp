#pragma once
#include <exception>
#include <core/types/primitives.hpp>
#include <core/platform/platform.hpp>

#define args_exception args::core::exception(__FILE__, __LINE__, __FULL_FUNC__)
#define args_exception_msg(msg) args::core::exception(msg, __FILE__, __LINE__, __FULL_FUNC__)

#define args_invalid_fetch_error args::core::invalid_fetch_error(__FILE__, __LINE__, __FULL_FUNC__)
#define args_invalid_fetch_msg(msg) args::core::invalid_fetch_error(msg, __FILE__, __LINE__, __FULL_FUNC__)

#define args_invalid_component_error args::core::invalid_component_error(__FILE__, __LINE__, __FULL_FUNC__)
#define args_invalid_component_msg(msg) args::core::invalid_component_error(msg, __FILE__, __LINE__, __FULL_FUNC__)
#define args_component_destroyed_error args::core::component_destroyed_error(__FILE__, __LINE__, __FUNC__)
#define args_component_destroyed_msg(msg) args::core::component_destroyed_error(msg, __FILE__, __LINE__, __FULL_FUNC__)

#define args_invalid_entity_error args::core::invalid_entity_error(__FILE__, __LINE__, __FULL_FUNC__)
#define args_invalid_entity_msg(msg) args::core::invalid_entity_error(msg, __FILE__, __LINE__, __FULL_FUNC__)
#define args_entity_not_found_error args::core::entity_not_found_error(__FILE__, __LINE__, __FUNC__)
#define args_entity_not_found_msg(msg) args::core::entity_not_found_error(msg, __FILE__, __LINE__, __FULL_FUNC__)

#define args_invalid_creation_error args::core::invalid_creation_error(__FILE__, __LINE__, __FULL_FUNC__)
#define args_invalid_creation_msg(msg) args::core::invalid_creation_error(msg, __FILE__, __LINE__, __FULL_FUNC__)
#define args_entity_exists_error args::core::entity_exists_error(__FILE__, __LINE__, __FUNC__)
#define args_entity_exists_msg(msg) args::core::entity_exists_error(msg, __FILE__, __LINE__, __FULL_FUNC__)
#define args_component_exists_error args::core::component_exists_error(__FILE__, __LINE__, __FUNC__)
#define args_component_exists_msg(msg) args::core::component_exists_error(msg, __FILE__, __LINE__, __FULL_FUNC__)

#define args_invalid_type_error args::core::invalid_type_error(__FILE__, __LINE__, __FULL_FUNC__)
#define args_invalid_type_msg(msg) args::core::invalid_type_error(msg, __FILE__, __LINE__, __FULL_FUNC__)
#define args_unknown_component_error args::core::unknown_component_error(__FILE__, __LINE__, __FULL_FUNC__)
#define args_unknown_component_msg(msg) args::core::unknown_component_error(msg, __FILE__, __LINE__, __FULL_FUNC__)
#define args_unknown_system_error args::core::unknown_system_error(__FILE__, __LINE__, __FULL_FUNC__)
#define args_unknown_system_msg(msg) args::core::unknown_system_error(msg, __FILE__, __LINE__, __FULL_FUNC__)

namespace args::core
{
	class exception : public std::exception
	{
	private:
		cstring file;
		uint line;
		cstring func;
		cstring message;

	public:
		exception(cstring file, uint line, cstring func) : std::exception(), file(file), line(line), func(func), message("Args generic exception occurred.") {}
		exception(cstring msg, cstring file, uint line, cstring func) : std::exception(), file(file), line(line), func(func), message(msg) {}

		virtual cstring what() const noexcept override { return message; }
		cstring get_file() const { return file; }
		uint get_line() const { return line; }
		cstring get_func() const { return func; }
	};

#pragma region invalid fetch
	class invalid_fetch_error : public exception
	{
	public:
		invalid_fetch_error(cstring file, uint line, cstring func) : exception("Invalid fetch occurred.", file, line, func) {}
		invalid_fetch_error(cstring msg, cstring file, uint line, cstring func) : exception(msg, file, line, func) {}
	};

	class invalid_component_error : public invalid_fetch_error
	{
	public:
		invalid_component_error(cstring file, uint line, cstring func) : invalid_fetch_error("Component invalid.", file, line, func) {}
		invalid_component_error(cstring msg, cstring file, uint line, cstring func) : invalid_fetch_error(msg, file, line, func) {}
	};

	class component_destroyed_error : public invalid_component_error
	{
	public:
		component_destroyed_error(cstring file, uint line, cstring func) : invalid_component_error("Component no longer exists.", file, line, func) {}
		component_destroyed_error(cstring msg, cstring file, uint line, cstring func) : invalid_component_error(msg, file, line, func) {}
	};

	class invalid_entity_error : public invalid_fetch_error
	{
	public:
		invalid_entity_error(cstring file, uint line, cstring func) : invalid_fetch_error("Entity invalid.", file, line, func) {}
		invalid_entity_error(cstring msg, cstring file, uint line, cstring func) : invalid_fetch_error(msg, file, line, func) {}
	};

	class entity_not_found_error : public invalid_entity_error
	{
	public:
		entity_not_found_error(cstring file, uint line, cstring func) : invalid_entity_error("Entity does not exist.", file, line, func) {}
		entity_not_found_error(cstring msg, cstring file, uint line, cstring func) : invalid_entity_error(msg, file, line, func) {}
	};
#pragma endregion

#pragma region invalid creation
	class invalid_creation_error : public exception
	{
	public:
		invalid_creation_error(cstring file, uint line, cstring func) : exception("Creation invalid.", file, line, func) {}
		invalid_creation_error(cstring msg, cstring file, uint line, cstring func) : exception(msg, file, line, func) {}
	};

	class entity_exists_error : public invalid_creation_error
	{
	public:
		entity_exists_error(cstring file, uint line, cstring func) : invalid_creation_error("Entity already exist.", file, line, func) {}
		entity_exists_error(cstring msg, cstring file, uint line, cstring func) : invalid_creation_error(msg, file, line, func) {}
	};

	class component_exists_error : public invalid_creation_error
	{
	public:
		component_exists_error(cstring file, uint line, cstring func) : invalid_creation_error("Component already exist.", file, line, func) {}
		component_exists_error(cstring msg, cstring file, uint line, cstring func) : invalid_creation_error(msg, file, line, func) {}
	};
#pragma endregion

#pragma region invalid type
	class invalid_type_error : public exception
	{
	public:
		invalid_type_error(cstring file, uint line, cstring func) : exception("Type invalid.", file, line, func) {}
		invalid_type_error(cstring msg, cstring file, uint line, cstring func) : exception(msg, file, line, func) {}
	};

	class unknown_component_error : public invalid_type_error
	{
	public:
		unknown_component_error(cstring file, uint line, cstring func) : invalid_type_error("Unknown component type.", file, line, func) {}
		unknown_component_error(cstring msg, cstring file, uint line, cstring func) : invalid_type_error(msg, file, line, func) {}
	};

	class unknown_system_error : public invalid_type_error
	{
	public:
		unknown_system_error(cstring file, uint line, cstring func) : invalid_type_error("Unknown sytem type.", file, line, func) {}
		unknown_system_error(cstring msg, cstring file, uint line, cstring func) : invalid_type_error(msg, file, line, func) {}
	};
#pragma endregion

}
