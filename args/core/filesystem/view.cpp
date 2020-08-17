#include <algorithm>
#include <core/filesystem/view.hpp>
#include "navigator.hpp"
#include "provider_registry.hpp"
#include "detail/strpath_manip.hpp"


namespace args::core::filesystem
{
    view::operator bool() const
    {
        return is_valid();
    }

    bool view::is_valid(bool deep_check) const
    {
        //check if path is non empty & if 
        if (m_path.empty()) return false;
        if (!provider_registry::has_domain(get_domain())) return false;

        //if deep checking also check if the path meets the requirements of the
        //navigator system
        if(deep_check)
        {
            const navigator n(m_path);
            if (n.find_solution().has_err()) return false;
        }

        return true;
    }

    file_traits view::file_info()
    {
        //get solution
        auto result = make_solution();

        if (result.has_err()) return invalid_file_t;
        else
        {
            //get resolver
            auto resolver = build();
            if (resolver == nullptr) return invalid_file_t;

            //get traits
            return resolver->get_traits();
        }

    }

    filesystem_traits view::filesystem_info()
    {
        //get solution
        auto result = make_solution();

        if (result.has_err()) return invalid_filesystem_t;
        else
        {
            //get resolver
            const auto resolver = build();
            if (resolver == nullptr) return invalid_filesystem_t;

            //get traits
            return resolver->get_fs_traits();
        }
    }

    std::string view::get_domain() const
    {
        //string magic to find the first : & substr
        const auto idx = m_path.find_first_of(':');
        return m_path.substr(0, idx + 1) + strpath_manip::separator() + strpath_manip::separator();
    }

    common::result_decay_more<basic_resource, fs_error> view::get()
    {
        using common::Err, common::Ok;

        //decay overloads the operator of ok_type and operator== for valid_t
        using decay = common::result_decay_more<basic_resource, fs_error>;

        //get solution
        auto result = make_solution();
        if (result.has_err()) Err(result.get_error());

        //get resolver of solution
        auto resolver = build();
        if (resolver == nullptr) return decay(Err(args_fs_error("unable to get required filesystem to get resource!")));

        //get & check traits
        const auto traits = resolver->get_traits();
        if (traits.is_valid && traits.exists && traits.can_be_read)
        {
            //wrap get in decay
            return decay(resolver->get());
        }
        return decay(Err(args_fs_error("invalid file traits: (not valid) or (does not exist) or (cannot be read)")));
    }

    common::result<void, fs_error> view::set(const basic_resource& resource)
    {
        using common::Ok, common::Err;

        //get solution
        auto result = make_solution();
        if (result.has_err()) return Err_of(result);

        //get resolver of solution
        auto resolver = build();
        if (resolver == nullptr) return Err(args_fs_error("unable to get required filesystem to set resource!"));

        //get & check traits
        const auto traits = resolver->get_traits();
        if (traits.is_valid && ((traits.can_be_written && !traits.is_directory) || traits.can_be_created))
        {
            //set
            return resolver->set(resource);
        }
        return Err(args_fs_error("invalid file traits: (not valid) or (not writeable or directory) or (not creatable)"));
    }

    view view::parent() const
    {
        //get parent path
        const auto p_path = strpath_manip::parent(m_path);
        return view(p_path);
    }

    view view::find(std::string_view identifier) const
    {
        //probably not necessarily necessary
        std::string sanitized = strpath_manip::sanitize(std::string(identifier));

        //basic bails
        if (sanitized == "..") return parent();
        if (sanitized == ".")  return *this;

        //subdirectorize and then sanitize
        sanitized = strpath_manip::subdir(m_path, sanitized);
        sanitized = strpath_manip::sanitize(sanitized, true);

        return view(sanitized);
    }

    view view::operator[](std::string_view identifier) const
    {
        return find(identifier);
    }

    std::string view::create_identifier(const navigator::solution::iterator& e)
    {
        //iterate through path and create the ident for the provider
        std::string result;
        for (auto iter = m_foundSolution.begin(); iter != e; ++iter)
        {
            result += iter->second;
        }
        return result;
    }

    //TODO(algo-ryth-mix): the navigator should probably return a more efficient
    //TODO(cont.)          representation to begin with
    std::shared_ptr<filesystem_resolver> view::build()
    {
        //first check if a solution even exists
        if (m_foundSolution.size() == 0)
        {
            return nullptr;
        }


        //second check if we even need to to resolution
        if (m_foundSolution.size() == 1)
        {
            auto& [r, path] = m_foundSolution.front();

            auto resolver = std::shared_ptr<filesystem_resolver>(r->make());
            resolver->set_target(path);

            return resolver;
        }

        //translate the solution into a resolution chain
        auto chain = translate_solution();

        if (!chain) return nullptr;

        //traverse resolution chain
        for (; chain->next != nullptr; chain = chain->next)
        {

            auto data = chain->provider->get();
            if (data.has_err()) return nullptr;

            //convert result -> resource -> data
            //and set as disk dat for subject
            chain->subject->set_disk_data(data.get().get());

        }

        //do it one last time for the last subject

        auto data = chain->provider->get();
        if (data.has_err()) return nullptr;

        //convert result -> resource -> data
        //and set as disk dat for subject
        chain->subject->set_disk_data(data.get().get());

        return chain->subject;

    }


    void view::make_inheritance()
    {

        //make all higher level fs inherit the traits from the lower level
        for (std::size_t i = 0; i < m_foundSolution.size() - 1; ++i)
        {
            m_foundSolution.at(i + 1).first->inherit(*m_foundSolution.at(i).first);
        }
    }

    std::shared_ptr<view::create_chain> view::translate_solution()
    {
        //this is a more approachable representation
         //of the solution
        std::shared_ptr<create_chain> chain = nullptr;

        make_inheritance();

        for (auto iter = m_foundSolution.rbegin(); iter != m_foundSolution.rend(); ++iter)
        {
            std::string identifier = create_identifier((iter + 1).base());
            auto& [resolver, resolver_path] = *iter;



            //we expect the first element to be valid no matter what if it isn't we have a deeper problem
            if (iter != m_foundSolution.rend() - 1)
            {
                auto* memory_resolver = dynamic_cast<mem_filesystem_resolver*>(resolver);

                //we were unable to create a complete chain as one of the resolver was not a memory
                //resolver when it should have been.
                if (memory_resolver == nullptr) return nullptr;
                auto shared_memory_resolver = std::shared_ptr<mem_filesystem_resolver>(memory_resolver->make_higher());
                shared_memory_resolver->set_target(resolver_path);
                shared_memory_resolver->set_identifier(identifier);
                std::shared_ptr<create_chain> previous = nullptr;

                //finish setting up previous link
                if (chain)
                {
                    chain->provider = shared_memory_resolver;
                    previous = chain;
                }

                //create new link
                chain = std::make_shared<create_chain>();
                chain->subject = shared_memory_resolver;
                chain->next = previous;


                //prewarming also does a in cache check, which we can use to see if
                //this resource is already loaded
                if (memory_resolver->prewarm())
                {

                    //if the resource was in cache we can use it to generate the higher levels
                    //and don't need to resolve all the way to the root

                    //we also need to skip this node, since it will be the root
                    chain = chain->next;
                    break;
                }

            }
            else
            {
                //all caches missed
                chain->provider = std::shared_ptr<filesystem_resolver>(resolver->make());
                chain->provider->set_identifier(resolver->get_identifier());
                chain->provider->set_target(resolver_path);
            }
        }
        return chain;
    }

    common::result<void, fs_error> view::make_solution()
    {
        using common::Ok;

        //check if a solution already exists
        if (m_foundSolution.empty())
        {

            //create solution using navigator
            const navigator n(m_path);
            auto solution = n.find_solution();

            if (solution.has_err())
                return Err_of(solution);

            m_foundSolution = solution.get();
        }
        //return empty ok
        return Ok();
    }
}
