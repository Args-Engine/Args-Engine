#pragma once
#include <core/filesystem/filesystem.hpp>

#include <iostream>

#include "core/filesystem/basic_resolver.hpp"

inline namespace {

using namespace ::args::core;
namespace fs = ::args::core::filesystem;

class mock_resolver final : public fs::mem_filesystem_resolver
{
public:
    mock_resolver(std::shared_ptr<byte_vec> v) : mem_filesystem_resolver(v){}
    mock_resolver(const mock_resolver& other) = default;
    mock_resolver(mock_resolver&& other) noexcept = default;
    mock_resolver& operator=(const mock_resolver& other) = default;
    mock_resolver& operator=(mock_resolver&& other) noexcept = default;
    ~mock_resolver() = default;


    [[nodiscard]] bool is_file() const noexcept override
    {
        if(get_target() == "test.txt") return true;
        return false;
    }
    [[nodiscard]] bool is_directory() const noexcept override { return false; }
    [[nodiscard]] bool is_valid() const noexcept override
    {
        if(get_target() == "test.txt") return true;
        return false;
    }
    [[nodiscard]] bool writeable() const noexcept override { return false; }
    [[nodiscard]] bool readable() const noexcept override {
        if(get_target() == "test.txt") return true;
        return false;
    }
    [[nodiscard]] bool creatable() const noexcept override { return false; }
    [[nodiscard]] bool exists() const noexcept override {
        if(get_target() == "test.txt") return true;
        return false;
    }
    [[nodiscard]] std::set<std::string> ls() const noexcept override { return {"test.txt"}; }

    common::result<filesystem::basic_resource, fs_error>
    get(interfaces::implement_signal_t) const noexcept override
    {
        
        using common::Err, common::Ok;

        if(!prewarm()) return Err(args_fs_error("was unable to cook the data!"));

        const auto result = get_data();

        if(get_target() == "test.txt") return Ok(filesystem::basic_resource(result));
        return Err(args_fs_error("this mock interface does not support file access"));
    }

    common::result<void,fs_error> set(interfaces::implement_signal_t, const filesystem::basic_resource& res) override
    {
        (void) res;
        return common::Err(args_fs_error("not implemented for the mock interface"));
    }
    void erase(interfaces::implement_signal_t) const noexcept override
    {
    }


    [[nodiscard]] mem_filesystem_resolver* make_higher() override
    {
        return new mock_resolver(nullptr);
    }
protected:
    void build_memory_representation(std::shared_ptr<const byte_vec> in, std::shared_ptr<byte_vec> out) const override
    {
        out->assign(in->begin(),in->end());
    }
};

}

std::ostream& operator<<(std::ostream& lhs, filesystem::basic_resource rhs)
{
    lhs << rhs.to_string();
    return lhs;
}

inline void test_filesystem()
{
    std::string  x = "always has been!";

    auto vector = std::make_shared<byte_vec>(x.begin(),x.end());


    fs::provider_registry::domain_create_resolver<fs::basic_resolver>("basic://","./assets");
    fs::provider_registry::domain_create_resolver<mock_resolver>("only_test_valid://",vector);
    fs::provider_registry::domain_create_resolver<mock_resolver>(".args-test",nullptr);



    auto nested_contents = fs::view("basic://config/test.args-test/test.txt").get();
    if(nested_contents == common::valid)
    {
        std::cout << "[arm-fs] -> found file !" << std::endl;
        std::cout << "[arm-fs] -> " << nested_contents << std::endl;
    }
    else
    {
        std::cout << "[arm-fs] -> " << nested_contents.get_error().what();
        __debugbreak();
    }

    auto traits = fs::view("only_test_valid://test.txt").file_info();
    if(traits != fs::invalid_file_t)
    {
        std::cout << "[arm-fs] -> " << fs::view("only_test_valid://test.txt").get() << std::endl;
    }
    else __debugbreak();

    auto contents = fs::view("basic://config/test.txt").get();

    if(contents == common::valid)
    {
        std::cout << "[arm-fs] -> found file !" << std::endl;
        std::cout <<"[arm-fs] -> " << contents << std::endl;
    }
    else
    {
        std::cout << "[arm-fs] -> " << contents.get_error().what();
        __debugbreak();
    }

    auto result = fs::view("basic://config/test2.txt").set(contents);
    if(result.has_err()) __debugbreak();

}

