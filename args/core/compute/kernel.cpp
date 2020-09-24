#include <core/compute/kernel.hpp>
#include <core/compute/program.hpp>
#include <core/logging/logging.hpp>

namespace args::core::compute
{
    Kernel& Kernel::build_buffer_names()
    {
        size_type size;
        cl_uint num_args;
        std::string container;

        //get the number of kernel arguments
        clGetKernelInfo(m_func, CL_KERNEL_NUM_ARGS, sizeof(cl_uint), &num_args, nullptr);

        for (cl_uint i = 0; i < num_args; ++i) {

            //get the length of the kernel argument
            if (clGetKernelArgInfo(m_func, i, CL_KERNEL_ARG_NAME, 0, nullptr, &size) == CL_KERNEL_ARG_INFO_NOT_AVAILABLE) {
                log::error("cannot get Argument name! Was the kernel built with -cl-kernel-arg-info ? ");
                continue;
            }

            //get the kernel argument
            container.resize(size, '\0');
            clGetKernelArgInfo(m_func, i, CL_KERNEL_ARG_NAME, size, reinterpret_cast<void*>(container.data()), nullptr);

            //remove trailing '\0'
            container.resize(container.size() - 1);

            //store it in the map
            m_paramsMap[container] = i;

        }
        return *this;
    }

    Kernel& Kernel::read_write_mode(buffer_type t)
    {
        m_default_mode = t;
        return *this;
    }

    Kernel& Kernel::set_buffer(Buffer buffer)
    {
        //get name from buffer
        if (buffer.m_name.empty())
            log::warn("Encountered unnamed buffer! binding to a Kernel-location will fail!");
        return set_buffer(buffer, buffer.m_name);
    }

    Kernel& Kernel::set_buffer(Buffer buffer, const std::string& name)
    {
        //translate name to index
        param_find([this, b = std::forward<Buffer>(buffer)](cl_uint index)
        {
            this->set_buffer(b, index);
        }, name);
        return *this;
    }

    Kernel& Kernel::set_buffer(Buffer buffer, cl_uint index)
    {
        //set kernel argument
        const cl_int ret = clSetKernelArg(m_func, index, sizeof(cl_mem), &buffer.m_memory_object);

        //check clSetKernelArg
        if (ret != CL_SUCCESS)
        {
            log::error("clSetKernelArg {}", ret);
        }
        return *this;
    }

    Kernel& Kernel::enqueue_buffer(Buffer buffer, block_mode blocking)
    {
        //get name from buffer
        if (buffer.m_name.empty())
            log::warn("Encountered unnamed buffer! binding to a Kernel-location will fail!");
        return enqueue_buffer(buffer, buffer.m_name, blocking);
    }

    Kernel& Kernel::enqueue_buffer(Buffer buffer, const std::string& name, block_mode blocking)
    {
        //translate name to index
        param_find([this, b = std::forward<Buffer>(buffer), blocking](cl_uint index)
        {
            this->enqueue_buffer(b, index, blocking);
        }, name);
        return *this;
    }

    Kernel& Kernel::enqueue_buffer(Buffer buffer, cl_uint index, block_mode blocking)
    {

        /*
            we are writing to the device buffer

            even though the buffer was r/w on
            the host it has now become read-only
            on the device

            //////////               //////////
            |  HOST  | -> write  ->  | DEVICE |
            | Buffer |               | Buffer |
            //////////               //////////

            the same happens in the other direction

            //////////               //////////
            |  HOST  |  <- read <-   | DEVICE |
            | Buffer |               | Buffer |
            //////////               //////////

            Here we are reading from the device
            so the buffer on the gpu was write-only
            (aka it can only write to the host buffer)
            but we still need to read it

        */
        cl_int ret;
        switch (buffer.m_type)
        {

            //buffer was read only
        case CL_MEM_READ_ONLY:
            //we "write" to a read only buffer because it is readonly for the kernel
            ret = clEnqueueWriteBuffer(m_queue, buffer.m_memory_object, static_cast<cl_bool>(blocking), 0, buffer.m_size, buffer.m_data, 0, nullptr, nullptr);
            break;

            //buffer was write only
        case CL_MEM_WRITE_ONLY:
            //similarly we read from a buffer that was write-only for the kernel 
            ret = clEnqueueReadBuffer(m_queue, buffer.m_memory_object, static_cast<cl_bool>(blocking), 0, buffer.m_size, buffer.m_data, 0, nullptr, nullptr);
            break;

            //buffer is read and write the default read/write mode needs to decide
        case CL_MEM_READ_WRITE:
            switch (m_default_mode)
            {
            case buffer_type::READ_BUFFER:
                //the mode was read on the host so we write to the kernel
                ret = clEnqueueWriteBuffer(m_queue, buffer.m_memory_object, static_cast<cl_bool>(blocking), 0, buffer.m_size, buffer.m_data, 0, nullptr, nullptr);
                break;

            case buffer_type::WRITE_BUFFER:
                //again the mode was write on the host so we read from it
                ret = clEnqueueReadBuffer(m_queue, buffer.m_memory_object, static_cast<cl_bool>(blocking), 0, buffer.m_size, buffer.m_data, 0, nullptr, nullptr);
                break;

            default:
                throw std::logic_error("Buffer enqueue type was neither read nor write ... I don't know what to do!");
            }
            break;
        default: throw std::logic_error("Buffer was neither read nor write nor readwrite");
        }

        if (ret != CL_SUCCESS)
            log::error("clEnqueueXXXXBuffer {}", ret);

        return *this;
    }

    Kernel& Kernel::set_and_enqueue_buffer(Buffer buffer, block_mode blocking)
    {
        //get name from buffer
        if (buffer.m_name.empty())
            log::warn("Encountered unnamed buffer! binding to a Kernel-location will fail!");
        return set_and_enqueue_buffer(buffer, buffer.m_name, blocking);
    }


    Kernel& Kernel::set_and_enqueue_buffer(Buffer buffer, const std::string& name, block_mode blocking)
    {
        //translate name to index
        param_find([this, b = std::forward<Buffer>(buffer), blocking](cl_uint index)
        {
            set_and_enqueue_buffer(b, index, blocking);
        }, name);
        return *this;
    }

    Kernel& Kernel::set_and_enqueue_buffer(Buffer buffer, cl_uint index, block_mode blocking)
    {
        //set and ... enqueue_buffer
        //nothing fun to see here
        enqueue_buffer(std::forward<Buffer>(buffer), index, blocking);
        return set_buffer(std::forward<Buffer>(buffer), index);
    }


    Kernel& Kernel::dispatch()
    {
        //enqueue the Kernel in the command queue
        cl_int ret = clEnqueueNDRangeKernel(
            m_queue,
            m_func,
            1,
            nullptr,
            &m_global_size,
            &m_local_size,
            0,
            nullptr,
            nullptr
        );

        //check if the enqueue was successful
        if (ret != CL_SUCCESS)
        {
            log::error("clEnqueueNDRangeKernel failed: {}", ret);
        }
        return *this;
    }

    void Kernel::finish() const
    {
        //execute all commands in the queue
        clFlush(m_queue);

        //probably not necessarily necessary
        clFinish(m_queue);

        //neccessary 
        clReleaseCommandQueue(m_queue);
    }


    Kernel::Kernel(Program* program, cl_kernel kernel) :
        m_default_mode(buffer_type::READ_BUFFER),
        m_prog(program),
        m_func(kernel),
        m_global_size(0),
        m_local_size(64)
    {
        cl_int ret;
        m_queue = program->make_cq();
    }

    Kernel& Kernel::local(size_t s)
    {
        m_local_size = s;
        return *this;
    }

    Kernel& Kernel::global(size_t s)
    {
        m_global_size = s;
        return *this;
    }
}
