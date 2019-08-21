#pragma once
#include "ICLmanager.h"
#include <memory>
#include <vector>
#include <optional>

namespace clde{
template <class T>
class CLDataStorage
{
private:
    std::shared_ptr<ICLmanager> m_context;
    unsigned int m_size;
    std::optional<cl_mem> m_data;

public:
    CLDataStorage(const std::shared_ptr<ICLmanager>& context):m_context(context), m_size(0)
    {

    }

    CLDataStorage(const unsigned int& size, const std::shared_ptr<ICLmanager>& context):m_context(context), m_size(size)
    {
        cl_int err = 0;
        m_data = clCreateBuffer(m_context->context(), CL_MEM_READ_WRITE, m_size * sizeof(T), nullptr, &err);
        if( err < 0)
            throw std::runtime_error("OpenCL: unable to create memory buffer");
    }

    CLDataStorage(const std::vector<T>& in, const std::shared_ptr<ICLmanager>& context):m_context(context), m_size(in.size())
    {
        cl_int err = 0;
        cl_event copy_event;
        m_data = clCreateBuffer(m_context->context(), CL_MEM_READ_WRITE , m_size * sizeof(T), nullptr, &err);
        if( err < 0)
            throw std::runtime_error("OpenCL: unable to create memory buffer");
        err =  clEnqueueWriteBuffer(m_context->command_queue(), *m_data, CL_FALSE, 0, m_size * sizeof(T), const_cast<T*>(in.data()), 0, nullptr, &copy_event );
        if( err < 0)
            throw std::runtime_error("OpenCL: unable to copy memory buffer");
        clWaitForEvents(1, &copy_event);
    }

    CLDataStorage(const CLDataStorage& in):m_context(in.m_context), m_size(in.m_size)
    {
        if(in.m_data){
            cl_int err = 0;
            cl_event copy_event;
            m_data = clCreateBuffer(m_context->context(), CL_MEM_READ_WRITE, m_size * sizeof(T), nullptr, &err);
            if( err < 0)
                throw std::runtime_error("OpenCL: unable to create memory buffer");
            err = clEnqueueCopyBuffer(m_context->command_queue(), *(in.m_data), *m_data, 0, 0, m_size * sizeof(T), 0, nullptr, &copy_event);
            if( err < 0)
                throw std::runtime_error("OpenCL: unable to copy memory buffer");
            clWaitForEvents(1, &copy_event);
        }
    }

    CLDataStorage(CLDataStorage&& in):m_context(in.m_context), m_size(in.m_size)
    {
        m_data = *(in.m_data);
        in.m_data.reset();
    }

    ~CLDataStorage()
    {
        if(m_data)
            clReleaseMemObject(*m_data);
    }

    CLDataStorage& operator=(const std::vector<T>& in)
    {
        cl_int err = 0;
        cl_event copy_event;
        if( in.size() != m_size)
        {
            if(m_data)
                clReleaseMemObject(*m_data);
            m_size = in.size();

            m_data = clCreateBuffer(m_context->context(), CL_MEM_READ_WRITE, m_size * sizeof(T), nullptr, &err);
            if( err < 0)
                throw std::runtime_error("OpenCL: unable to create memory buffer");

        }
        err =  clEnqueueWriteBuffer(m_context->command_queue(), *m_data, CL_FALSE, 0, m_size * sizeof(T), const_cast<T*>(in.data()), 0, nullptr, &copy_event );
        if( err < 0)
            throw std::runtime_error("OpenCL: unable to copy memory buffer");
        clWaitForEvents(1, &copy_event);
        return *this;
    }
    CLDataStorage& operator=(const CLDataStorage& in)
    {
        if(in.m_data)
        {
            cl_event copy_event;
            cl_int err = 0;
            if(m_size != in.m_size)
            {
                if(m_data)
                    clReleaseMemObject(*m_data);
                m_size = in.m_size;
                m_data = clCreateBuffer(m_context->context(), CL_MEM_READ_WRITE, m_size * sizeof(T), nullptr, &err);
                if( err < 0)
                    throw std::runtime_error("OpenCL: unable to create memory buffer");
            }
            err = clEnqueueCopyBuffer(m_context->command_queue(), *(in.m_data), *m_data, 0, 0, m_size * sizeof(T), 0, nullptr, &copy_event);
            if( err < 0)
                throw std::runtime_error("OpenCL: unable to copy memory buffer");
            clWaitForEvents(1, &copy_event);

        }
        else if(m_data){
            clReleaseMemObject(*m_data);
            m_data.reset();
            m_size = 0;
        }
        return *this;
    }

    const cl_mem& data() const
    {
        if(!m_data)
            throw std::runtime_error("trying to access non existing memory object");
        return *m_data;
    }

    std::vector<T> read() const
    {
        std::vector<T> result;
        if (m_data && m_size > 0)
        {
            cl_int err = 0;
            result.resize(m_size);
            err = clEnqueueReadBuffer(m_context->command_queue(), *m_data, CL_TRUE, 0, m_size * sizeof(T), result.data(), 0, nullptr, nullptr);
            if( err < 0 )
                throw std::runtime_error("OpenCL: unable to read memory buffer");
        }
        return result;
    }

    const unsigned int& size() const
    {
        return m_size;
    }
};
};
