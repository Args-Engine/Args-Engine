#include <core/scheduling/processchain.hpp>
#include <core/scheduling/process.hpp>
#include <core/scheduling/scheduler.hpp>
#include <core/common/exception.hpp>
#include <thread>

#include <core/logging/logging.hpp>

namespace args::core::scheduling
{
    void ProcessChain::threadedRun(ProcessChain* chain)
    {
        log::info("Chain started.");
        try
        {
            while (!chain->m_exit->load(std::memory_order_acquire)) // Check for exit flag.
            {
                chain->runInCurrentThread(); // Execute all processes.

                if (chain->m_scheduler->syncRequested()) // Sync if requested.
                    chain->m_scheduler->waitForProcessSync();

                if (chain->m_low_power)
                    std::this_thread::sleep_for(std::chrono::milliseconds(1));
            }

            chain->m_scheduler->reportExit(chain->m_threadId); // Mark Exit.
        }
        catch (const args::core::exception& e)
        {
            chain->m_scheduler->reportExitWithError(chain->m_name, std::this_thread::get_id(), e); // Mark error.
        }
        catch (const std::exception& e)
        {
            chain->m_scheduler->reportExitWithError(chain->m_name, std::this_thread::get_id(), e); // Mark error.
        }
    }

    bool ProcessChain::run(bool low_power)
    {
        m_low_power = low_power;
        m_exit->store(false, std::memory_order_release);
        std::thread::id threadId = m_scheduler->getChainThreadId(m_nameHash);
        if (threadId != std::thread::id())
        {
            m_threadId = threadId;

            m_scheduler->sendCommand(threadId, [](void* param) { ProcessChain::threadedRun(reinterpret_cast<ProcessChain*>(param)); }, this);

            return true;
        }
        return false;
    }

    void ProcessChain::exit()
    {
        m_exit->store(true, std::memory_order_release);
    }

    void ProcessChain::runInCurrentThread()
    {
        hashed_sparse_set<id_type> finishedProcesses;
        async::readonly_guard guard(m_processesLock); // Hooking more processes whilst executing isn't allowed.
        do
        {
            for (auto [id, process] : m_processes)
                if (!finishedProcesses.contains(id))
                    if (process->execute(m_scheduler->getTimeScale())) // If the process wasn't finished then execute it and check if it's finished now.
                        finishedProcesses.insert(id);

        } while (finishedProcesses.size() != m_processes.size());
    }

    void ProcessChain::addProcess(Process* process)
    {
        async::readwrite_guard guard(m_processesLock);
        if (m_processes.insert(process->id(), process).second)
            process->m_hooks.insert(m_nameHash);
    }

    void ProcessChain::removeProcess(Process* process)
    {
        async::readwrite_guard guard(m_processesLock);
        if (m_processes.erase(process->id()))
            process->m_hooks.erase(m_nameHash);
    }
}
