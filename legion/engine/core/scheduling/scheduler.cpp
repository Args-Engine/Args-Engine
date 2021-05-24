#include <core/scheduling/scheduler.hpp>

namespace legion::core::scheduling
{
    constexpr size_type reserved_threads = 1; // this, OS

    sparse_map<id_type, ProcessChain> Scheduler::m_processChains;

    multicast_delegate<Scheduler::frame_callback_type> Scheduler::m_onFrameStart;
    multicast_delegate<Scheduler::frame_callback_type> Scheduler::m_onFrameEnd;

    const size_type Scheduler::m_maxThreadCount = reserved_threads >= std::thread::hardware_concurrency() ? 0 : std::thread::hardware_concurrency() - reserved_threads;
    size_type Scheduler::m_availableThreads = m_maxThreadCount;

    Scheduler::per_thread_map<std::thread> Scheduler::m_threads;

    Scheduler::per_thread_map<async::rw_lock_pair<async::runnables_queue>> Scheduler::m_commands;
    async::rw_lock_pair<async::job_queue> Scheduler::m_jobs;
    size_type Scheduler::m_jobPoolSize = 0;

    std::atomic<bool> Scheduler::m_exit = { false };
    std::atomic<bool> Scheduler::m_start = { false };
    int Scheduler::m_exitCode = 0;

    std::atomic<float> Scheduler::m_pollTime = { 0.1f };

    void Scheduler::threadMain(bool lowPower, std::string name)
    {
        log::info("Thread {} assigned.", std::this_thread::get_id());
        async::set_thread_name(name.c_str());

        OPTICK_THREAD(name.c_str());

        while (!m_start.load(std::memory_order_relaxed))
            std::this_thread::yield();

        time::timer clock;
        time::span timeBuffer;
        time::span sleepTime;

        while (!m_exit.load(std::memory_order_relaxed))
        {
            bool noWork = true;

            {
                auto& [lock, commandQueue] = m_commands.at(std::this_thread::get_id());

                async::readonly_guard guard(lock);

                if (!commandQueue.empty())
                {
                    noWork = false;
                    commandQueue.front()->execute();
                    commandQueue.pop();
                }
            }

            {
                auto& [lock, jobQueue] = m_jobs;
                async::readonly_guard guard(lock);
                if (!jobQueue.empty())
                {
                    noWork = false;

                    auto jobPoolPtr = jobQueue.front();
                    if (jobPoolPtr->is_done())
                    {
                        async::readwrite_guard wguard(lock);
                        if (!jobQueue.empty() && jobQueue.front()->is_done())
                            jobQueue.pop();
                    }
                    else
                        jobPoolPtr->complete_job();
                }
            }

            if (noWork)
            {
                timeBuffer += clock.restart();

                if (lowPower || LEGION_CONFIGURATION == LEGION_DEBUG_VALUE)
                {
                    OPTICK_EVENT("Sleep");
                    std::this_thread::sleep_for(std::chrono::microseconds(1));
                }
                else if (timeBuffer >= sleepTime * m_pollTime.load(std::memory_order_relaxed))
                {
                    OPTICK_EVENT("Sleep");
                    timeBuffer -= sleepTime;

                    time::timer sleepTimer;
                    sleepTimer.start();
                    std::this_thread::sleep_for(std::chrono::nanoseconds(1));
                    sleepTime = sleepTimer.elapsed_time();
                }
                else
                    std::this_thread::yield();
            }
            else
            {
                timeBuffer = 0.f;
                clock.start();
            }
        }
    }

    void Scheduler::tryCompleteJobPool()
    {
        auto& [lock, jobQueue] = m_jobs;
        async::readwrite_guard wguard(lock);
        if (!jobQueue.empty())
        {
            if (jobQueue.front()->is_done())
            {
                jobQueue.pop();
            }
        }
    }

    void Scheduler::doTick(Clock::span_type deltaTime)
    {
        OPTICK_FRAME("Main thread");

        time::span dt{ deltaTime };
        m_onFrameStart(dt, time::span(Clock::elapsedSinceTickStart()));

        for (auto [_, chain] : m_processChains)
            chain.runInCurrentThread(dt);

        m_onFrameEnd(dt, time::span(Clock::elapsedSinceTickStart()));
    }

    pointer<std::thread> Scheduler::getThread(std::thread::id id)
    {
        return { &m_threads.at(id) };
    }

    int Scheduler::run(bool lowPower, size_type minThreads)
    {
        bool m_lowPower = lowPower;

        size_type m_minThreads = minThreads < 1 ? 1 : minThreads;

        if (m_maxThreadCount < m_minThreads)
            m_lowPower = true;

        if (m_availableThreads < m_minThreads)
            m_availableThreads = m_minThreads;

        pointer<std::thread> ptr;
        std::string name = "Worker ";

        {
            async::readwrite_guard guard(log::impl::thread_names_lock);
            while ((ptr = createThread(Scheduler::threadMain, m_lowPower, name + std::to_string(m_jobPoolSize))) != nullptr)
            {
                log::impl::thread_names[ptr->get_id()] = name + std::to_string(m_jobPoolSize++);
                m_commands.try_emplace(ptr->get_id());
            }
            log::impl::thread_names[std::this_thread::get_id()] = "Main thread";
        }

        m_start.store(true, std::memory_order_release);

        Clock::subscribeToTick(doTick);

        float prevPolltime = m_pollTime.load(std::memory_order_relaxed);

        size_type solidness = 1;

        while (!m_exit.load(std::memory_order_relaxed))
        {
            Clock::update();

            float deltaTime = static_cast<float>(Clock::lastTickDuration());
            static float avgDeltaTime = deltaTime;
            static size_type framecount = 1;
            static float totalTime = deltaTime;
            static float timeSpentAboveAvg = 0;

            uint fps = math::uround(1.f / deltaTime);

            totalTime += deltaTime;
            framecount++;
            avgDeltaTime = totalTime / static_cast<float>(framecount);

            uint avg = math::uround(1.f / avgDeltaTime);

            if (fps < avg)
            {
                m_pollTime.store(math::clamp(((prevPolltime * static_cast<float>(solidness)) + math::linearRand(-0.1f, 0.1f)) / static_cast<float>(solidness + 1), 0.f, 1.f), std::memory_order_relaxed);
                solidness = math::max<size_type>(solidness - 1, 0);

                timeSpentAboveAvg = 0.f;
            }
            else if (fps > avg)
            {
                timeSpentAboveAvg += deltaTime;

                solidness = math::min(solidness + 1, std::numeric_limits<size_type>::max());

                prevPolltime = m_pollTime.load(std::memory_order_relaxed);
            }

            if (timeSpentAboveAvg > 0.5f)
            {
                framecount = 1;
                totalTime = deltaTime;
                timeSpentAboveAvg = 0.f;
            }

            //log::debug("trend {: f}, poll percentage {:6f}, fps {}", m_pollTime.load(std::memory_order_relaxed) - prevPolltime, m_pollTime.load(std::memory_order_relaxed), avg);

            if (m_lowPower)
                std::this_thread::yield();
            else
                L_PAUSE_INSTRUCTION();
        }

        for (auto& [_, thread] : m_threads)
            if (thread.joinable())
                thread.join();

        return m_exitCode;
    }

    void Scheduler::exit(int exitCode)
    {
        m_exitCode = exitCode;
        m_exit.store(true, std::memory_order_release);
    }

    size_type Scheduler::jobPoolSize() noexcept
    {
        return m_jobPoolSize;
    }

    pointer<ProcessChain> Scheduler::createProcessChain(cstring name)
    {
        id_type id = nameHash(name);
        return { &m_processChains.emplace(id, name, id).first.value() };
    }

    pointer<ProcessChain> Scheduler::getChain(id_type id)
    {
        if (m_processChains.contains(id))
            return { &m_processChains.at(id) };
        return { nullptr };
    }

    pointer<ProcessChain> Scheduler::getChain(cstring name)
    {
        id_type id = nameHash(name);
        if (m_processChains.contains(id))
            return { &m_processChains.at(id) };
        return { nullptr };
    }

    void Scheduler::subscribeToChainStart(id_type chainId, const chain_callback_delegate& callback)
    {
        m_processChains.at(chainId).subscribeToChainStart(callback);
    }

    void Scheduler::subscribeToChainStart(cstring chainName, const chain_callback_delegate& callback)
    {
        id_type chainId = nameHash(chainName);
        m_processChains.at(chainId).subscribeToChainStart(callback);
    }

    void Scheduler::unsubscribeFromChainStart(id_type chainId, const chain_callback_delegate& callback)
    {
        m_processChains.at(chainId).unsubscribeFromChainStart(callback);
    }

    void Scheduler::unsubscribeFromChainStart(cstring chainName, const chain_callback_delegate& callback)
    {
        id_type chainId = nameHash(chainName);
        m_processChains.at(chainId).unsubscribeFromChainStart(callback);
    }

    void Scheduler::subscribeToChainEnd(id_type chainId, const chain_callback_delegate& callback)
    {
        m_processChains.at(chainId).subscribeToChainEnd(callback);
    }

    void Scheduler::subscribeToChainEnd(cstring chainName, const chain_callback_delegate& callback)
    {
        id_type chainId = nameHash(chainName);
        m_processChains.at(chainId).subscribeToChainEnd(callback);
    }

    void Scheduler::subscribeToFrameStart(const frame_callback_delegate& callback)
    {
        m_onFrameStart.push_back(callback);
    }

    void Scheduler::unsubscribeFromFrameStart(const frame_callback_delegate& callback)
    {
        m_onFrameStart.erase(callback);
    }

    void Scheduler::subscribeToFrameEnd(const frame_callback_delegate& callback)
    {
        m_onFrameEnd.push_back(callback);
    }

    void Scheduler::unsubscribeFromFrameEnd(const frame_callback_delegate& callback)
    {
        m_onFrameEnd.erase(callback);
    }

    void Scheduler::unsubscribeFromChainEnd(id_type chainId, const chain_callback_delegate& callback)
    {
        m_processChains.at(chainId).unsubscribeFromChainEnd(callback);
    }

    void Scheduler::unsubscribeFromChainEnd(cstring chainName, const chain_callback_delegate& callback)
    {
        id_type chainId = nameHash(chainName);
        m_processChains.at(chainId).unsubscribeFromChainEnd(callback);
    }

    bool Scheduler::hookProcess(cstring chainName, Process& process)
    {
        id_type chainId = nameHash(chainName);

        if (m_processChains.contains(chainId))
        {
            m_processChains[chainId].addProcess(process);
            return true;
        }

        return false;
    }

    bool Scheduler::hookProcess(cstring chainName, pointer<Process> process)
    {
        id_type chainId = nameHash(chainName);

        if (m_processChains.contains(chainId))
        {
            m_processChains[chainId].addProcess(process);
            return true;
        }

        return false;
    }

    bool Scheduler::unhookProcess(cstring chainName, Process& process)
    {
        id_type chainId = nameHash(chainName);

        if (m_processChains.contains(chainId))
            return m_processChains[chainId].removeProcess(process);

        return false;
    }

    bool Scheduler::unhookProcess(cstring chainName, pointer<Process> process)
    {
        id_type chainId = nameHash(chainName);

        if (m_processChains.contains(chainId))
            return m_processChains[chainId].removeProcess(process);

        return false;
    }

    bool Scheduler::unhookProcess(id_type chainId, pointer<Process> process)
    {
        if (m_processChains.contains(chainId))
            return m_processChains[chainId].removeProcess(process);

        return false;
    }

}
