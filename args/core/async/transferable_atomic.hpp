#pragma once
#include <atomic>
#include <core/async/readonly_rw_spinlock.hpp>

/** @file transferable_atomic.hpp
 */

namespace args::core::async
{
	/**@class transferable_atomic
	 * @brief Copyable wrapper for std::atomic.
	 */
	template <typename T>
	struct transferable_atomic
	{
	private:
		std::atomic<T> m_atomic;
		mutable readonly_rw_spinlock m_lock;
	public:
		transferable_atomic() noexcept = default;

		constexpr transferable_atomic(T val) noexcept : m_atomic(val), m_lock() {}

		transferable_atomic(const std::atomic<T>& other) : m_atomic(other.load(std::memory_order_acquire)), m_lock() {}

		transferable_atomic(const transferable_atomic<T>& other) : m_atomic(other->load(std::memory_order_acquire)), m_lock() {}

		readonly_rw_spinlock& get_lock()
		{
			return m_lock;
		}

		transferable_atomic<T>& operator=(const transferable_atomic<T>& other)
		{
			m_atomic.store(other->load(std::memory_order_acquire), std::memory_order_release);
			return *this;
		}

		transferable_atomic<T>& operator=(const std::atomic<T>& other)
		{
			m_atomic.store(other.load(std::memory_order_acquire), std::memory_order_release);
			return *this;
		}

		transferable_atomic<T>& copy(const transferable_atomic<T>& other, std::memory_order loadOrder = std::memory_order_acquire, std::memory_order storeOrder = std::memory_order_release)
		{
			m_atomic.store(other->load(loadOrder), storeOrder);
			return *this;
		}

		transferable_atomic<T>& move(const transferable_atomic<T>& other, std::memory_order loadOrder = std::memory_order_acquire, std::memory_order storeOrder = std::memory_order_release)
		{
			m_atomic.store(other->load(loadOrder), storeOrder);
			other->store(T(), storeOrder);
			return *this;
		}

		std::atomic<T>& get()
		{
			return m_atomic;
		}

		const std::atomic<T>& get() const
		{
			return m_atomic;
		}

		std::atomic<T>* operator->()
		{
			return &m_atomic;
		}

		const std::atomic<T>* operator->() const
		{
			return &m_atomic;
		}
	};
}