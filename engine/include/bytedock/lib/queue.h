/*
* Copyright (C) 2025 ByteDance and/or its affiliates
*
* This program is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#pragma once

#include <condition_variable>
#include <cassert>
#include <mutex>
#include <queue>

#include "bytedock/lib/utility.h"

namespace bytedock {

template <typename>
struct is_pair: std::false_type {};

template <class A, class B>
struct is_pair<std::pair<A, B> >: std::true_type {};

/**
 * This blocking queue is designed for the scenario where:
 * - there are multiple producers
 * - there are multiple consumers
 * - message number is finite
 * - queue size is limited for memory size
 *
 * Demo code snippet for producer is shared as below:
 * ```c++
 *   size_t capacity = 8;
 *   std:string eoq = "ANY UNIQUE VALUE";
 *   size_t nreaders = 4;
 *   block_queue<std:string> queue(capacity, eoq, nreaders);
 *   for (size_t i = 0; i < 100) {
 *       queue.push(std::to_string(i));
 *   }
 *   queue.close();  // Notify consumers to quit
 * ``` 
 *
 * Demo code snippet for consumer is shared as below:
 * ```c++
 *   std::string path = std::move(queue.pop());
 *   while (!queue.is_eoq(path)) {
 *     path = std::move(queue.pop());
 *   }
 *   queue.close();  // Must give back the EOQ flag
 * ```
 */
template<class T>
class blocking_queue {
public:
    blocking_queue(const size_t capacity,
                   const T& eoq,
                   const int nproducers = 1
                   ) : capacity_(capacity), eoq_(eoq), nproducers_(nproducers) {
        assert(capacity_ > 0);
        assert(nproducers_ > 0);
    }

    DISABLE_COPY_AND_ASSIGN(blocking_queue);

    void push(T& data) {
        std::unique_lock<std::mutex> lock(mtx_);
        while (queue_.size() >= capacity_) {
            full_.wait(lock);
        }
        assert(queue_.size() < capacity_);
        queue_.push(data);
        empty_.notify_all();
    }

    void push(T&& data) {
        std::unique_lock<std::mutex> lock(mtx_);
        while (queue_.size() >= capacity_) {
            full_.wait(lock);
        }
        assert(queue_.size() < capacity_);
        queue_.push(std::move(data));
        empty_.notify_all();
    }

    T pop() {
        std::unique_lock<std::mutex> lock(mtx_);
        while(queue_.empty()) {
            empty_.wait(lock);
        }
        assert(!queue_.empty());
        T out = std::move(queue_.front());
        queue_.pop();
        full_.notify_all();
        return out;
    }

    void close() {
        std::lock_guard<std::mutex> lock(mtx_);
        if (nproducers_ > 0) {
            --nproducers_;
        }
        if (nproducers_ < 1) {
            queue_.push(eoq_);
            empty_.notify_all();
        }
    }

    bool is_eoq(const T& item) const {
        if constexpr (is_pair<T>::value) {
            return std::get<0>(eoq_) == std::get<0>(item);
        } else {
            return eoq_ == item;
        }
    }

private:
    mutable std::mutex mtx_;
    std::condition_variable full_;
    std::condition_variable empty_;
    std::queue<T> queue_;
    const size_t capacity_;
    const T eoq_;  // End of queue
    int nproducers_;
};

}
