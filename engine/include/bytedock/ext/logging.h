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

#include <boost/log/common.hpp>
#include <boost/log/sinks.hpp>
#include <boost/log/trivial.hpp>

#include "bytedock/ext/counter.h"

namespace bytedock {

namespace sinks = boost::log::sinks;
namespace src = boost::log::sources;
namespace trivial = boost::log::trivial;

using backend_t = sinks::text_ostream_backend;
using sink_t = sinks::asynchronous_sink<backend_t>;

BOOST_LOG_INLINE_GLOBAL_LOGGER_DEFAULT(global_logger, src::severity_logger_mt<trivial::severity_level>)

#define LOG_TRACE BOOST_LOG_SEV(global_logger::get(), trivial::trace)
#define LOG_DEBUG BOOST_LOG_SEV(global_logger::get(), trivial::debug)
#define LOG_INFO BOOST_LOG_SEV(global_logger::get(), trivial::info)
#define LOG_WARNING BOOST_LOG_SEV(global_logger::get(), trivial::warning)
#define LOG_ERROR BOOST_LOG_SEV(global_logger::get(), trivial::error)
#define LOG_FATAL BOOST_LOG_SEV(global_logger::get(), trivial::fatal)

#define LOG_TELL_COUNTER(category) BOOST_LOG_SCOPED_THREAD_TAG((category), get_sequential_id((category)));
#define LOG_TELL_THREAD_ID(...) LOG_TELL_COUNTER("ThreadID")

boost::shared_ptr<sink_t> enable_global_logging(const std::string& path, int verbosity);

void disable_sink(boost::shared_ptr<sink_t> sink, bool deregister = true);

}
