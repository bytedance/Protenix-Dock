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

#include "bytedock/ext/logging.h"

#include <iostream>

#include <boost/log/expressions.hpp>
#include <boost/log/attributes.hpp>
#include <boost/log/support/date_time.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>

#include "bytedock/ext/pfile.h"
#include "bytedock/lib/error.h"

namespace bytedock {

namespace attrs = boost::log::attributes;
namespace expr = boost::log::expressions;

boost::shared_ptr<sink_t> enable_global_logging(const std::string& path, int verbosity) {
    auto backend = boost::make_shared<backend_t>();
    if (path == "-") {
        boost::shared_ptr<std::ostream> cerr(&std::clog, [](auto* p){});
        backend->add_stream(cerr);
    } else {
        boost::shared_ptr<std::ostream> flog(open_for_append(path).release());
        if (!flog->good()) {
            throw file_system_error("Failed to open a text file for logging!");
        }
        backend->add_stream(flog);
    }
    backend->auto_flush(true);

    auto sink = boost::make_shared<sink_t>(backend);
    sink->set_formatter(expr::stream
        << "[" << expr::format_date_time<boost::posix_time::ptime>("TimeStamp", "%Y-%m-%d %H:%M:%S")
        << " - " << std::setw(3) << std::setfill('0') << expr::attr<size_t>("ThreadID")
        << " - " << expr::attr<trivial::severity_level>("Severity")
        << "] " << expr::smessage);

    auto core = boost::log::core::get();
    core->add_global_attribute("TimeStamp", attrs::local_clock());
    // Atribute "ThreadID" must be told by macro LOG_TELL_THREAD_ID()
    core->reset_filter();
    core->set_filter(trivial::severity > static_cast<trivial::severity_level>(2 - verbosity));
    core->add_sink(sink);
    return sink;
}

void disable_sink(boost::shared_ptr<sink_t > sink, bool deregister) {
    if (deregister) boost::log::core::get()->remove_sink(sink);
    sink->stop();
    sink->flush();
    sink.reset();
}

}
