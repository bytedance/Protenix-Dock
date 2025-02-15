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

#include "bytedock/core/source.h"

#include "bytedock/ext/logging.h"
#include "bytedock/ext/pfile.h"

#include <fstream>

namespace bytedock {

void index_file_source::fill(blocking_queue<std::string>& queue) {
    auto infile = open_for_read(path_);
    std::string line;
    size_t line_no = 0;
    while (std::getline(*infile, line)) {
        line_no += 1;
        if (queue.is_eoq(line)) {
            LOG_WARNING << "Skip illegal line#" << line_no
                        << " with empty text.";
        } else if (!is_file(line)) {
            LOG_WARNING << "Skip illegal line#" << line_no
                        << " with ghost file: " << line;
        } else {
            queue.push(std::move(line));
        }
    }
    queue.close();
}

}
