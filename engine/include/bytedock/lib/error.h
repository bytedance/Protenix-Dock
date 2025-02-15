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

#include <exception>
#include <string>

namespace bytedock {

class message_based_error : public std::exception {
public:
    explicit message_based_error(const std::string& message) : message_(message) {}
    explicit message_based_error(std::string&& message) : message_(std::move(message)) {}

    const char* what() const throw () {
        return message_.c_str();
    }

protected:
    const std::string message_;
};

class failed_conf_error : public message_based_error {
public:
    explicit failed_conf_error(const std::string& message) : message_based_error(message) {}
    explicit failed_conf_error(std::string&& message) : message_based_error(std::move(message)) {}
};

class file_system_error : public message_based_error {
public:
    explicit file_system_error(const std::string& message) : message_based_error(message) {}
    explicit file_system_error(std::string&& message) : message_based_error(std::move(message)) {}
};

}
