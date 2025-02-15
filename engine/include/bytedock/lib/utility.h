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

#define DISABLE_COPY_AND_ASSIGN(className) \
    className(const className&) = delete; \
    className& operator=(const className&) = delete

#define CHECK_PARAMETER_POSITIVE(x) ((x) > 1e-5)
#define CHECK_PARAMETER_NON_POSITIVE(x) ((x) < 1e-5)
#define CHECK_PARAMETER_NON_ZERO(x) ((x) > 1e-5 || (x) < -1e-5)
