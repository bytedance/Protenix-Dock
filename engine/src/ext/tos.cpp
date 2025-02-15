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

#ifdef ENABLE_VETOS_ACCESS
#include "bytedock/ext/tos.h"

#include <cstdlib>
#include <sstream>
#include <thread>

#include <boost/filesystem/operations.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <boost/property_tree/json_parser.hpp>

#include "bytedock/lib/error.h"

namespace pt = boost::property_tree;
namespace fs = boost::filesystem;

namespace bytedock {

toutiao_object_storage& toutiao_object_storage::singleton() {
    static toutiao_object_storage instance;
    return instance;
}

inline std::string get_mc_or_rclone_config_file() {
    const char* value = std::getenv("HOME");
    if (!value) {
        throw failed_conf_error("ENV variable \"HOME\" is not defined!");
    }
    fs::path home(value);
    auto file = home / ".mc" / "config.json";
    if (fs::exists(file)) {
        std::cerr << "Loading credentials from config file of MinIO client ..."
                  << std::endl;
        return file.string();
    }
    file = home / ".config" / "rclone" / "rclone.conf";
    if (fs::exists(file)) {
        std::cerr << "Loading credentials from config file of RClone tool ..."
                  << std::endl;
        return file.string();
    }
    throw failed_conf_error(
        "Config file of neither MinIO client nor RClone tool could be found!"
    );
}

std::string get_region(const std::string& endpoint) {
    static const char* regions[] = {
        "cn-beijing", "cn-shanghai", "cn-guangzhou"
    };
    constexpr size_t nallowed = sizeof(regions) / sizeof(regions[0]);
    for (size_t i = 0; i < nallowed; i++) {
        auto& candidate = regions[i];
        if (endpoint.find(candidate) != std::string::npos) return candidate;
    }
    throw failed_conf_error("No valid region found in TOS URL: " + endpoint);
}

inline std::string s3_to_standard(const std::string& endpoint) {
    size_t start = endpoint.find("-s3");
    if (start == std::string::npos) return endpoint;
    return endpoint.substr(0, start) + endpoint.substr(start+3);
}

std::unique_ptr<VolcengineTos::TosClientV2> create_ve_tos_client(
    std::string& config_file
) {
    pt::ptree properties;
    std::string url, ak, sk;
    if (config_file.compare(config_file.size() - 5, 5, ".json") == 0) {
        pt::read_json(config_file, properties);
        url = properties.get<std::string>("aliases.tos.url");
        ak = properties.get<std::string>("aliases.tos.accessKey");
        sk = properties.get<std::string>("aliases.tos.secretKey");
    } else {
        pt::read_ini(config_file, properties);
        url = properties.get<std::string>("tos.endpoint");
        ak = properties.get<std::string>("tos.access_key_id");
        sk = properties.get<std::string>("tos.secret_access_key");
    }
    std::string region = get_region(url);
    VolcengineTos::ClientConfig conf;
    conf.endPoint = s3_to_standard(url);
    conf.enableCRC = true;
    return std::make_unique<VolcengineTos::TosClientV2>(region, ak, sk, conf);
}

toutiao_object_storage::toutiao_object_storage() {
    VolcengineTos::InitializeClient();
    auto config_file = get_mc_or_rclone_config_file();
    client_ = create_ve_tos_client(config_file);
}

toutiao_object_storage::~toutiao_object_storage() {
    VolcengineTos::CloseClient();
}

int64_t toutiao_object_storage::get_file_length(
    const std::string& bucket, const std::string& key
) {
    VolcengineTos::HeadObjectV2Input input(bucket, key);
    auto output = client_->headObject(input);
    if (output.isSuccess()) {
        return output.result().getContentLength();
    } else {
        std::ostringstream oss;
        oss << "[HeadObject] tos:" << bucket << "/" << key << std::endl
                << output.error().String();
        throw file_system_error(oss.str());
    }
}

std::shared_ptr<std::istream> toutiao_object_storage::open_for_read(
    const std::string& bucket, const std::string& key
) {
    VolcengineTos::GetObjectV2Input input(bucket, key);
    auto output = client_->getObject(input);
    if (output.isSuccess()) {
        auto derived = output.result().getContent();
        return std::static_pointer_cast<std::istream>(derived);
    } else {
        std::ostringstream oss;
        oss << "[GetObject] tos:" << bucket << "/" << key << std::endl
            << output.error().String();
        throw file_system_error(oss.str());
    }
}

std::shared_ptr<std::ostream> toutiao_object_storage::open_for_write(
    const std::string& bucket, const std::string& key
) {
    auto derived = std::shared_ptr<std::stringstream>(
        new std::stringstream, [=](std::stringstream* p) {
        /**
         * Do not delete `p` in flusher because its consumer `input` is
         * destructed later.
         */
        std::shared_ptr<std::stringstream> owned(p);
        VolcengineTos::PutObjectV2Input input(bucket, key, owned);
        auto output = this->client_->putObject(input);
        if (!output.isSuccess()) {
            std::ostringstream oss;
            oss << "[PutObject] tos:" << bucket << "/" << key << std::endl
                << output.error().String();
            /**
             * It is OK to throw exception in destructor because our purpose is
             * to commit rather than clean.
             */
            throw file_system_error(oss.str());
        }
    });
    return std::static_pointer_cast<std::ostream>(derived);
}

class otosbuf : public std::stringbuf {
public:
    otosbuf(
        const VolcengineTos::TosClientV2* client,
        const std::string& bucket, const std::string& key,
        int64_t offset, uint64_t hash
    ): std::stringbuf(std::ios_base::out),
       client_(client), bucket_(bucket), key_(key),
       offset_(offset), hash_(hash) {};

    /**
     * Since fields of `offset` and `hash` are calculated before instantiation,
     * it is only valid for the first appending request. Thus, `otosbuf` can
     * not be copy-constructible which prevents any copy appending text with
     * wrong values of offset & hash.
     */
    DISABLE_COPY_AND_ASSIGN(otosbuf);

    virtual ~otosbuf() {
        sync();
    }

protected:
    virtual int sync() {
        auto buffer = str();
        if (buffer.size() > 0) {
            auto stream = std::make_shared<std::stringstream>(
                std::move(buffer)
            );
            VolcengineTos::AppendObjectV2Input input(
                bucket_, key_, stream, offset_
            );
            input.setPreHashCrc64Ecma(hash_);
            auto output = client_->appendObject(input);
            if (!output.isSuccess()) {
                std::ostringstream oss;
                oss << "[AppendObject] tos:" << bucket_ << "/" << key_
                    << std::endl << output.error().String();
                throw file_system_error(oss.str());
            }
            offset_ = output.result().getNextAppendOffset();
            hash_ = output.result().getHashCrc64ecma();
            str("");  // Clear buffer after commitment to TOS
        }
        return 0;  // Success
    }

private:
    const VolcengineTos::TosClientV2* client_;
    const std::string bucket_;
    const std::string key_;
    int64_t offset_;
    uint64_t hash_;
};

class otstream : public std::ostream {
public:
    otstream(
        const VolcengineTos::TosClientV2* client,
        const std::string& bucket, const std::string& key,
        int64_t offset, uint64_t hash
    ): sb_(client, bucket, key, offset, hash) {
        init(&sb_);
    }

private:
    otosbuf sb_;  // This field prevents `otstream` being copy-constructible
};

std::unique_ptr<std::ostream> toutiao_object_storage::open_for_append(
    const std::string& bucket, const std::string& key
) {
    int64_t offset = 0;
    uint64_t hash = 0;
    VolcengineTos::HeadObjectV2Input input(bucket, key);
    auto output = client_->headObject(input);
    if (output.isSuccess()) {
        auto& result = output.result();
        if (result.getObjectType() != "Appendable") {
            std::ostringstream oss;
            oss << "Existed object [tos:" << bucket << "/" << key
                << "] is not appendable!";
            throw file_system_error(oss.str());
        }
        offset = result.getContentLength();
        hash = result.getHashCrc64Ecma();
    } else if (output.error().getStatusCode() != 404) {
        std::ostringstream oss;
        oss << "[HeadObject] tos:" << bucket << "/" << key << std::endl
                << output.error().String();
        throw file_system_error(oss.str());
    }
    /**
     * Since `toutiao_object_storage` is a static field, it is destructed
     * after `otstream` on heap is destroyed. Thus, passing a raw pointer
     * managed by `unique_ptr` is safe.
     */
    auto raw = new otstream(client_.get(), bucket, key, offset, hash);
    return std::unique_ptr<std::ostream>(raw);
}

inline void adjust_upload_parallelism(
    const std::string& local_path, size_t& user_value
) {
    if (user_value == 0) {
        uintmax_t nthreads = std::thread::hardware_concurrency();
        uintmax_t nparts = fs::file_size(local_path) / (4*1024*1024);
        user_value = nparts > 0 ? std::min(nthreads, nparts) : nthreads;
    }
}

void toutiao_object_storage::upload_file(
    const std::string& bucket, const std::string& key,
    const std::string& local_path, size_t num_threads
) {
    VolcengineTos::UploadFileV2Input input(bucket, key);
    adjust_upload_parallelism(local_path, num_threads);
    input.setTaskNum(num_threads);
    input.setFilePath(local_path);
    auto output = client_->uploadFile(input);
    if (!output.isSuccess()) {
        std::ostringstream oss;
        oss << "[UploadFile] tos:" << bucket << "/" << key << std::endl
            << output.error().String();
        throw file_system_error(oss.str());
    }
}

bool toutiao_object_storage::is_file(const std::string& bucket,
                                     const std::string& key) {
    VolcengineTos::HeadObjectV2Input input(bucket, key);
    auto output = client_->headObject(input);
    if (output.isSuccess()) {
        return true;
    } else if (output.error().getStatusCode() == 404) {
        return false;
    } else {
        std::ostringstream oss;
        oss << "[HeadObject] tos:" << bucket << "/" << key << std::endl
            << output.error().String();
        throw file_system_error(oss.str());
    }
}

}
#endif
