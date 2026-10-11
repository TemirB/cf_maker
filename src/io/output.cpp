#include "io/output.h"

#include <filesystem>
#include <memory>
#include <stdexcept>
#include <string>

#include <TDirectory.h>
#include <TFile.h>
#include <TObject.h>

namespace
{
class CheckedOutputFile final : public TFile
{
  public:
    explicit CheckedOutputFile(const std::string& path) : TFile(path.c_str(), "RECREATE") {}

    [[nodiscard]] bool close_failed() const
    {
        return close_failed_;
    }

  private:
    Int_t SysClose(Int_t descriptor) override
    {
        const Int_t result = TFile::SysClose(descriptor);
        if (result != 0) {
            close_failed_ = true;
        }
        return result;
    }

    bool close_failed_ = false;
};

void check_output_file(const TFile& file)
{
    if (file.IsZombie() || !file.IsOpen() || !file.IsWritable() ||
        file.TestBit(TFile::kWriteError)) {
        throw std::runtime_error("output file is not writable: " + std::string(file.GetName()));
    }
}
} // namespace

std::unique_ptr<TFile> create_output_file(const std::string& path)
{
    if (std::filesystem::is_directory(path)) {
        throw std::runtime_error("output file path is a directory: " + path);
    }
    auto file = std::make_unique<CheckedOutputFile>(path);
    check_output_file(*file);
    return file;
}

void write_output_object(TFile& file, TObject& object, const char* name, int options)
{
    check_output_file(file);
    TDirectory::TContext context(&file);
    if (object.Write(name, options) <= 0 || file.TestBit(TFile::kWriteError)) {
        throw std::runtime_error("cannot write object " +
                                 std::string(name ? name : object.GetName()) + " to " +
                                 file.GetName());
    }
}

void finish_output_file(TFile& file)
{
    check_output_file(file);
    // Explicitly written, detached objects need not contribute to TFile::Write's
    // byte count, so zero is valid when the directory itself has no objects.
    if (file.Write() < 0 || file.TestBit(TFile::kWriteError)) {
        throw std::runtime_error("cannot finalize output file: " + std::string(file.GetName()));
    }
    file.Flush();
    if (file.TestBit(TFile::kWriteError)) {
        throw std::runtime_error("cannot flush output file: " + std::string(file.GetName()));
    }
    file.Close();
    const auto* checked = dynamic_cast<const CheckedOutputFile*>(&file);
    if (file.IsOpen() || file.TestBit(TFile::kWriteError) || (checked && checked->close_failed())) {
        throw std::runtime_error("cannot close output file: " + std::string(file.GetName()));
    }
}
