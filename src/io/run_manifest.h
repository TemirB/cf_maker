#pragma once

#include <memory>

struct Config;

// Construct after preflight and write_run_config. A fatal exception leaves the
// manifest in running state, so consumers cannot reuse a previous successful run.
class RunManifest
{
  public:
    explicit RunManifest(const Config& cfg);
    ~RunManifest();
    RunManifest(const RunManifest&) = delete;
    RunManifest& operator=(const RunManifest&) = delete;
    RunManifest(RunManifest&&) = delete;
    RunManifest& operator=(RunManifest&&) = delete;

    // Exit codes preserve the CLI contract: 0 is complete, 2 has unusable fits.
    void finish(const Config& cfg, int exit_code);

  private:
    struct State;
    std::unique_ptr<State> state_;
};
